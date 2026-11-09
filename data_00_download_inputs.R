# Download the initial input archive. Run from the repository root after setting
# FIGSHARE_INPUT_URL to the Figshare download link for the ZIP file.

download_inputs <- function() {
  if (!file.exists("data_00.a_Handle_AHOI_data.R") ||
      !dir.exists("functions")) {
    stop("Run this script from the repository root.")
  }

  url <- Sys.getenv("FIGSHARE_INPUT_URL", unset = "")
  if (!nzchar(url)) {
    stop("Set FIGSHARE_INPUT_URL to the Figshare ZIP download link first; ",
         "see README.md.")
  }
  if (!grepl("^(https|file)://", url, ignore.case = TRUE)) {
    stop("FIGSHARE_INPUT_URL must be an HTTPS or local file:// URL.")
  }

  archive <- tempfile(fileext = ".zip")
  staging <- tempfile(pattern = "wkfisheco-inputs-")
  if (!dir.create(staging)) stop("Could not create a staging directory.")
  on.exit(unlink(c(archive, staging), recursive = TRUE), add = TRUE)

  old_timeout <- getOption("timeout")
  options(timeout = max(1800, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  message("Downloading the input archive...")
  status <- utils::download.file(url, archive, mode = "wb", quiet = TRUE)
  if (status != 0L) stop("Download failed (status ", status, ").")

  # R's internal unzip cannot reliably handle individual files >= 4 GiB.
  # Use a ZIP64-capable Info-ZIP executable when available.
  unzip_program <- Sys.which("unzip")
  unzip_method <- if (nzchar(unzip_program)) unzip_program else "internal"
  if (unzip_method == "internal") {
    con <- file(archive, "rb")
    archive_size <- file.info(archive)$size
    seek(con, where = max(0, archive_size - 65577), origin = "start")
    tail_bytes <- readBin(con, "raw", n = min(archive_size, 65577))
    close(con)
    if (length(grepRaw(as.raw(c(0x50, 0x4b, 0x06, 0x07)),
                       tail_bytes, fixed = TRUE))) {
      stop("This ZIP uses ZIP64. Install a ZIP64-capable Info-ZIP unzip ",
           "executable before downloading these inputs.")
    }
  }
  zip_contents <- withCallingHandlers(
    utils::unzip(archive, list = TRUE, unzip = unzip_method),
    warning = function(w) stop("Cannot read ZIP archive: ", conditionMessage(w))
  )
  paths <- zip_contents$Name
  if (!length(paths) || anyNA(paths) || anyDuplicated(paths)) {
    stop("The ZIP archive is empty or has duplicate/invalid paths.")
  }
  if (any(grepl("\\\\|(^/)|(^|/)\\.\\.?(/|$)|^[A-Za-z]:|[?*\\[\"\r\n]", paths))) {
    stop("The ZIP archive contains unsafe paths.")
  }

  is_dir <- endsWith(paths, "/")
  files <- paths[!is_dir]
  file_sizes <- zip_contents$Length[!is_dir]
  if (!length(files) || anyNA(file_sizes) || any(file_sizes <= 0)) {
    stop("The ZIP archive contains no files or contains empty files.")
  }

  ahoi_dir <- "data/Env/AHOI/v23.05/"
  rect_dir <- "data/ICES_stat_rect/"
  rect_name <- "StatRec_map_Areas_Full_20170124"
  required <- c(
    paste0(rect_dir, rect_name, c(".shp", ".shx", ".dbf", ".prj")),
    "data/stock/Area_NSCOD3_table_6a_inshore_offshore.xlsx",
    "data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata",
    "data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp45_BottomTemp_raster_2006_2099.RData",
    "data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp85_BottomTemp_raster_2006_2099.RData"
  )
  missing <- setdiff(required, files)
  ahoi <- files[startsWith(files, ahoi_dir) &
                  !grepl("/", substring(files, nchar(ahoi_dir) + 1L), fixed = TRUE) &
                  endsWith(files, ".nc")]
  if (length(missing) || !length(ahoi)) {
    stop("The archive is missing required inputs:\n",
         paste(c(missing, if (!length(ahoi)) paste0(ahoi_dir, "*.nc")),
               collapse = "\n"))
  }

  # Only source files belong in the archive, not generated analysis products.
  rect_sidecar <- startsWith(files, paste0(rect_dir, rect_name, ".")) &
    !grepl("/", substring(files, nchar(rect_dir) + 1L), fixed = TRUE)
  unexpected <- files[!files %in% c(required, ahoi) & !rect_sidecar]
  if (length(unexpected)) {
    stop("Unexpected files in input archive:\n", paste(unexpected, collapse = "\n"))
  }

  # data_00.a takes the first run of digits in each filename as YYYYMM.
  date_text <- regmatches(basename(ahoi), regexpr("[0-9]+", basename(ahoi)))
  if (length(date_text) != length(ahoi) ||
      any(!grepl("^[0-9]{6}$", date_text))) {
    stop("Each AHOI filename must have YYYYMM as its first run of digits.")
  }
  dates <- as.Date(paste0(date_text, "01"), format = "%Y%m%d")
  if (anyNA(dates) || anyDuplicated(dates)) {
    stop("AHOI filenames must identify distinct, valid YYYYMM months.")
  }
  expected_months <- seq(min(dates), max(dates), by = "month")
  if (length(dates) != length(expected_months)) {
    warning("The AHOI filenames have gaps in their monthly coverage; ",
            "review the input before running the analysis.")
  }
  if (min(dates) != as.Date("1945-01-01") ||
      max(dates) != as.Date("2024-12-01")) {
    warning("AHOI filename coverage is ", format(min(dates), "%Y-%m"), " to ",
            format(max(dates), "%Y-%m"),
            "; downstream scripts expect a 1945-2024 raster. ",
            "Review coverage and filenames as described in README.md.")
  }

  message("Extracting and validating inputs in a temporary directory...")
  if (unzip_method == "internal") {
    message("Using R's internal unzip; individual archive members must be smaller than 4 GiB.")
    withCallingHandlers(
      utils::unzip(archive, files = files, exdir = staging),
      warning = function(w) stop("Could not extract ZIP archive: ", conditionMessage(w))
    )
  } else {
    # Batches avoid exceeding command-line length limits with monthly files.
    batches <- split(files, ceiling(seq_along(files) / 50))
    for (batch in batches) {
      output <- suppressWarnings(system2(
        unzip_program,
        c("-oq", shQuote(archive), shQuote(batch), "-d", shQuote(staging)),
        stdout = TRUE, stderr = TRUE
      ))
      status <- attr(output, "status")
      if (!is.null(status) && status != 0L) {
        stop("Could not extract ZIP archive:\n", paste(output, collapse = "\n"))
      }
    }
  }
  staged <- file.path(staging, files)
  staged_info <- file.info(staged)
  if (anyNA(staged_info$size) || any(staged_info$isdir) ||
      any(staged_info$size != file_sizes) ||
      any(nzchar(Sys.readlink(staged)))) {
    stop("Extracted inputs are missing, incomplete, or not regular files.")
  }

  # Refuse symlinks in existing destination paths and never overwrite data.
  for (path in files) {
    parts <- strsplit(path, "/", fixed = TRUE)[[1]]
    for (i in seq_along(parts)) {
      component <- do.call(file.path, as.list(parts[seq_len(i)]))
      link <- Sys.readlink(component)
      if (!is.na(link) && nzchar(link)) {
        stop("Destination contains a symbolic link: ", component)
      }
    }
  }
  existing <- file.exists(files)
  if (any(existing)) {
    same_size <- file.info(files[existing])$size == file_sizes[existing]
    same <- rep(FALSE, sum(existing))
    if (any(same_size)) {
      same[same_size] <- unname(tools::md5sum(files[existing][same_size])) ==
        unname(tools::md5sum(staged[existing][same_size]))
    }
    if (anyNA(same) || !all(same)) {
      stop("Existing input files differ from the archive; no files were overwritten:\n",
           paste(files[existing][is.na(same) | !same], collapse = "\n"))
    }
  }

  to_copy <- which(!existing)
  installed <- character()
  success <- FALSE
  on.exit(if (!success) unlink(installed), add = TRUE)
  for (i in to_copy) {
    parent <- dirname(files[i])
    if (!dir.exists(parent) &&
        !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create directory: ", parent)
    }
    if (!file.copy(staged[i], files[i], overwrite = FALSE)) {
      stop("Could not install input file: ", files[i])
    }
    installed <- c(installed, files[i])
  }
  success <- TRUE
  message("Inputs ready: ", length(to_copy), " installed, ",
          sum(existing), " already present and identical.")
}

download_inputs()
rm(download_inputs)
