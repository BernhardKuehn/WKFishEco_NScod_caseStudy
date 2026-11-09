# Download the external inputs. Run from the repository root.

download_inputs <- function() {
  if (!file.exists("data_00.a_Handle_AHOI_data.R") ||
      !dir.exists("functions")) {
    stop("Run this script from the repository root.")
  }
  
  selection <- Sys.getenv("DOWNLOAD_INPUTS", unset = "all")
  if (!selection %in% c("all", "figshare", "ices")) {
    stop("DOWNLOAD_INPUTS must be 'all', 'figshare', or 'ices'.")
  }
  
  figshare_files <- c(
    FIGSHARE_RCP45_URL = "data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp45_BottomTemp_raster_2006_2099.RData",
    FIGSHARE_RCP85_URL = "data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp85_BottomTemp_raster_2006_2099.RData",
    FIGSHARE_AHOI_RASTER_URL = "data/Env/processed/AHOI_BottomTemp_raster_1945_2024.tif",
    FIGSHARE_STOCK_AREAS_URL = "data/stock/Area_NSCOD3_table_6a_inshore_offshore.xlsx",
    FIGSHARE_STOCK_OBJECTS_URL = "data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata"
  )
  
  figshare_urls <- Sys.getenv(names(figshare_files), unset = "")
  if (selection != "ices" && any(!nzchar(figshare_urls))) {
    stop(
      "Set the missing Figshare file download URLs: ",
      paste(names(figshare_files)[!nzchar(figshare_urls)], collapse = ", "),
      ". For ICES only, set DOWNLOAD_INPUTS=ices."
    )
  }
  
  ices_url <- Sys.getenv(
    "ICES_RECTANGLES_URL",
    unset = "https://gis.ices.dk/shapefiles/ICES_rectangles.zip"
  )
  
  urls <- c(
    if (selection != "ices") figshare_urls,
    if (selection != "figshare") ices_url
  )
  if (any(!grepl("^(https|file)://", urls, ignore.case = TRUE))) {
    stop("Download URLs must use HTTPS or local file:// URLs.")
  }
  
  staging <- tempfile(pattern = "wkfisheco-inputs-")
  if (!dir.create(staging)) stop("Could not create a staging directory.")
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  
  old_timeout <- getOption("timeout")
  options(timeout = max(1800, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  
  download <- function(url, path, label) {
    message("Downloading ", label, "...")
    status <- utils::download.file(url, path, mode = "wb", quiet = TRUE)
    if (status != 0L || !file.exists(path) ||
        is.na(file.info(path)$size) || file.info(path)$size <= 0) {
      stop("Download failed or produced an empty file: ", label)
    }
  }
  
  staged <- character()
  destinations <- character()
  
  if (selection != "ices") {
    for (i in seq_along(figshare_files)) {
      path <- file.path(staging, paste0("figshare-", i))
      download(figshare_urls[i], path, basename(figshare_files[i]))
      staged <- c(staged, path)
      destinations <- c(destinations, unname(figshare_files[i]))
    }
  }
  
  if (selection != "figshare") {
    archive <- file.path(staging, "ices.zip")
    download(ices_url, archive, "ICES statistical rectangles")
    
    unzip_program <- Sys.which("unzip")
    unzip_method <- if (nzchar(unzip_program)) unzip_program else "internal"
    contents <- withCallingHandlers(
      utils::unzip(archive, list = TRUE, unzip = unzip_method),
      warning = function(w) stop("Cannot read ICES ZIP: ", conditionMessage(w))
    )
    
    paths <- contents$Name
    if (!length(paths) || anyNA(paths) || anyDuplicated(paths) ||
        any(grepl("\\\\|(^/)|(^|/)\\.\\.?(/|$)|^[A-Za-z]:|[?*\\[\"\r\n]", paths))) {
      stop("The ICES ZIP contains empty, duplicate, or unsafe paths.")
    }
    
    files <- paths[!endsWith(paths, "/")]
    for (path in files) {
      if (any(startsWith(paths, paste0(path, "/")))) {
        stop("The ICES ZIP contains a file used as a directory: ", path)
      }
    }
    
    dbfs <- files[grepl("\\.dbf$", files, ignore.case = TRUE)]
    if (!length(dbfs)) {
      stop("The ICES ZIP contains no shapefile attribute table.")
    }
    
    extract <- function(members, directory) {
      if (unzip_method == "internal") {
        withCallingHandlers(
          utils::unzip(archive, files = members, exdir = directory),
          warning = function(w) {
            stop("Could not extract ICES ZIP: ", conditionMessage(w))
          }
        )
      } else {
        output <- suppressWarnings(system2(
          unzip_program,
          c("-oq", shQuote(archive), shQuote(members), "-d", shQuote(directory)),
          stdout = TRUE,
          stderr = TRUE
        ))
        status <- attr(output, "status")
        if (!is.null(status) && status != 0L) {
          stop("Could not extract ICES ZIP:\n", paste(output, collapse = "\n"))
        }
      }
      
      extracted <- file.path(directory, members)
      expected <- contents$Length[match(members, paths)]
      info <- file.info(extracted)
      links <- Sys.readlink(extracted)
      if (anyNA(info$size) || any(info$isdir) || anyNA(expected) ||
          any(expected <= 0) || any(info$size != expected) ||
          any(!is.na(links) & nzchar(links))) {
        stop("ICES ZIP contains missing, incomplete, or non-regular files.")
      }
      extracted
    }
    
    # Read DBF field names without requiring a spatial R package.
    dbf_fields <- function(path) {
      con <- file(path, "rb")
      on.exit(close(con))
      
      header <- readBin(con, "raw", n = 32L)
      if (length(header) != 32L) return(character())
      
      header_length <- as.integer(header[9]) + 256L * as.integer(header[10])
      if (header_length < 33L || header_length > file.info(path)$size ||
          (header_length - 33L) %% 32L != 0L) {
        return(character())
      }
      
      count <- (header_length - 33L) / 32L
      fields <- character(count)
      for (i in seq_len(count)) {
        descriptor <- readBin(con, "raw", n = 32L)
        if (length(descriptor) != 32L) return(character())
        name <- descriptor[1:11]
        zero <- match(as.raw(0), name, nomatch = 12L)
        fields[i] <- rawToChar(name[seq_len(zero - 1L)])
      }
      fields
    }
    
    probe_dir <- file.path(staging, "dbf")
    dir.create(probe_dir)
    dbf_paths <- extract(dbfs, probe_dir)
    matches <- vapply(dbf_paths, function(path) {
      all(c("ICESNAME", "Area_27") %in% dbf_fields(path))
    }, logical(1))
    
    if (sum(matches) != 1L) {
      stop(
        "Expected exactly one ICES shapefile with ICESNAME and Area_27 ",
        "fields; found ", sum(matches), ". The ICES layer may not be ",
        "compatible with data_01_prepare_bottomT_by_stubstock.R."
      )
    }
    
    selected_dbf <- dbfs[matches]
    stem <- substr(selected_dbf, 1L, nchar(selected_dbf) - 4L)
    same_stem <- startsWith(files, paste0(stem, "."))
    extensions <- tolower(substring(files, nchar(stem) + 2L))
    allowed <- c("shp", "shx", "dbf", "prj", "cpg", "sbn", "sbx", "qix")
    members <- files[same_stem & extensions %in% allowed]
    member_ext <- extensions[same_stem & extensions %in% allowed]
    
    if (anyDuplicated(member_ext) ||
        !all(c("shp", "shx", "dbf", "prj") %in% member_ext)) {
      stop(
        "The compatible ICES shapefile is missing required sidecars or ",
        "contains duplicate extensions."
      )
    }
    
    rect_dir <- file.path(staging, "rectangles")
    dir.create(rect_dir)
    staged <- c(staged, extract(members, rect_dir))
    destinations <- c(destinations, paste0(
      "data/ICES_stat_rect/StatRec_map_Areas_Full_20170124.", member_ext
    ))
  }
  
  for (path in destinations) {
    parts <- strsplit(path, "/", fixed = TRUE)[[1]]
    for (i in seq_along(parts)) {
      component <- do.call(file.path, as.list(parts[seq_len(i)]))
      link <- Sys.readlink(component)
      if (!is.na(link) && nzchar(link)) {
        stop("Destination contains a symbolic link: ", component)
      }
    }
  }
  
  existing <- file.exists(destinations)
  if (any(existing)) {
    same <- unname(tools::md5sum(destinations[existing])) ==
      unname(tools::md5sum(staged[existing]))
    if (anyNA(same) || !all(same)) {
      stop(
        "Existing input files differ; no files were overwritten:\n",
        paste(destinations[existing][is.na(same) | !same], collapse = "\n")
      )
    }
  }
  
  installed <- character()
  success <- FALSE
  on.exit(if (!success) unlink(installed), add = TRUE)
  
  for (i in which(!existing)) {
    parent <- dirname(destinations[i])
    if (!dir.exists(parent) &&
        !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
      stop("Could not create directory: ", parent)
    }
    if (!file.copy(staged[i], destinations[i], overwrite = FALSE)) {
      stop("Could not install input file: ", destinations[i])
    }
    installed <- c(installed, destinations[i])
  }
  
  success <- TRUE
  message(
    "Inputs ready: ", length(installed), " installed, ",
    sum(existing), " already present and identical."
  )
}

download_inputs()
rm(download_inputs)

# download data from ICES (statistical rectangles)
Sys.setenv(DOWNLOAD_INPUTS = "ices")
source("data_00_download_inputs.R")

# Download figshare data
Sys.unsetenv("DOWNLOAD_INPUTS")
Sys.setenv(
  FIGSHARE_RCP45_URL = "https://figshare.com/ndownloader/files/RCP45_FILE_ID",
  FIGSHARE_RCP85_URL = "https://figshare.com/ndownloader/files/RCP85_FILE_ID",
  FIGSHARE_AHOI_RASTER_URL = "https://figshare.com/ndownloader/files/AHOI_FILE_ID",
  FIGSHARE_STOCK_AREAS_URL = "https://figshare.com/ndownloader/files/AREAS_FILE_ID",
  FIGSHARE_STOCK_OBJECTS_URL = "https://figshare.com/ndownloader/files/STOCK_FILE_ID"
)
source("data_00_download_inputs.R")


