# Download the external inputs. Run from the repository root.

download_inputs <- function(
    figshare_doi = "https://doi.org/10.6084/m9.figshare.34961302",
    ices_url = "https://gis.ices.dk/shapefiles/ICES_StatRec_mapto_ICES_Areas.zip",
    figshare = TRUE, ices = TRUE, figshare_api_url = NULL) {
  
  library(jsonlite)
  
  if (!file.exists("data_00.a_Handle_AHOI_data.R") ||
      !dir.exists("functions")) {
    stop("Run this script from the repository root.")
  }
  
  if (!is.logical(figshare) || length(figshare) != 1L || is.na(figshare) ||
      !is.logical(ices) || length(ices) != 1L || is.na(ices) ||
      (!figshare && !ices)) {
    stop("Select at least one source using figshare = TRUE or ices = TRUE.")
  }
  
  if (figshare && !requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "The jsonlite package is required to read Figshare metadata. ",
      "Install it with install.packages('jsonlite')."
    )
  }
  
  if (figshare && is.null(figshare_api_url)) {
    if (length(figshare_doi) != 1L || is.na(figshare_doi) ||
        !grepl("^https://doi\\.org/10\\.6084/m9\\.figshare\\.[0-9]+$",
               figshare_doi)) {
      stop("Expected an unversioned Figshare DOI URL.")
    }
    figshare_api_url <- paste0(
      "https://api.figshare.com/v2/articles/",
      sub(".*\\.", "", figshare_doi)
    )
  }
  
  urls <- c(if (figshare) figshare_api_url, if (ices) ices_url)
  if (anyNA(urls) || any(!grepl("^(https|file)://", urls, ignore.case = TRUE))) {
    stop("Source URLs must use HTTPS or local file:// URLs.")
  }
  
  staging <- tempfile(pattern = "wkfisheco-inputs-")
  if (!dir.create(staging)) stop("Could not create a staging directory.")
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  
  old_timeout <- getOption("timeout")
  options(timeout = max(1800, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  
  download <- function(url, path, label, method = "auto") {
    message("Downloading ", label, "...")
    status <- utils::download.file(
      url, path, mode = "wb", quiet = TRUE, method = method
    )
    if (status != 0L || !file.exists(path) ||
        is.na(file.info(path)$size) || file.info(path)$size <= 0) {
      stop("Download failed or produced an empty file: ", label)
    }
  }
  
  staged <- character()
  destinations <- character()
  
  if (figshare) {
    figshare_files <- c(
      "data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp45_BottomTemp_raster_2006_2099.RData",
      "data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp85_BottomTemp_raster_2006_2099.RData",
      "data/Env/processed/AHOI_BottomTemp_raster_1945_2024.tif",
      "data/stock/Area_NSCOD3_table_6a_inshore_offshore.xlsx",
      "data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata"
    )
    
    aliases <- lapply(basename(figshare_files), identity)
    aliases[[3]] <- c(
      aliases[[3]],
      "AHOI_BottomTempRaster_1945_2024.tif"
    )
    
    metadata <- file.path(staging, "figshare.json")
    download(figshare_api_url, metadata, "Figshare article metadata")
    
    article <- tryCatch(
      jsonlite::fromJSON(metadata, simplifyVector = FALSE),
      error = function(e) {
        stop("Could not read Figshare metadata: ", conditionMessage(e))
      }
    )
    
    if (!is.list(article) || !is.list(article$files) || !length(article$files)) {
      stop("Figshare metadata contains no file list.")
    }
    
    file_names <- vapply(article$files, function(file) {
      if (is.list(file) && is.character(file$name) &&
          length(file$name) == 1L && !is.na(file$name)) {
        file$name
      } else {
        NA_character_
      }
    }, character(1))
    
    if (anyNA(file_names) || any(!nzchar(file_names))) {
      stop("Figshare metadata contains an invalid filename.")
    }
    
    selected <- lapply(seq_along(figshare_files), function(i) {
      matches <- which(tolower(file_names) %in% tolower(aliases[[i]]))
      
      if (length(matches) != 1L) {
        stop(
          "Expected exactly one Figshare file matching ",
          paste(aliases[[i]], collapse = " or "),
          "; found ", length(matches), "."
        )
      }
      
      file <- article$files[[matches]]
      if (!is.character(file$download_url) ||
          length(file$download_url) != 1L ||
          is.na(file$download_url) ||
          !grepl("^(https|file)://", file$download_url, ignore.case = TRUE) ||
          !is.numeric(file$size) || length(file$size) != 1L ||
          is.na(file$size) || !is.finite(file$size) || file$size <= 0 ||
          !is.character(file$computed_md5) ||
          length(file$computed_md5) != 1L ||
          is.na(file$computed_md5) ||
          !grepl("^[0-9a-fA-F]{32}$", file$computed_md5)) {
        stop(
          "Figshare metadata has an invalid URL, size or checksum for ",
          file_names[matches], "."
        )
      }
      file
    })
    
    for (i in seq_along(figshare_files)) {
      path <- file.path(staging, paste0("figshare-", i))
      download(selected[[i]]$download_url, path, basename(figshare_files[i]))
      
      checksum <- unname(tools::md5sum(path))
      if (is.na(checksum) || file.info(path)$size != selected[[i]]$size ||
          tolower(checksum) != tolower(selected[[i]]$computed_md5)) {
        stop(
          "Downloaded Figshare file failed size or checksum validation: ",
          basename(figshare_files[i])
        )
      }
      
      staged <- c(staged, path)
      destinations <- c(destinations, figshare_files[i])
    }
  }
  
  if (ices) {
    archive <- file.path(staging, "ices.zip")
    download(ices_url, archive, "ICES statistical rectangles",
             method = "libcurl")
    
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