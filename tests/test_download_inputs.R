# Run from the repository root with Rscript tests/test_download_inputs.R.
# Uses base R and the system zip utility; no real input data are needed.
stopifnot(file.exists("data_00_download_inputs.R"), nzchar(Sys.which("zip")))
script <- normalizePath("data_00_download_inputs.R")
root <- tempfile("download-input-tests-")
dir.create(root)

make_workspace <- function(name) {
  path <- file.path(root, name)
  dir.create(path)
  dir.create(file.path(path, "functions"))
  file.create(file.path(path, "data_00.a_Handle_AHOI_data.R"))
  file.copy(script, file.path(path, basename(script)))
  path
}

source_dir <- file.path(root, "source")
dir.create(source_dir)
rect <- "data/ICES_stat_rect/StatRec_map_Areas_Full_20170124"
files <- c(
  "data/Env/AHOI/v23.05/194501.nc",
  "data/Env/AHOI/v23.05/194502.nc",
  paste0(rect, c(".shp", ".shx", ".dbf", ".prj", ".cpg")),
  "data/stock/Area_NSCOD3_table_6a_inshore_offshore.xlsx",
  "data/stock/cod46a7d20_FLstockObjects_assessment_2026.Rdata",
  "data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp45_BottomTemp_raster_2006_2099.RData",
  "data/Env/POLCOMS-ERSEM/POLCOM_ERSEM_GreaterNorthSea_rcp85_BottomTemp_raster_2006_2099.RData"
)
for (path in files) {
  full <- file.path(source_dir, path)
  dir.create(dirname(full), recursive = TRUE, showWarnings = FALSE)
  writeLines(paste("fixture", path), full)
}

make_zip <- function(name, contents) {
  archive <- file.path(root, paste0(name, ".zip"))
  old <- setwd(source_dir)
  on.exit(setwd(old))
  status <- utils::zip(archive, contents, flags = "-q")
  stopifnot(status == 0L)
  archive
}

run_downloader <- function(workspace, archive, without_unzip = FALSE) {
  old <- setwd(workspace)
  on.exit(setwd(old))
  env <- paste0("FIGSHARE_INPUT_URL=file://", archive)
  if (without_unzip) env <- c(env, paste0("PATH=", root))
  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    "data_00_download_inputs.R",
    stdout = TRUE, stderr = TRUE,
    env = env
  ))
  status <- attr(result, "status")
  list(status = if (is.null(status)) 0L else status,
       output = paste(result, collapse = "\n"))
}

test <- function(condition, description) {
  if (!isTRUE(condition)) stop("Failed: ", description)
  cat("OK:", description, "\n")
}

archive <- make_zip("complete", files)
workspace <- make_workspace("fresh")
result <- run_downloader(workspace, archive)
test(result$status == 0L && all(file.exists(file.path(workspace, files))),
     "fresh workspace receives every input")
result <- run_downloader(workspace, archive)
test(result$status == 0L && grepl("already present and identical", result$output),
     "repeat run preserves identical inputs")

workspace_internal <- make_workspace("internal-unzip")
result <- run_downloader(workspace_internal, archive, without_unzip = TRUE)
test(result$status == 0L && all(file.exists(file.path(workspace_internal, files))),
     "base R extraction works when unzip is unavailable")

zip64_marker <- file.path(root, "zip64-marker.zip")
stopifnot(file.copy(archive, zip64_marker))
con <- file(zip64_marker, "ab")
writeBin(as.raw(c(0x50, 0x4b, 0x06, 0x07)), con)
close(con)
workspace <- make_workspace("zip64-without-unzip")
result <- run_downloader(workspace, zip64_marker, without_unzip = TRUE)
test(result$status != 0L && grepl("ZIP64", result$output) &&
       !dir.exists(file.path(workspace, "data")),
     "ZIP64 is rejected when only base R extraction is available")

months <- format(seq(as.Date("1945-01-01"), by = "month", length.out = 60), "%Y%m")
many_ahoi <- paste0("data/Env/AHOI/v23.05/", months, ".nc")
for (path in setdiff(many_ahoi, files)) {
  writeLines("fixture", file.path(source_dir, path))
}
many_files <- c(many_ahoi, files[-c(1, 2)])
workspace_many <- make_workspace("many-months")
result <- run_downloader(workspace_many, make_zip("many-months", many_files))
test(result$status == 0L && all(file.exists(file.path(workspace_many, many_files))),
     "multiple extraction batches install all monthly files")

workspace <- file.path(root, "fresh")
changed <- file.path(workspace, files[1])
absent <- file.path(workspace, files[2])
writeLines("local change", changed)
unlink(absent)
result <- run_downloader(workspace, archive)
test(result$status != 0L && identical(readLines(changed), "local change") &&
       !file.exists(absent), "conflicts stop before installing or overwriting files")

workspace <- make_workspace("incomplete")
result <- run_downloader(workspace, make_zip("incomplete", files[-3]))
test(result$status != 0L && grepl("missing required inputs", result$output) &&
       !dir.exists(file.path(workspace, "data")),
     "missing required input fails before installation")

workspace <- make_workspace("unexpected")
extra <- "data/generated.rda"
writeLines("generated", file.path(source_dir, extra))
result <- run_downloader(workspace, make_zip("unexpected", c(files, extra)))
test(result$status != 0L && grepl("Unexpected files", result$output) &&
       !dir.exists(file.path(workspace, "data")),
     "generated or unrelated files are rejected")

workspace <- make_workspace("invalid-date")
bad_date <- "data/Env/AHOI/v23.05/version2_194501.nc"
writeLines("fixture", file.path(source_dir, bad_date))
result <- run_downloader(workspace, make_zip("invalid-date", c(files[-1], bad_date)))
test(result$status != 0L && grepl("first run of digits", result$output),
     "AHOI filename must match the processing script's date convention")

workspace <- make_workspace("duplicate-month")
duplicate <- "data/Env/AHOI/v23.05/194501_other.nc"
writeLines("fixture", file.path(source_dir, duplicate))
result <- run_downloader(workspace, make_zip("duplicate-month", c(files, duplicate)))
test(result$status != 0L && grepl("distinct, valid YYYYMM", result$output),
     "duplicate AHOI months are rejected")

workspace <- make_workspace("bad-zip")
bad_zip <- file.path(root, "bad.zip")
writeLines("not a zip", bad_zip)
result <- run_downloader(workspace, bad_zip)
test(result$status != 0L && !dir.exists(file.path(workspace, "data")),
     "invalid ZIP does not install files")

workspace <- make_workspace("unsafe")
writeLines("escape", file.path(root, "escape"))
result <- run_downloader(workspace, make_zip("unsafe", c(files, "../escape")))
test(result$status != 0L && grepl("unsafe paths", result$output) &&
       !file.exists(file.path(workspace, "escape")),
     "archive path traversal is rejected")

workspace <- make_workspace("symlink")
outside <- file.path(root, "outside")
dir.create(outside)
if (file.symlink(outside, file.path(workspace, "data"))) {
  result <- run_downloader(workspace, archive)
  test(result$status != 0L && grepl("symbolic link", result$output) &&
         !length(list.files(outside, all.files = FALSE)),
       "destination symlinks are rejected")
}

unlink(root, recursive = TRUE)
cat("All downloader tests passed.\n")
