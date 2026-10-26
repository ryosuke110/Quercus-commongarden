#!/usr/bin/env Rscript

# Build plot-level Qmon/Qser presence--absence data from Japan NFI files.
# Raw NFI files are supplied by the user and are never redistributed.

packages <- c("readr", "dplyr", "stringr", "tibble", "tidyr")
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall with: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript 00_build_occurrence_input.R ",
    "06_style3_large.csv 01_plot.csv [output_directory]",
    call. = FALSE
  )
}
tree_file <- args[[1]]
plot_file <- args[[2]]
output_dir <- if (length(args) >= 3L) args[[3]] else "data_work"

if (!file.exists(tree_file)) stop("Tree CSV not found: ", tree_file, call. = FALSE)
if (!file.exists(plot_file)) stop("Plot CSV not found: ", plot_file, call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

detect_encoding <- function(path) {
  guess <- suppressWarnings(readr::guess_encoding(path, n_max = 10000))
  if (!nrow(guess) || is.na(guess$encoding[[1]])) "UTF-8" else guess$encoding[[1]]
}

read_header <- function(path, encoding) {
  names(readr::read_csv(
    path,
    n_max = 0,
    locale = readr::locale(encoding = encoding),
    show_col_types = FALSE,
    name_repair = "minimal"
  ))
}

find_column <- function(header, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% header]
  if (length(hit)) return(hit[[1]])
  if (!required) return(NA_character_)
  stop(
    label, " column not found. Candidates: ", paste(candidates, collapse = ", "),
    "\nColumns present: ", paste(header, collapse = " | "),
    call. = FALSE
  )
}

clean_text <- function(x) {
  x <- stringr::str_replace_all(as.character(x), "\\u3000", " ")
  stringr::str_squish(x)
}

clean_id <- function(x) {
  x <- clean_text(x)
  x[x == ""] <- NA_character_
  x
}

tree_encoding <- detect_encoding(tree_file)
plot_encoding <- detect_encoding(plot_file)
tree_header <- read_header(tree_file, tree_encoding)
plot_header <- read_header(plot_file, plot_encoding)

id_candidates <- c("格子点ＩＤ", "格子点ID", "格子点ＩD", "格子点IＤ")
tree_id_col <- find_column(tree_header, id_candidates, "Tree plot ID")
species_col <- find_column(tree_header, c("樹種", "樹種名"), "Species")
plot_id_col <- find_column(plot_header, id_candidates, "Plot ID")
latitude_col <- find_column(
  plot_header, c("世界北緯１", "世界北緯1", "緯度", "latitude"), "Latitude"
)
longitude_col <- find_column(
  plot_header, c("世界東経１", "世界東経1", "経度", "longitude"), "Longitude"
)
land_use_col <- find_column(plot_header, c("土地利用", "土地利用区分"), "Land use")

trees <- readr::read_csv(
  tree_file,
  locale = readr::locale(encoding = tree_encoding),
  col_select = dplyr::all_of(c(tree_id_col, species_col)),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  progress = interactive(),
  name_repair = "minimal"
)
names(trees) <- c("plot_id", "species_jp")
trees <- trees |>
  dplyr::mutate(
    plot_id = clean_id(.data$plot_id),
    species_jp = clean_text(.data$species_jp)
  )

plots <- readr::read_csv(
  plot_file,
  locale = readr::locale(encoding = plot_encoding),
  col_select = dplyr::all_of(c(
    plot_id_col, latitude_col, longitude_col, land_use_col
  )),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  name_repair = "minimal"
)
names(plots) <- c("plot_id", "latitude", "longitude", "land_use")
plots <- plots |>
  dplyr::mutate(
    plot_id = clean_id(.data$plot_id),
    land_use = clean_text(.data$land_use),
    latitude = suppressWarnings(as.numeric(.data$latitude)),
    longitude = suppressWarnings(as.numeric(.data$longitude))
  ) |>
  dplyr::filter(.data$land_use == "森林", !is.na(.data$plot_id))

duplicate_ids <- plots |>
  dplyr::count(.data$plot_id, name = "n") |>
  dplyr::filter(.data$n > 1L)
if (nrow(duplicate_ids)) {
  stop(
    "Duplicated forest plot IDs in plot CSV. Examples: ",
    paste(utils::head(duplicate_ids$plot_id, 10), collapse = ", "),
    call. = FALSE
  )
}

focal_records <- trees |>
  dplyr::filter(.data$species_jp %in% c("ミズナラ", "コナラ")) |>
  dplyr::mutate(
    species_code = dplyr::if_else(.data$species_jp == "ミズナラ", "Qmon", "Qser")
  ) |>
  dplyr::select("plot_id", "species_jp", "species_code")

presence <- focal_records |>
  dplyr::filter(!is.na(.data$plot_id)) |>
  dplyr::distinct(.data$plot_id, .data$species_code) |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(
    names_from = "species_code",
    values_from = "present",
    values_fill = 0L
  )
if (!"Qmon" %in% names(presence)) presence$Qmon <- integer(nrow(presence))
if (!"Qser" %in% names(presence)) presence$Qser <- integer(nrow(presence))

# The forest-plot table is the master, so plots where neither focal species was
# recorded are retained as valid double absences.
occurrence <- plots |>
  dplyr::select("plot_id", "longitude", "latitude") |>
  dplyr::left_join(presence, by = "plot_id") |>
  dplyr::mutate(
    Qmon = dplyr::coalesce(as.integer(.data$Qmon), 0L),
    Qser = dplyr::coalesce(as.integer(.data$Qser), 0L)
  ) |>
  dplyr::arrange(.data$plot_id) |>
  dplyr::select("plot_id", "Qmon", "Qser", "longitude", "latitude")

if (anyDuplicated(occurrence$plot_id)) stop("Internal error: plot_id is not unique.")
if (any(!occurrence$Qmon %in% 0:1) || any(!occurrence$Qser %in% 0:1)) {
  stop("Internal error: Qmon/Qser contains values other than 0 and 1.")
}

focal_ids <- unique(focal_records$plot_id[!is.na(focal_records$plot_id)])
unmatched_ids <- setdiff(focal_ids, occurrence$plot_id)
bad_coordinates <- is.na(occurrence$longitude) | is.na(occurrence$latitude) |
  occurrence$longitude < 122 | occurrence$longitude > 154 |
  occurrence$latitude < 20 | occurrence$latitude > 46.5

audit <- tibble::tribble(
  ~metric, ~value, ~detail,
  "tree_file", NA_real_, normalizePath(tree_file),
  "plot_file", NA_real_, normalizePath(plot_file),
  "tree_encoding", NA_real_, tree_encoding,
  "plot_encoding", NA_real_, plot_encoding,
  "forest_plots", nrow(occurrence), "master analysis population",
  "Qmon_tree_records", sum(focal_records$species_code == "Qmon"), "ミズナラ records",
  "Qser_tree_records", sum(focal_records$species_code == "Qser"), "コナラ records",
  "Qmon_presence_plots", sum(occurrence$Qmon), "plot-level presence",
  "Qser_presence_plots", sum(occurrence$Qser), "plot-level presence",
  "cooccurrence_plots", sum(occurrence$Qmon == 1L & occurrence$Qser == 1L), "same 0.1-ha plot",
  "neither_plots", sum(occurrence$Qmon == 0L & occurrence$Qser == 0L), "retained absences",
  "focal_ids_not_in_forest_master", length(unmatched_ids),
  paste(utils::head(unmatched_ids, 20), collapse = ";"),
  "invalid_or_missing_coordinates", sum(bad_coordinates), "outside Japan bounds or NA"
)

readr::write_csv(occurrence, file.path(output_dir, "Qmon_Qser_occurrence.csv"), na = "")
readr::write_csv(audit, file.path(output_dir, "occurrence_audit.csv"), na = "NA")

message("Occurrence input created: ",
        normalizePath(file.path(output_dir, "Qmon_Qser_occurrence.csv")))
message("Forest plots: ", nrow(occurrence))
message("Qmon/Qser/co-occurrence: ", sum(occurrence$Qmon), "/",
        sum(occurrence$Qser), "/",
        sum(occurrence$Qmon == 1L & occurrence$Qser == 1L))

