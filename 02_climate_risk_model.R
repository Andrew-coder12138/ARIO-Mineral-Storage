# Standalone climate risk ARIO simulation.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
project_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
Sys.setenv(ARIO_PROJECT_ROOT = project_dir)
# Input processing, material constraints and simulation methods.

# --- Initialization and shared utilities ---
options(scipen = 999, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(matrixStats)
  library(readxl)
  library(stringr)
})

project_root <- normalizePath(Sys.getenv("ARIO_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
base_data_path <- file.path(project_root, "data")
z_path <- file.path(base_data_path, "io", "transactions.csv")
y_path <- file.path(base_data_path, "io", "final_demand.csv")
v_path <- file.path(base_data_path, "io", "value_added.csv")
tq_path <- file.path(base_data_path, "io", "emissions.csv")
battery_forecast_path <- file.path(base_data_path, "demand", "storage_demand.csv")
main_scenario_path <- file.path(base_data_path, "production")
out_root <- file.path(project_root, "outputs")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
task_log_file <- file.path(out_root, "task_status.csv")
error_log_file <- file.path(out_root, "run_log.txt")
# Each model stores results in its own output directory.

model_cfg <- fread(file.path(project_root, "config", "model_parameters.csv"))
stopifnot(all(c("Parameter", "Value") %in% names(model_cfg)))
get_param <- function(nm) {
  x <- model_cfg[Parameter == nm, Value]
  if (length(x) != 1L) stop("Missing or duplicate model parameter: ", nm)
  v <- suppressWarnings(as.numeric(x))
  if (!is.finite(v)) stop("Non-numeric model parameter: ", nm)
  v
}
START_YEAR <- as.integer(get_param("start_year"))
END_YEAR <- as.integer(get_param("end_year"))
years <- seq.int(START_YEAR, END_YEAR)
years_cols <- as.character(years)
TT <- length(years)
CAPACITY_FACTOR <- get_param("capacity_factor")
STOCK_MONTHS <- get_param("stock_months")
STOCK_PARA_VALUE <- STOCK_MONTHS / 12
ENABLE_MINERAL_SUPPLY_EXPANSION <- TRUE
STOP_IF_BATTERY_FORECAST_UNMATCHED <- TRUE
USE_GWH_DIRECTLY_AS_DEMAND <- FALSE
stopifnot(TT > 1L, CAPACITY_FACTOR > 0, STOCK_MONTHS > 0)

log_line <- function(text) {
  line <- sprintf("[%s] %s", Sys.time(), text)
  cat(line, "\n")
  cat(line, "\n", file = error_log_file, append = TRUE)
}

append_task_status <- function(macro, tech, status, message = "", mode = "baseline") {
  dt <- data.table(Time = as.character(Sys.time()), Model = mode, Macro = macro,
                   Tech = tech, Status = status, Message = message)
  fwrite(dt, task_log_file, append = file.exists(task_log_file))
}

check_file_exists <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing input file: %s", path))
}

safe_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) stop("Matrix contains non-finite values")
  x
}

# --- Input data and validation ---
for (p in c(z_path, y_path, v_path, tq_path, battery_forecast_path)) check_file_exists(p)
if (!dir.exists(main_scenario_path)) stop(sprintf("Missing scenario directory: %s", main_scenario_path))


log_line("Loading input-output tables")

data_Z <- fread(z_path, header = TRUE)
labels_raw <- trimws(as.character(data_Z[[1]]))
input_cols <- trimws(names(data_Z)[-1])
if (!identical(labels_raw, input_cols)) stop("Input-output matrix row and column labels do not match")
Z_raw <- safe_numeric_matrix(data_Z[, -1, with = FALSE])
rm(data_Z); gc()

Y_tab <- fread(y_path)
if (!identical(labels_raw, as.character(Y_tab[[1]]))) stop("Final demand sector order mismatch")
Y_raw <- safe_numeric_matrix(Y_tab[, -1, with = FALSE])

V_tab <- fread(v_path)
if (!identical(labels_raw, names(V_tab)[-1])) stop("Value added sector order mismatch")
V_tmp <- safe_numeric_matrix(V_tab[, -1, with = FALSE])
V_raw <- matrix(colSums(V_tmp, na.rm = TRUE), nrow = 1)
rm(V_tmp); gc()

TQ_tab <- fread(tq_path)
if (!identical(labels_raw, names(TQ_tab)[-1])) stop("Emissions sector order mismatch")
TQ_tmp <- safe_numeric_matrix(TQ_tab[, -1, with = FALSE])
TQ_raw <- matrix(colSums(TQ_tmp, na.rm = TRUE), nrow = 1)
rm(TQ_tmp); gc()

# --- Input consistency checks ---
if (length(labels_raw) != length(unique(labels_raw))) stop("Duplicate sector labels")
if (nrow(Z_raw) != ncol(Z_raw)) stop("Transaction matrix must be square")
if (nrow(Y_raw) != length(labels_raw)) stop("Final demand dimension mismatch")
if (ncol(V_raw) != length(labels_raw) || ncol(TQ_raw) != length(labels_raw)) stop("Satellite account dimension mismatch")
if (any(Z_raw < 0) || any(Y_raw < 0) || any(V_raw < 0) || any(TQ_raw < 0)) stop("Input matrices contain negative values")
raw_parts <- strsplit(labels_raw, ".", fixed = TRUE)
if (any(lengths(raw_parts) < 2)) stop("Sector labels must include a region prefix")
if (length(unique(vapply(raw_parts, `[`, character(1), 1))) < 2L) stop("Insufficient number of regions")

sector_meta <- fread(file.path(project_root, "config", "sector_dictionary.csv"))
model_sector_raw <- sub("^[^.]+\\.", "", labels_raw)
if (!all(unique(model_sector_raw) %in% sector_meta$Sector_Code)) stop("Sector dictionary is incomplete")

# --- Sector classification and mapping ---
battery_sectors <- c("BAT_LIB", "BAT_RFB", "BAT_LAB", "BAT_MAB", "BAT_TSB")

mineral_sectors <- c(
  "S36.LI", "S36.CO", "S34.NI", "S77.MN", "S77.V", "S77.CR", "S40.GR",
  "S31.CU", "S34.AL", "S77.FE", "S34.PB", "S36.NA", "S36.REE",
  "S34.SIL", "S34.PGM", "S34.NB", "S34.W", "S36.BE", "S34.GA",
  "S34.GE", "S34.HF"
)

electrified_sectors <- c("ELEC_ROAD", "ELEC_RAIL", "ELEC_POWER", "ELEC_BUILDING", "ELEC_INDUSTRY")
traditional_sectors <- c("CONV_ROAD", "CONV_RAIL", "CONV_POWER", "CONV_BUILDING", "CONV_INDUSTRY")
sector_export <- c(electrified_sectors, traditional_sectors)

keep_sectors <- c(mineral_sectors, battery_sectors, electrified_sectors, traditional_sectors)

critical_mineral_to_battery <- list(
  BAT_LIB = c("S36.LI", "S36.CO", "S34.NI", "S77.MN", "S40.GR", "S31.CU", "S34.AL"),
  BAT_RFB = c("S77.V", "S36.LI", "S40.GR", "S31.CU"),
  BAT_LAB = c("S34.PB", "S31.CU", "S34.AL"),
  BAT_MAB = c("S36.NA", "S77.MN", "S40.GR", "S31.CU", "S34.AL"),
  BAT_TSB = c("S36.LI", "S36.NA", "S36.REE", "S34.NI", "S31.CU", "S34.AL")
)

safe_num <- function(x) {
  x <- as.numeric(x)
  x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
  x
}

parse_region_dot <- function(labels) {
  sapply(strsplit(labels, "\\."), `[`, 1)
}

parse_sector_after_first_dot <- function(labels) {
  sapply(strsplit(labels, "\\."), function(x) {
    if (length(x) <= 1) return(NA_character_)
    paste(x[-1], collapse = ".")
  })
}

parse_last_token <- function(labels) {
  sapply(strsplit(labels, "[\\._]"), function(x) tail(x, 1))
}

extract_sector <- function(labels) {
  sapply(strsplit(labels, "\\."), function(x) paste(x[-1], collapse = "."))
}

normalize_macro_scenario <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- gsub("/", "_", x, fixed = TRUE)
  x <- gsub("\\\\", "_", x)
  x <- gsub("_results_output$", "", x, ignore.case = TRUE)
  x <- gsub("_result_output$", "", x, ignore.case = TRUE)
  x <- gsub("_results$", "", x, ignore.case = TRUE)
  x <- gsub("_result$", "", x, ignore.case = TRUE)
  trimws(x)
}

normalize_tech_scenario <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- gsub("\\s*\\(", " (", x)
  x <- gsub("\\)\\s*", ")", x)
  trimws(x)
}

normalize_to_model_sector <- function(full_sector, last_token) {
  out <- full_sector

  idx_bat <- last_token %in% battery_sectors
  out[idx_bat] <- last_token[idx_bat]

  idx_elec <- last_token %in% electrified_sectors
  out[idx_elec] <- last_token[idx_elec]

  idx_trad <- last_token %in% traditional_sectors
  out[idx_trad] <- last_token[idx_trad]

  for (ms in mineral_sectors) {
    tk <- tail(strsplit(ms, "\\.")[[1]], 1)
    idx <- !(full_sector %in% mineral_sectors) & last_token == tk
    out[idx] <- ms
  }

  out
}

safe_clean_filename <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[^[:alnum:]_. -]+", "_", x, perl = TRUE)
  x <- gsub("[[:space:]]+", "_", x, perl = TRUE)
  x <- gsub("_+", "_", x, perl = TRUE)
  x <- gsub("^_+|_+$", "", x, perl = TRUE)
  ifelse(nchar(x) == 0, "unnamed", substr(x, 1, 160))
}

get_tech_code <- function(x) {
  code <- stringr::str_extract(as.character(x), "S[0-9]+")
  ifelse(is.na(code), as.character(x), code)
}

mineral_code_alias <- c(
  "LI" = "S36.LI", "LITHIUM" = "S36.LI",
  "CO" = "S36.CO", "COBALT" = "S36.CO",
  "NI" = "S34.NI", "NICKEL" = "S34.NI",
  "MN" = "S77.MN", "MANGANESE" = "S77.MN",
  "V" = "S77.V", "VANADIUM" = "S77.V",
  "CR" = "S77.CR", "CHROMIUM" = "S77.CR",
  "GR" = "S40.GR", "GRAPHITE" = "S40.GR",
  "CU" = "S31.CU", "COPPER" = "S31.CU",
  "AL" = "S34.AL", "ALUMINUM" = "S34.AL", "ALUMINIUM" = "S34.AL",
  "FE" = "S77.FE", "IRON" = "S77.FE",
  "PB" = "S34.PB", "LEAD" = "S34.PB",
  "NA" = "S36.NA", "SODIUM" = "S36.NA",
  "REE" = "S36.REE", "RAREEARTH" = "S36.REE",
  "SIL" = "S34.SIL", "SILICON" = "S34.SIL",
  "PGM" = "S34.PGM",
  "NB" = "S34.NB", "NIOBIUM" = "S34.NB",
  "W" = "S34.W", "TUNGSTEN" = "S34.W",
  "BE" = "S36.BE", "BERYLLIUM" = "S36.BE",
  "GA" = "S34.GA", "GALLIUM" = "S34.GA",
  "GE" = "S34.GE", "GERMANIUM" = "S34.GE",
  "HF" = "S34.HF", "HAFNIUM" = "S34.HF"
)

standardize_mapped_code <- function(x) {
  x0 <- as.character(x)
  x0 <- trimws(x0)
  x0 <- gsub("\\s+", "", x0)
  x0 <- ifelse(grepl("^S[0-9]+_", x0, ignore.case = TRUE), sub("_", ".", x0, fixed = TRUE), x0)
  is_full <- grepl("^S[0-9]+\\.", x0, ignore.case = TRUE)
  x_upper <- toupper(gsub("[^A-Za-z0-9]", "", x0))
  alias_hit <- mineral_code_alias[x_upper]
  out <- ifelse(is_full, x0, ifelse(!is.na(alias_hit), alias_hit, x0))
  as.character(out)
}

valid_growth <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x) | x < 0] <- 1
  x
}

force_rrnn_vector <- function(x, rrnn, nn, name = "object") {
  if (is.null(x)) stop(sprintf("%s is empty", name))

  if (is.matrix(x) || is.data.frame(x) || is.data.table(x)) {
    x <- as.matrix(x)
    if (nrow(x) == 1 && ncol(x) == rrnn) return(as.numeric(x[1, ]))
    if (ncol(x) == 1 && nrow(x) == rrnn) return(as.numeric(x[, 1]))
    if (length(x) == rrnn) return(as.numeric(x))
    stop(sprintf("%s has dimension %s; expected length %d", name, paste(dim(x), collapse = " x "), rrnn))
  }

  x <- as.numeric(x)
  if (length(x) != rrnn) stop(sprintf("%s has length %d; expected %d", name, length(x), rrnn))
  x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
  x
}

safe_assign_row <- function(mat, t, vec, mat_name, vec_name) {
  if (length(vec) != ncol(mat)) {
    stop(sprintf("%s row %d assignment mismatch: %s length=%d; %s columns=%d", mat_name, t, vec_name, length(vec), mat_name, ncol(mat)))
  }
  mat[t, ] <- vec
  mat
}

# --- Input-output aggregation ---
raw_region <- parse_region_dot(labels_raw)
raw_full_sector <- parse_sector_after_first_dot(labels_raw)
raw_last_token <- parse_last_token(labels_raw)
raw_model_sector <- normalize_to_model_sector(raw_full_sector, raw_last_token)

RegionList_vec <- unique(raw_region)
RR <- length(RegionList_vec)

SectorList_vec <- c(keep_sectors, "OTHER")
NN <- length(SectorList_vec)
RRNN <- RR * NN

region_id <- match(raw_region, RegionList_vec)
sector_id <- match(raw_model_sector, SectorList_vec)
sector_id[is.na(sector_id)] <- NN

new_node_id <- (region_id - 1) * NN + sector_id

G_node <- sparseMatrix(i = new_node_id, j = seq_along(labels_raw), x = 1, dims = c(RRNN, length(labels_raw)))

IOZ_0 <- as.matrix(G_node %*% Z_raw %*% t(G_node))
IOF_0 <- as.matrix(G_node %*% Y_raw)
IOV_0 <- as.matrix(V_raw %*% t(G_node))
IOTQ_0 <- as.matrix(TQ_raw %*% t(G_node))

rm(Z_raw, Y_raw, V_raw, TQ_raw, G_node); gc()

new_labels <- paste(rep(RegionList_vec, each = NN), rep(SectorList_vec, times = RR), sep = ".")
rownames(IOZ_0) <- new_labels
colnames(IOZ_0) <- new_labels

IOX_0_vec <- rowSums(IOZ_0, na.rm = TRUE) + rowSums(IOF_0, na.rm = TRUE)
IOX_0_vec_safe <- IOX_0_vec
IOX_0_vec_safe[is.na(IOX_0_vec_safe) | is.nan(IOX_0_vec_safe) | is.infinite(IOX_0_vec_safe) | IOX_0_vec_safe <= 0] <- 1e-12

I_Sum_sparse <- Matrix::Matrix(matrix(rep(diag(NN), RR), nrow = NN), sparse = TRUE)

IOZ_C <- as.matrix(I_Sum_sparse %*% IOZ_0)
IOF_C <- as.matrix(I_Sum_sparse %*% IOF_0)

E_CZ <- sweep(IOZ_C, 2, IOX_0_vec_safe, "/")
E_CZ[is.na(E_CZ) | is.nan(E_CZ) | is.infinite(E_CZ)] <- 0

E_VA_vec_raw <- as.numeric(IOV_0) / IOX_0_vec_safe
E_VA_vec_raw[is.na(E_VA_vec_raw) | is.nan(E_VA_vec_raw) | is.infinite(E_VA_vec_raw)] <- 0

E_TQ_vec <- as.numeric(IOTQ_0) / IOX_0_vec_safe
E_TQ_vec[is.na(E_TQ_vec) | is.nan(E_TQ_vec) | is.infinite(E_TQ_vec)] <- 0

mat.key <- E_CZ <= 0
mat.key[is.na(mat.key)] <- FALSE

for (bat in names(critical_mineral_to_battery)) {
  bat_local <- which(SectorList_vec == bat)
  min_local <- which(SectorList_vec %in% critical_mineral_to_battery[[bat]])
  if (length(bat_local) == 0 || length(min_local) == 0) next
  bat_global_cols_tmp <- unlist(lapply(0:(RR - 1), function(r) r * NN + bat_local))
  mat.key[min_local, bat_global_cols_tmp] <- FALSE
}

StockPara_0 <- matrix(STOCK_PARA_VALUE, nrow = NN, ncol = RRNN)

battery_local_idx <- which(SectorList_vec %in% battery_sectors)
battery_global_cols <- unlist(lapply(0:(RR - 1), function(r) r * NN + battery_local_idx))
mineral_local_idx <- which(SectorList_vec %in% mineral_sectors)

if (any(IOX_0_vec[battery_global_cols] <= 0)) stop("Base-year storage output must be positive")
battery_node_dt <- data.table(
  Node_Index = battery_global_cols,
  Node_Name = new_labels[battery_global_cols],
  Region = parse_region_dot(new_labels[battery_global_cols]),
  Sector = extract_sector(new_labels[battery_global_cols]),
  IOX_0 = IOX_0_vec[battery_global_cols],
  IOV_0 = as.numeric(IOV_0)[battery_global_cols],
  E_VA_raw = E_VA_vec_raw[battery_global_cols]
)

fwrite(battery_node_dt, file.path(out_root, "Node_Reference.csv"))

log_line(sprintf("Aggregation complete: RR=%d, NN=%d, RRNN=%d, BatteryNodes=%d", RR, NN, RRNN, length(battery_global_cols)))

# --- Storage demand projections ---
log_line("Loading demand projections")

battery_data_global <- fread(battery_forecast_path)

required_cols <- c("Scenario", "Technology", "Storage_Type", "Year", "Shipment_GWh")
missing_cols <- setdiff(required_cols, names(battery_data_global))
if (length(missing_cols) > 0) stop(sprintf("Missing demand columns: %s", paste(missing_cols, collapse = ", ")))

battery_data_global[, Macro_raw := trimws(as.character(Scenario))]
battery_data_global[, Macro_norm := normalize_macro_scenario(Macro_raw)]
battery_data_global[, Tech_raw := trimws(as.character(Technology))]
battery_data_global[, Tech_norm := normalize_tech_scenario(Tech_raw)]
battery_data_global[, Battery_Class_raw := trimws(as.character(Storage_Type))]
battery_data_global[, Battery_Code := Battery_Class_raw]
if (any(!battery_data_global$Battery_Code %in% battery_sectors)) {
  stop("Unknown storage technology code")
}
battery_data_global[, Year := as.integer(Year)]
battery_data_global[, Shipment_GWH := suppressWarnings(as.numeric(Shipment_GWh))]
if (any(!is.finite(battery_data_global$Shipment_GWH)) || any(battery_data_global$Shipment_GWH < 0)) {
  stop("Storage demand contains invalid values")
}
if (anyDuplicated(battery_data_global, by = c("Scenario", "Technology", "Storage_Type", "Year"))) {
  stop("Duplicate demand scenario, technology, or year records")
}

fwrite(unique(battery_data_global[, .(Macro_raw, Macro_norm, Tech_raw, Tech_norm, Battery_Class_raw, Battery_Code)]),
       file.path(out_root, "Demand_Mapping.csv"))

make_battery_demand_matrix <- function(macro_raw, macro_norm, tech_raw, tech_norm) {
  batt_scene_data <- battery_data_global[Macro_norm == macro_norm & Tech_norm == tech_norm]

  if (nrow(batt_scene_data) == 0) {
    msg <- sprintf("Demand matching failed: macro='%s', normalized='%s'; technology='%s', normalized='%s'", macro_raw, macro_norm, tech_raw, tech_norm)
    if (STOP_IF_BATTERY_FORECAST_UNMATCHED) stop(msg)
    log_line(paste("Warning:", msg))
  }

  Battery_Demand_TT <- matrix(0, nrow = TT, ncol = RRNN)
  colnames(Battery_Demand_TT) <- new_labels
  check_list <- list()

  for (bat_code in battery_sectors) {
    bat_local <- which(SectorList_vec == bat_code)
    if (length(bat_local) == 0) next

    bat_nodes <- unlist(lapply(0:(RR - 1), function(r) r * NN + bat_local))
    sub <- batt_scene_data[Battery_Code == bat_code]

    if (nrow(sub) == 0) {
      multiplier <- rep(1, TT)
      shipment_vec <- rep(NA_real_, TT)
      names(shipment_vec) <- years_cols
      matched <- FALSE
    } else {
      shipment_by_year <- sub[, .(Shipment_GWH = sum(Shipment_GWH, na.rm = TRUE)), by = Year][order(Year)]
      shipment_vec <- rep(NA_real_, TT)
      names(shipment_vec) <- years_cols

      for (yy in shipment_by_year$Year) {
        if (as.character(yy) %in% years_cols) {
          shipment_vec[as.character(yy)] <- shipment_by_year[Year == yy, Shipment_GWH][1]
        }
      }

      if (all(is.na(shipment_vec))) {
        multiplier <- rep(1, TT)
      } else {
        for (k in seq_along(shipment_vec)) {
          if (is.na(shipment_vec[k])) {
            prev <- which(!is.na(shipment_vec[1:k]))
            if (length(prev) > 0) shipment_vec[k] <- shipment_vec[max(prev)]
          }
        }
        for (k in rev(seq_along(shipment_vec))) {
          if (is.na(shipment_vec[k])) {
            nxt <- which(!is.na(shipment_vec[k:length(shipment_vec)])) + k - 1
            if (length(nxt) > 0) shipment_vec[k] <- shipment_vec[min(nxt)]
          }
        }
        shipment_vec[is.na(shipment_vec)] <- 0

        base_2024 <- shipment_vec[as.character(2024)]
        if (is.na(base_2024) || base_2024 <= 0) {
          positive_vals <- shipment_vec[shipment_vec > 0]
          base_2024 <- ifelse(length(positive_vals) > 0, positive_vals[1], 1)
        }

        multiplier <- shipment_vec / base_2024
        multiplier[is.na(multiplier) | is.nan(multiplier) | is.infinite(multiplier)] <- 1
        multiplier <- pmax(multiplier, 0)
      }
      matched <- TRUE
    }

    base_econ <- IOX_0_vec[bat_nodes]
    total_base_econ <- sum(base_econ, na.rm = TRUE)

    if (USE_GWH_DIRECTLY_AS_DEMAND) {
      demand_total_by_year <- shipment_vec
      demand_total_by_year[is.na(demand_total_by_year)] <- 0
      shares <- if (total_base_econ > 0) base_econ / total_base_econ else rep(1 / length(bat_nodes), length(bat_nodes))
      demand_mat <- outer(demand_total_by_year, shares)
    } else {
      demand_mat <- outer(multiplier, base_econ)
    }

    demand_mat[is.na(demand_mat) | is.nan(demand_mat) | is.infinite(demand_mat)] <- 0
    Battery_Demand_TT[, bat_nodes] <- demand_mat

    check_list[[length(check_list) + 1]] <- data.table(
      Battery_Code = bat_code,
      Battery_Class_Matched = paste(sort(unique(na.omit(sub$Battery_Class_raw))), collapse = "|"),
      Matched = matched,
      Matched_Row_N = nrow(sub),
      Node_N = length(bat_nodes),
      Base_IOX0_Total = total_base_econ,
      Shipment_2024 = shipment_vec[as.character(2024)],
      Shipment_2030 = shipment_vec[as.character(2030)],
      Shipment_2050 = shipment_vec[as.character(2050)],
      Shipment_2100 = shipment_vec[as.character(2100)],
      Multiplier_2030 = multiplier[which(years == 2030)],
      Multiplier_2050 = multiplier[which(years == 2050)],
      Multiplier_2100 = multiplier[which(years == 2100)],
      Multiplier_Unique_N = uniqueN(round(multiplier, 12)),
      Demand_2024_Total = sum(Battery_Demand_TT[which(years == 2024), bat_nodes], na.rm = TRUE),
      Demand_2030_Total = sum(Battery_Demand_TT[which(years == 2030), bat_nodes], na.rm = TRUE),
      Demand_2050_Total = sum(Battery_Demand_TT[which(years == 2050), bat_nodes], na.rm = TRUE),
      Demand_2100_Total = sum(Battery_Demand_TT[which(years == 2100), bat_nodes], na.rm = TRUE)
    )
  }

  Battery_Demand_Check <- rbindlist(check_list, fill = TRUE)
  Storage_Demand_TS <- data.table(
    Year = years,
    Storage_Demand_Total = rowSums(Battery_Demand_TT[, battery_global_cols, drop = FALSE], na.rm = TRUE)
  )

  if (nrow(Battery_Demand_Check) == 0 || any(!Battery_Demand_Check$Matched)) {
    stop(sprintf(
      "Scenario %s / %s is missing a storage technology demand series",
      macro_raw, tech_raw
    ))
  }

  list(
    Battery_Demand_TT = Battery_Demand_TT,
    Battery_Demand_Check = Battery_Demand_Check,
    Storage_Demand_TS = Storage_Demand_TS
  )
}

# --- Production and mineral supply paths ---
scenario_table <- fread(file.path(project_root, "config", "scenario_manifest.csv"))
required <- c("Scenario_ID", "Production_File")
if (!all(required %in% names(scenario_table))) stop("Missing scenario manifest columns")
if (anyDuplicated(scenario_table$Scenario_ID)) stop("Duplicate scenario IDs")
if (nrow(scenario_table) != 3L) stop("Three macro scenarios are required")
scenario_table[, Macro_raw := Scenario_ID]
scenario_table[, Macro_norm := normalize_macro_scenario(Macro_raw)]
scenario_table[, Scenario_File := file.path(main_scenario_path, Production_File)]
if (any(!file.exists(scenario_table$Scenario_File))) stop("Missing production scenario workbook")
prod_list <- list()
for (i in seq_len(nrow(scenario_table))) {
  nm <- scenario_table$Macro_norm[i]
  prod_list[[nm]] <- as.data.table(read_excel(scenario_table$Scenario_File[i], sheet = "Production"))
}
fwrite(scenario_table, file.path(out_root, "scenario_index.csv"))

prepare_production_wide <- function(prod_df) {
  prod_df <- as.data.table(copy(prod_df))
  names(prod_df) <- trimws(as.character(names(prod_df)))

  if (!all(c("Region", "mapped_code") %in% names(prod_df))) {
    return(data.table())
  }

  prod_df[, Region := trimws(as.character(Region))]
  prod_df[, mapped_code_raw := trimws(as.character(mapped_code))]
  prod_df[, mapped_code_std := standardize_mapped_code(mapped_code_raw)]
  prod_df[, full_label_raw := paste(Region, mapped_code_raw, sep = ".")]
  prod_df[, full_label_std := paste(Region, mapped_code_std, sep = ".")]

  for (yc in years_cols) {
    if (!yc %in% names(prod_df)) {
      prod_df[, (yc) := NA_real_]
    } else {
      prod_df[, (yc) := suppressWarnings(as.numeric(get(yc)))]
    }
  }

  prod_df
}

sum_production_by_year <- function(rows_dt) {
  vals <- rep(NA_real_, TT)
  names(vals) <- years_cols
  if (is.null(rows_dt) || nrow(rows_dt) == 0) return(vals)

  for (yc in years_cols) {
    if (!yc %in% names(rows_dt)) {
      vals[yc] <- NA_real_
    } else {
      v <- suppressWarnings(as.numeric(rows_dt[[yc]]))
      if (length(v) == 0 || all(is.na(v))) {
        vals[yc] <- NA_real_
      } else {
        vals[yc] <- sum(v, na.rm = TRUE)
      }
    }
  }
  vals
}

fill_growth_locf <- function(g) {
  g <- as.numeric(g)
  names(g) <- years_cols

  g[is.nan(g) | is.infinite(g) | g < 0] <- NA_real_

  if (is.na(g["2024"]) || !is.finite(g["2024"])) g["2024"] <- 1

  for (i in seq_along(g)) {
    if (is.na(g[i]) || !is.finite(g[i])) {
      if (i == 1) {
        g[i] <- 1
      } else {
        g[i] <- g[i - 1]
      }
    }
  }

  g[is.na(g) | is.nan(g) | is.infinite(g)] <- 1
  g
}

make_growth_from_rows <- function(rows_dt, fallback_growth = NULL, apply_strength = FALSE) {
  vals <- sum_production_by_year(rows_dt)
  base <- vals["2024"]

  if (is.finite(base) && !is.na(base) && base > 0) {
    raw_g <- vals / base
    raw_g[is.nan(raw_g) | is.infinite(raw_g)] <- NA_real_
    g <- fill_growth_locf(raw_g)
  } else {
    if (!is.null(fallback_growth)) {
      g <- fill_growth_locf(fallback_growth)
    } else {
      g <- rep(1, TT)
      names(g) <- years_cols
    }
  }

  if (apply_strength) {
    g <- valid_growth(g)
  } else {
    g[is.na(g) | is.nan(g) | is.infinite(g) | g < 0] <- 1
  }

  names(g) <- years_cols
  g
}

make_macro_growth_matrix <- function(macro_norm) {
  growth_matrix <- matrix(1, nrow = RRNN, ncol = TT)
  mineral_supply_growth_matrix <- matrix(1, nrow = RRNN, ncol = TT)
  colnames(growth_matrix) <- years_cols
  colnames(mineral_supply_growth_matrix) <- years_cols
  rownames(growth_matrix) <- new_labels
  rownames(mineral_supply_growth_matrix) <- new_labels

  empty_check <- data.table(
    Mineral_Code = mineral_sectors,
    Direct_Region_N = 0L,
    Fallback_Region_N = RR,
    Missing_Model_Node_N = 0L,
    Global_Growth_2025 = 1,
    Global_Growth_2030 = 1,
    Global_Growth_2050 = 1,
    Global_Growth_2100 = 1,
    Last_Year_Growth = 1,
    Available_Year_Max = NA_integer_
  )

  if (!macro_norm %in% names(prod_list)) {
    log_line(sprintf("Warning: macro scenario %s not found; default growth is 1", macro_norm))
    return(list(
      growth_matrix = growth_matrix,
      mineral_supply_growth_matrix = mineral_supply_growth_matrix,
      mineral_growth_check = empty_check
    ))
  }

  prod_df <- prepare_production_wide(prod_list[[macro_norm]])

  if (nrow(prod_df) == 0) {
    log_line(sprintf("Warning: Production sheet lacks Region/mapped_code or is empty: %s", macro_norm))
    return(list(
      growth_matrix = growth_matrix,
      mineral_supply_growth_matrix = mineral_supply_growth_matrix,
      mineral_growth_check = empty_check
    ))
  }

  year_has_data <- sapply(years_cols, function(yc) {
    v <- suppressWarnings(as.numeric(prod_df[[yc]]))
    any(is.finite(v) & !is.na(v))
  })
  available_years <- as.integer(years_cols[year_has_data])
  available_year_max <- if (length(available_years) > 0) max(available_years, na.rm = TRUE) else NA_integer_

  unique_labels <- unique(prod_df$full_label_std)
  unique_labels <- unique_labels[!is.na(unique_labels) & unique_labels != ""]

  for (fl in unique_labels) {
    idx <- match(fl, new_labels)
    if (is.na(idx)) next
    rows_fl <- prod_df[full_label_std == fl]
    growth_matrix[idx, ] <- make_growth_from_rows(rows_fl, fallback_growth = rep(1, TT), apply_strength = FALSE)
  }

  mineral_check_list <- list()

  for (ms in mineral_sectors) {
    rows_ms <- prod_df[mapped_code_std == ms]

    global_growth <- make_growth_from_rows(rows_ms, fallback_growth = rep(1, TT), apply_strength = TRUE)

    direct_region_n <- 0L
    fallback_region_n <- 0L
    missing_model_node_n <- 0L

    for (reg in RegionList_vec) {
      target_idx <- match(paste(reg, ms, sep = "."), new_labels)
      if (is.na(target_idx)) {
        missing_model_node_n <- missing_model_node_n + 1L
        next
      }

      rows_reg_ms <- rows_ms[Region == reg]
      vals_reg <- sum_production_by_year(rows_reg_ms)
      base_reg <- vals_reg["2024"]

      if (nrow(rows_reg_ms) > 0 && is.finite(base_reg) && !is.na(base_reg) && base_reg > 0) {
        reg_growth <- make_growth_from_rows(rows_reg_ms, fallback_growth = global_growth, apply_strength = TRUE)
        mineral_supply_growth_matrix[target_idx, ] <- reg_growth
        direct_region_n <- direct_region_n + 1L
      } else {
        mineral_supply_growth_matrix[target_idx, ] <- global_growth
        fallback_region_n <- fallback_region_n + 1L
      }
    }

    mineral_check_list[[length(mineral_check_list) + 1]] <- data.table(
      Mineral_Code = ms,
      Direct_Region_N = direct_region_n,
      Fallback_Region_N = fallback_region_n,
      Missing_Model_Node_N = missing_model_node_n,
      Global_Growth_2025 = global_growth[as.character(2025)],
      Global_Growth_2030 = global_growth[as.character(2030)],
      Global_Growth_2050 = global_growth[as.character(2050)],
      Global_Growth_2100 = global_growth[as.character(2100)],
      Last_Year_Growth = global_growth[length(global_growth)],
      Available_Year_Max = available_year_max
    )
  }

  mineral_growth_check <- rbindlist(mineral_check_list, fill = TRUE)

  list(
    growth_matrix = growth_matrix,
    mineral_supply_growth_matrix = mineral_supply_growth_matrix,
    mineral_growth_check = mineral_growth_check
  )
}

make_mineral_supply_multiplier_matrix <- function(mineral_supply_growth_matrix, t) {
  mult <- matrix(1, nrow = NN, ncol = RRNN)

  dest_regions <- parse_region_dot(new_labels)

  for (ms in mineral_sectors) {
    row_id <- which(SectorList_vec == ms)
    if (length(row_id) == 0) next

    idx_by_dest_region <- match(paste(dest_regions, ms, sep = "."), new_labels)
    m <- rep(1, RRNN)
    valid <- !is.na(idx_by_dest_region)
    m[valid] <- mineral_supply_growth_matrix[idx_by_dest_region[valid], t]
    m[!is.finite(m) | m < 0] <- 0
    mult[row_id, ] <- m
  }

  mult
}

# --- Technology coefficients ---
TECH_MINERAL_INTENSITY_DT <- fread(file.path(project_root, "config", "technology_coefficients.csv"))
required <- c("Tech_Code", "Battery_Code", "Mineral_Code", "Multiplier")
if (!all(required %in% names(TECH_MINERAL_INTENSITY_DT))) stop("Missing technology coefficient columns")
if (anyDuplicated(TECH_MINERAL_INTENSITY_DT, by = c("Tech_Code", "Battery_Code", "Mineral_Code"))) stop("Duplicate technology coefficients")
if (any(!is.finite(TECH_MINERAL_INTENSITY_DT$Multiplier) | TECH_MINERAL_INTENSITY_DT$Multiplier <= 0)) stop("Technology multipliers must be positive")
if (!all(TECH_MINERAL_INTENSITY_DT$Battery_Code %in% battery_sectors)) stop("Invalid battery sector code")
if (!all(TECH_MINERAL_INTENSITY_DT$Mineral_Code %in% mineral_sectors)) stop("Invalid mineral sector code")
tech_table <- fread(file.path(project_root, "config", "technology_scenarios.csv"))
if (!all(c("Tech_Code", "Technology") %in% names(tech_table))) stop("Missing technology scenario columns")
if (nrow(tech_table) != 4L || anyDuplicated(tech_table$Tech_Code)) stop("Four unique technology scenarios are required")
tech_scenarios <- tech_table$Technology
expected_demand <- CJ(Scenario = scenario_table$Scenario_ID, Technology = tech_scenarios,
                     Storage_Type = battery_sectors, Year = years)
actual_demand <- battery_data_global[, .(Scenario, Technology, Storage_Type, Year)]
if (nrow(actual_demand) != nrow(expected_demand)) stop("Incorrect number of storage demand records")
if (nrow(fsetdiff(expected_demand, actual_demand)) || nrow(fsetdiff(actual_demand, expected_demand))) stop("Incomplete storage demand scenario combinations")
if (any(!TECH_MINERAL_INTENSITY_DT$Tech_Code %in% tech_table$Tech_Code)) stop("Unknown technology scenario in coefficient table")
for (tc in tech_table$Tech_Code[tech_table$Tech_Code != "S0"]) {
  if (!nrow(TECH_MINERAL_INTENSITY_DT[Tech_Code == tc])) stop("Missing mineral coefficients for technology scenario: ", tc)
}

make_tech_mineral_intensity_matrix <- function(tech_raw) {
  if (get_tech_code(tech_raw) != "S0") stop("Technology model is not initialized")
  matrix(1, nrow = NN, ncol = RRNN)
}

# --- Inventory and sector balance ---
sync_order_rows_to_target <- function(OrderMat, target_rows, target_vec, final_demand_start_col) {
  if (length(target_rows) == 0) return(OrderMat)
  target_vec <- as.numeric(target_vec)
  target_vec[is.na(target_vec) | is.nan(target_vec) | is.infinite(target_vec) | target_vec < 0] <- 0
  old_sum <- rowSums(OrderMat[target_rows, , drop = FALSE], na.rm = TRUE)
  old_sum[is.na(old_sum) | is.nan(old_sum) | is.infinite(old_sum) | old_sum < 0] <- 0

  for (ii in seq_along(target_rows)) {
    rr <- target_rows[ii]
    target <- target_vec[ii]
    if (!is.finite(target) || target < 0) target <- 0

    if (old_sum[ii] > 0) {
      OrderMat[rr, ] <- OrderMat[rr, ] * (target / old_sum[ii])
    } else {
      OrderMat[rr, ] <- 0
      OrderMat[rr, final_demand_start_col] <- target
    }
  }
  OrderMat[is.na(OrderMat) | is.nan(OrderMat) | is.infinite(OrderMat)] <- 0
  OrderMat
}

make_et_sector_annual_check <- function(IOX_TT_Res, Demand_TT_Res, Gap_TT_Res, labels, years, sector_vec) {
  sec <- extract_sector(labels)
  out <- rbindlist(lapply(sector_vec, function(s) {
    cols <- which(sec == s)
    if (length(cols) == 0) {
      return(data.table(Year = years, Sector = s, Node_N = 0L, Output = 0, Demand = 0, Gap = 0))
    }
    data.table(
      Year = years,
      Sector = s,
      Node_N = length(cols),
      Output = rowSums(IOX_TT_Res[, cols, drop = FALSE], na.rm = TRUE),
      Demand = rowSums(Demand_TT_Res[, cols, drop = FALSE], na.rm = TRUE),
      Gap = rowSums(Gap_TT_Res[, cols, drop = FALSE], na.rm = TRUE)
    )
  }), fill = TRUE)
  out[, Sector_Type := ifelse(grepl("^ELEC_", Sector), "Electrified", "Traditional")]
  out[]
}

# --- Climate risk adjustments ---
build_risk_effect <- function(risk_profile, capacity_growth, supply_growth) {
  needed <- c("Year", "Region", "Mineral_Code", "Loss_Rate")
  if (!all(needed %in% names(risk_profile))) stop("Risk input columns are incomplete")
  shock <- data.table::as.data.table(data.table::copy(risk_profile))
  if (anyDuplicated(shock, by = c("Year", "Region", "Mineral_Code"))) {
    stop("Duplicate risk input keys")
  }
  if (any(!is.finite(shock$Loss_Rate)) ||
      any(shock$Loss_Rate < 0 | shock$Loss_Rate >= 1)) stop("Risk loss rates must be in [0,1)")
  if (any(!shock$Region %in% RegionList_vec)) stop("Unknown risk region")
  if (any(!shock$Mineral_Code %in% mineral_sectors)) stop("Unknown risk mineral")
  if (any(!shock$Year %in% years)) stop("Risk event year outside model period")
  if (nrow(shock) == 0) stop("Risk input is empty")
  shock[, Node := paste(Region, Mineral_Code, sep = ".")]
  shock[, Node_ID := match(Node, new_labels)]
  shock[, Year_ID := match(Year, years)]
  if (anyNA(shock$Node_ID) || anyNA(shock$Year_ID)) stop("Risk node mapping failed")
  shock[, `:=`(
    Base_Capacity_Growth = capacity_growth[cbind(Node_ID, Year_ID)],
    Base_Supply_Growth = supply_growth[cbind(Node_ID, Year_ID)]
  )]
  for (i in seq_len(nrow(shock))) {
    j <- shock$Node_ID[i]
    t <- shock$Year_ID[i]
    factor <- 1 - shock$Loss_Rate[i]
    capacity_growth[j, t] <- capacity_growth[j, t] * factor
    supply_growth[j, t] <- supply_growth[j, t] * factor
  }
  shock[, `:=`(
    Adjusted_Capacity_Growth = capacity_growth[cbind(Node_ID, Year_ID)],
    Adjusted_Supply_Growth = supply_growth[cbind(Node_ID, Year_ID)]
  )]
  list(capacity_growth = capacity_growth, supply_growth = supply_growth,
       diagnostics = shock[, .(Year, Region, Mineral_Code, Loss_Rate,
                                Base_Capacity_Growth, Adjusted_Capacity_Growth,
                                Base_Supply_Growth, Adjusted_Supply_Growth)])
}

# --- ARIO simulation ---
run_single_task <- function(macro_raw, macro_norm, scenario_file, tech_raw,
                            model_mode = "baseline", risk_profile = NULL) {
  if (!model_mode %in% c("baseline", "climate_risk", "technology")) {
    stop("Unknown model mode: ", model_mode)
  }
  if (model_mode == "climate_risk" && is.null(risk_profile)) {
    stop("Climate risk model requires an annual loss input")
  }
  tech_norm <- normalize_tech_scenario(tech_raw)
  out_dir <- file.path(out_root, model_mode, "scenarios",
                       safe_clean_filename(macro_raw), safe_clean_filename(get_tech_code(tech_raw)))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  append_task_status(macro_raw, tech_raw, "started", "", model_mode)
  log_line(sprintf("Starting: %s | %s", macro_raw, tech_raw))

  tryCatch({
    mg <- make_macro_growth_matrix(macro_norm)
    growth_matrix <- mg$growth_matrix
    mineral_supply_growth_matrix <- mg$mineral_supply_growth_matrix
    if (model_mode == "climate_risk") {
      impacts <- build_risk_effect(risk_profile, growth_matrix, mineral_supply_growth_matrix)
      growth_matrix <- impacts$capacity_growth
      mineral_supply_growth_matrix <- impacts$supply_growth
      fwrite(impacts$diagnostics, file.path(out_dir, "Risk_Applied.csv"))
    }
    fwrite(mg$mineral_growth_check, file.path(out_dir, "Mineral_Growth.csv"))
    save(mineral_supply_growth_matrix, file = file.path(out_dir, "Mineral_Supply_Growth_TT_Results.RData"))

    bd <- make_battery_demand_matrix(macro_raw, macro_norm, tech_raw, tech_norm)
    Battery_Demand_TT <- bd$Battery_Demand_TT
    fwrite(bd$Battery_Demand_Check, file.path(out_dir, "Storage_Demand_Check.csv"))
    fwrite(bd$Storage_Demand_TS, file.path(out_dir, "Storage_Demand_Series.csv"))
    save(Battery_Demand_TT, file = file.path(out_dir, "Battery_Demand_TT_Results.RData"))

    for (bat in battery_sectors) {
      bat_local <- which(SectorList_vec == bat)
      if (length(bat_local) == 0) next
      bat_nodes <- unlist(lapply(0:(RR - 1), function(r) r * NN + bat_local))
      total_base <- sum(Battery_Demand_TT[1, bat_nodes], na.rm = TRUE)
      if (total_base > 0) {
        mult <- rowSums(Battery_Demand_TT[, bat_nodes, drop = FALSE], na.rm = TRUE) / total_base
        mult[is.na(mult) | is.nan(mult) | is.infinite(mult)] <- 1
        mult <- pmax(mult, 0)
        growth_matrix[bat_nodes, ] <- matrix(mult, nrow = length(bat_nodes), ncol = TT, byrow = TRUE)
      }
    }

    IOX_TT_Res <- matrix(0, nrow = TT, ncol = RRNN)
    Demand_TT_Res <- matrix(0, nrow = TT, ncol = RRNN)
    Demand_Output_Gap_TT_Res <- matrix(0, nrow = TT, ncol = RRNN)
    Carbon_Footprint_Res <- matrix(0, nrow = TT, ncol = RRNN)
    MinZX_TT_Res <- matrix(0, nrow = TT, ncol = RRNN)
    VX_TT_Res <- matrix(0, nrow = TT, ncol = RRNN)
    OX_TT_Res <- matrix(0, nrow = TT, ncol = RRNN)
    Capacity_TT_Res <- matrix(0, nrow = TT, ncol = RRNN)

    colnames(IOX_TT_Res) <- new_labels
    colnames(Demand_TT_Res) <- new_labels
    colnames(Demand_Output_Gap_TT_Res) <- new_labels
    colnames(Carbon_Footprint_Res) <- new_labels
    colnames(MinZX_TT_Res) <- new_labels
    colnames(VX_TT_Res) <- new_labels
    colnames(OX_TT_Res) <- new_labels
    colnames(Capacity_TT_Res) <- new_labels

    StockMat <- StockPara_0 * IOZ_C
    StockAdd <- IOZ_C

    tech_intensity_mat <- make_tech_mineral_intensity_matrix(tech_raw)
    E_CZ_eff <- E_CZ * tech_intensity_mat
    E_CZ_eff[is.na(E_CZ_eff) | is.nan(E_CZ_eff) | is.infinite(E_CZ_eff) | E_CZ_eff < 0] <- 0
    fwrite(TECH_MINERAL_INTENSITY_DT[Tech_Code == get_tech_code(tech_raw)], file.path(out_dir, "Technology_Coefficients.csv"))

    IOV_0_vec <- as.numeric(IOV_0)
    IOV_0_vec[is.na(IOV_0_vec) | is.nan(IOV_0_vec) | is.infinite(IOV_0_vec)] <- 0

    E_VA_vec <- E_VA_vec_raw
    E_VA_vec[E_VA_vec <= 0 | is.na(E_VA_vec) | is.nan(E_VA_vec) | is.infinite(E_VA_vec)] <- 1e-12

    Labor_vec <- rep(1, RRNN)

    for (t in seq_len(TT)) {
      current_capacity_max <- CAPACITY_FACTOR * IOX_0_vec * growth_matrix[, t]
      current_capacity_max <- force_rrnn_vector(current_capacity_max, RRNN, NN, "current_capacity_max")
      current_capacity_max[current_capacity_max < 0] <- 0

      capacity_mat <- matrix(rep(current_capacity_max, each = NN), nrow = NN, ncol = RRNN)

      final_demand_t <- sweep(IOF_0, 1, growth_matrix[, t], "*")
      final_demand_t[is.na(final_demand_t) | is.nan(final_demand_t) | is.infinite(final_demand_t)] <- 0

      OrderMat <- cbind(IOZ_0, final_demand_t)

      if (length(battery_global_cols) > 0) {
        OrderMat <- sync_order_rows_to_target(
          OrderMat = OrderMat,
          target_rows = battery_global_cols,
          target_vec = Battery_Demand_TT[t, battery_global_cols],
          final_demand_start_col = ncol(IOZ_0) + 1L
        )
      }

      OX <- rowSums(OrderMat, na.rm = TRUE)
      OX <- force_rrnn_vector(OX, RRNN, NN, "OX_after_order_sync")
      OX[OX < 0] <- 0

      OX[battery_global_cols] <- Battery_Demand_TT[t, battery_global_cols]
      OX <- force_rrnn_vector(OX, RRNN, NN, "OX_after_battery_safeguard")

      IOV_t_vec <- IOV_0_vec * Labor_vec * growth_matrix[, t]
      IOV_t_vec[is.na(IOV_t_vec) | is.nan(IOV_t_vec) | is.infinite(IOV_t_vec) | IOV_t_vec < 0] <- 0

      VX <- IOV_t_vec / E_VA_vec
      zero_va <- E_VA_vec_raw <= 0 | !is.finite(E_VA_vec_raw)
      VX[zero_va] <- current_capacity_max[zero_va]
      VX[!is.finite(VX) | VX < 0] <- 0
      VX <- force_rrnn_vector(VX, RRNN, NN, "VX")

      if (ENABLE_MINERAL_SUPPLY_EXPANSION) {
        mineral_supply_mult <- make_mineral_supply_multiplier_matrix(mineral_supply_growth_matrix, t)
        mineral_supply_cap <- IOZ_C * mineral_supply_mult
        StockAdd[mineral_local_idx, ] <- pmin(
          StockAdd[mineral_local_idx, , drop = FALSE],
          mineral_supply_cap[mineral_local_idx, , drop = FALSE]
        )
      }

      temp_max_cons <- sweep(E_CZ_eff, 2, current_capacity_max, "*")
      temp_max_cons[is.na(temp_max_cons) | is.nan(temp_max_cons) | is.infinite(temp_max_cons) | temp_max_cons < 0] <- 0

      StockCons <- pmin(StockAdd + StockMat, temp_max_cons, na.rm = TRUE)
      StockCons[is.na(StockCons) | is.nan(StockCons) | is.infinite(StockCons)] <- 0

      ZX <- StockCons / E_CZ_eff
      ZX[is.na(ZX) | is.nan(ZX) | is.infinite(ZX)] <- Inf
      ZX[mat.key] <- capacity_mat[mat.key]

      for (bat in names(critical_mineral_to_battery)) {
        bat_local <- which(SectorList_vec == bat)
        critical_rows <- which(SectorList_vec %in% critical_mineral_to_battery[[bat]])
        if (length(bat_local) == 0 || length(critical_rows) == 0) next

        bat_cols <- unlist(lapply(0:(RR - 1), function(r) r * NN + bat_local))
        tmp_zx <- ZX[critical_rows, bat_cols, drop = FALSE]
        tmp_cap <- capacity_mat[critical_rows, bat_cols, drop = FALSE]
        tmp_ecz <- E_CZ_eff[critical_rows, bat_cols, drop = FALSE]
        zero_coef_mask <- tmp_ecz <= 0 | is.na(tmp_ecz) | is.nan(tmp_ecz) | is.infinite(tmp_ecz)
        tmp_zx[zero_coef_mask] <- tmp_cap[zero_coef_mask]
        ZX[critical_rows, bat_cols] <- tmp_zx
      }

      min_ZX <- matrixStats::colMins(ZX, na.rm = TRUE)
      invalid_stock <- !is.finite(min_ZX)
      min_ZX[invalid_stock] <- current_capacity_max[invalid_stock]
      min_ZX <- force_rrnn_vector(min_ZX, RRNN, NN, "min_ZX")

      IOX_t <- pmin(min_ZX, VX, OX, current_capacity_max, na.rm = TRUE)
      IOX_t[is.na(IOX_t) | is.nan(IOX_t) | is.infinite(IOX_t) | IOX_t < 0] <- 0
      IOX_t <- force_rrnn_vector(IOX_t, RRNN, NN, "IOX_t")

      Demand_gap <- pmax(OX - IOX_t, 0)
      Demand_gap <- force_rrnn_vector(Demand_gap, RRNN, NN, "Demand_gap")

      Carbon_t <- IOX_t * E_TQ_vec
      Carbon_t[is.na(Carbon_t) | is.nan(Carbon_t) | is.infinite(Carbon_t)] <- 0
      Carbon_t <- force_rrnn_vector(Carbon_t, RRNN, NN, "Carbon_t")

      IOX_TT_Res <- safe_assign_row(IOX_TT_Res, t, IOX_t, "IOX_TT_Res", "IOX_t")
      Demand_TT_Res <- safe_assign_row(Demand_TT_Res, t, OX, "Demand_TT_Res", "OX")
      Demand_Output_Gap_TT_Res <- safe_assign_row(Demand_Output_Gap_TT_Res, t, Demand_gap, "Demand_Output_Gap_TT_Res", "Demand_gap")
      Carbon_Footprint_Res <- safe_assign_row(Carbon_Footprint_Res, t, Carbon_t, "Carbon_Footprint_Res", "Carbon_t")
      MinZX_TT_Res <- safe_assign_row(MinZX_TT_Res, t, min_ZX, "MinZX_TT_Res", "min_ZX")
      VX_TT_Res <- safe_assign_row(VX_TT_Res, t, VX, "VX_TT_Res", "VX")
      OX_TT_Res <- safe_assign_row(OX_TT_Res, t, OX, "OX_TT_Res", "OX")
      Capacity_TT_Res <- safe_assign_row(Capacity_TT_Res, t, current_capacity_max, "Capacity_TT_Res", "current_capacity_max")

      actual_cons <- sweep(E_CZ_eff, 2, IOX_t, "*")
      actual_cons[is.na(actual_cons) | is.nan(actual_cons) | is.infinite(actual_cons) | actual_cons < 0] <- 0

      StockMat <- pmax(StockMat + StockAdd - actual_cons, 0)

      replenishment_driver <- pmax(IOX_t / IOX_0_vec_safe, 0)
      StockAdd <- sweep(IOZ_C, 2, replenishment_driver, "*")
      StockAdd[is.na(StockAdd) | is.nan(StockAdd) | is.infinite(StockAdd)] <- 0
    }

    Scenario_Metrics_Summary <- data.table(
      Node_Index = seq_len(RRNN),
      Node_Name = new_labels,
      Region = parse_region_dot(new_labels),
      Sector = extract_sector(new_labels),
      IOX_0 = IOX_0_vec,
      IOV_0 = as.numeric(IOV_0),
      E_VA_raw = E_VA_vec_raw,
      E_TQ = E_TQ_vec
    )

    Constraint_Diagnostics_TT_Res <- list(
      OX_TT_Res = OX_TT_Res,
      MinZX_TT_Res = MinZX_TT_Res,
      VX_TT_Res = VX_TT_Res,
      Capacity_TT_Res = Capacity_TT_Res,
      IOX_TT_Res = IOX_TT_Res
    )

    save(IOX_TT_Res, file = file.path(out_dir, "IOX_TT_Results.RData"))
    save(Demand_TT_Res, file = file.path(out_dir, "Demand_TT_Results.RData"))
    save(Demand_Output_Gap_TT_Res, file = file.path(out_dir, "Demand_Output_Gap_TT_Results.RData"))
    save(Carbon_Footprint_Res, file = file.path(out_dir, "Carbon_Footprint_Results.RData"))
    save(Constraint_Diagnostics_TT_Res, file = file.path(out_dir, "Constraint_Diagnostics_TT_Results.RData"))
    fwrite(Scenario_Metrics_Summary, file.path(out_dir, "Scenario_Metrics_Summary.csv"))

    storage_check <- data.table(
      Year = years,
      Storage_Demand = rowSums(Demand_TT_Res[, battery_global_cols, drop = FALSE], na.rm = TRUE),
      Storage_Output = rowSums(IOX_TT_Res[, battery_global_cols, drop = FALSE], na.rm = TRUE),
      Storage_Gap = rowSums(Demand_Output_Gap_TT_Res[, battery_global_cols, drop = FALSE], na.rm = TRUE),
      Battery_Demand_Exogenous = rowSums(Battery_Demand_TT[, battery_global_cols, drop = FALSE], na.rm = TRUE),
      Storage_MinZX = rowSums(MinZX_TT_Res[, battery_global_cols, drop = FALSE], na.rm = TRUE),
      Storage_VX = rowSums(VX_TT_Res[, battery_global_cols, drop = FALSE], na.rm = TRUE),
      Storage_Capacity = rowSums(Capacity_TT_Res[, battery_global_cols, drop = FALSE], na.rm = TRUE)
    )
    fwrite(storage_check, file.path(out_dir, "Storage_Balance.csv"))

    et_cols <- which(extract_sector(new_labels) %in% sector_export)
      if (length(et_cols) > 0) {
        Sector_IOX_TT_Res <- IOX_TT_Res[, et_cols, drop = FALSE]
        Sector_Demand_TT_Res <- Demand_TT_Res[, et_cols, drop = FALSE]
        Sector_Gap_TT_Res <- Demand_Output_Gap_TT_Res[, et_cols, drop = FALSE]
        save(Sector_IOX_TT_Res, Sector_Demand_TT_Res, Sector_Gap_TT_Res,
             file = file.path(out_dir, "Sector_Annual_Matrices.RData"))

        et_annual_check <- make_et_sector_annual_check(
          IOX_TT_Res = IOX_TT_Res,
          Demand_TT_Res = Demand_TT_Res,
          Gap_TT_Res = Demand_Output_Gap_TT_Res,
          labels = new_labels,
          years = years,
          sector_vec = sector_export
        )
        fwrite(et_annual_check, file.path(out_dir, "Sector_Annual_Check.csv"))
      }

    append_task_status(macro_raw, tech_raw, "success", "completed", model_mode)
    log_line(sprintf("Finished: %s | %s", macro_raw, tech_raw))
    TRUE
  }, error = function(e) {
    msg <- e$message
    append_task_status(macro_raw, tech_raw, "failed", msg, model_mode)
    log_line(sprintf("Run failed: %s | %s -> %s", macro_raw, tech_raw, msg))
    FALSE
  })
}


# Run one model family and check accounting identities.
run_model_family <- function(mode, risk_data = NULL) {
  if (!mode %in% c("baseline", "climate_risk", "technology")) stop("Invalid model family")
  family_dir <- file.path(out_root, mode)
  dir.create(family_dir, recursive = TRUE, showWarnings = FALSE)
  task_log_file <<- file.path(family_dir, "task_status.csv")
  error_log_file <<- file.path(family_dir, "run_log.txt")
  if (file.exists(task_log_file)) file.remove(task_log_file)
  if (file.exists(error_log_file)) file.remove(error_log_file)

  selected <- if (mode == "technology") tech_table[Tech_Code != "S0"] else tech_table[Tech_Code == "S0"]
  expected_runs <- nrow(scenario_table) * nrow(selected)
  required_runs <- if (mode == "technology") 9L else 3L
  if (expected_runs != required_runs) stop("Incorrect task count")
  result_list <- list()
  sector_list <- list()

  for (i in seq_len(nrow(scenario_table))) {
    for (j in seq_len(nrow(selected))) {
      macro <- scenario_table$Macro_raw[i]
      tech <- selected$Technology[j]
      good <- run_single_task(macro, scenario_table$Macro_norm[i],
                              scenario_table$Scenario_File[i], tech,
                              model_mode = mode, risk_profile = risk_data)
      if (!isTRUE(good)) stop(sprintf("Model simulation failed: %s / %s / %s", mode, macro, tech))
      scenario_dir <- file.path(family_dir, "scenarios", macro, selected$Tech_Code[j])
      balance <- fread(file.path(scenario_dir, "Storage_Balance.csv"))
      balance[, `:=`(Model = mode, Scenario = macro, Technology = selected$Tech_Code[j])]
      if (nrow(balance) != TT || any(balance$Year != years)) stop("Scenario annual coverage failed")
      err <- abs(balance$Storage_Demand - balance$Storage_Output - balance$Storage_Gap)
      if (any(!is.finite(err)) || any(err > 1e-7 * pmax(1, balance$Storage_Demand))) {
        stop("Storage accounting check failed: ", macro, "/", tech)
      }
      if (any(balance$Storage_Output < -1e-8) ||
          any(balance$Storage_Output > balance$Storage_Demand + 1e-7)) stop("Invalid production bounds")
      result_list[[length(result_list) + 1L]] <- balance
      sector <- fread(file.path(scenario_dir, "Sector_Annual_Check.csv"))
      sector[, `:=`(Model = mode, Scenario = macro, Technology = selected$Tech_Code[j])]
      sector_list[[length(sector_list) + 1L]] <- sector
    }
  }
  results <- rbindlist(result_list, fill = TRUE)
  setcolorder(results, c("Model", "Scenario", "Technology", "Year"))
  fwrite(results, file.path(family_dir, "scenario_summary.csv"))
  sector <- rbindlist(sector_list, fill = TRUE)
  setcolorder(sector, c("Model", "Scenario", "Technology", "Year", "Sector"))
  fwrite(sector, file.path(family_dir, "sector_balance.csv"))
  fwrite(results[Year %in% c(2024L, 2030L, 2050L, 2075L, 2100L),
                 .(Model, Scenario, Technology, Year, Storage_Demand, Storage_Output, Storage_Gap)],
         file.path(family_dir, "scenario_comparison.csv"))
  status <- fread(task_log_file)
  if (nrow(status[Status == "success"]) != expected_runs || any(status$Status == "failed")) {
    stop("Unexpected task status")
  }
  if (uniqueN(results, by = c("Scenario", "Technology")) != expected_runs) stop("Missing scenario combination")
  if (nrow(results) != expected_runs * TT) stop("Incomplete scenario output")
  cat(sprintf("PASS: %s (%d simulations, %d annual records)\n", mode, expected_runs, nrow(results)))
  invisible(results)
}

# Verify cross-model comparability after the three runs have completed.
run_model_tests <- function() {
  inputs <- data.table::fread(z_path)
  if (!identical(as.character(inputs[[1]]), names(inputs)[-1])) stop("Transaction matrix labels mismatch")
  if (length(new_labels) != 74L || RR != 2L || NN != 37L) stop("Unexpected example IO dimensions")
  for (mode in c("baseline", "climate_risk", "technology")) {
    p <- file.path(out_root, mode, "scenario_summary.csv")
    if (!file.exists(p)) stop("Missing model results for ", mode)
    a <- fread(p)
    target_n <- if (mode == "technology") 9L else 3L
    if (nrow(a) != target_n * TT) stop("Wrong row count for ", mode)
    if (any(abs(a$Storage_Demand - a$Storage_Output - a$Storage_Gap) >
            1e-7 * pmax(1, a$Storage_Demand))) stop("Invalid accounting in ", mode)
  }
  base <- fread(file.path(out_root, "baseline", "scenario_summary.csv"))
  risk <- fread(file.path(out_root, "climate_risk", "scenario_summary.csv"))
  tech <- fread(file.path(out_root, "technology", "scenario_summary.csv"))
  comparison <- merge(base[, .(Scenario, Year, Demand_Base = Storage_Demand,
                              Output_Base = Storage_Output)],
                      risk[, .(Scenario, Year, Demand_Risk = Storage_Demand,
                              Output_Risk = Storage_Output)], by = c("Scenario", "Year"))
  if (nrow(comparison) != 3L * TT) stop("Baseline/risk comparison incomplete")
  risk_effects <- comparison[, .(Output_Difference = sum(abs(Output_Base - Output_Risk))),
                             by = Scenario]
  if (any(risk_effects$Output_Difference <= 1e-7)) {
    stop("Climate risk did not alter computed output in every macro scenario")
  }
  technology_effects <- merge(
    tech[, .(Scenario, Technology, Year, Storage_Output)],
    base[, .(Scenario, Year, Base_Output = Storage_Output)],
    by = c("Scenario", "Year"))[, .(Output_Difference = sum(abs(Storage_Output - Base_Output))),
                                by = .(Scenario, Technology)]
  if (nrow(technology_effects) != 9L || any(technology_effects$Output_Difference <= 1e-7)) {
    stop("Technology inputs did not alter computed outputs")
  }
  if (any(abs(comparison$Demand_Base - comparison$Demand_Risk) > 1e-7)) {
    stop("Risk baseline must use the same reference demand")
  }
  applied <- rbindlist(lapply(scenario_table$Scenario_ID, function(id) {
    p <- file.path(out_root, "climate_risk", "scenarios", id, "S0", "Risk_Applied.csv")
    if (!file.exists(p)) stop("Missing climate risk diagnostic: ", id)
    fread(p)
  }))
  if (!nrow(applied) || !any(applied$Loss_Rate > 0)) stop("No climate loss applied")
  if (any(applied$Adjusted_Supply_Growth > applied$Base_Supply_Growth + 1e-10) ||
      !any(applied$Adjusted_Capacity_Growth < applied$Base_Capacity_Growth - 1e-10)) {
    stop("Climate risk did not change production inputs")
  }
  if (uniqueN(tech, by = c("Scenario", "Technology")) != 9L) stop("Technology cases incomplete")
  for (tid in c("S1", "S2", "S3")) {
    if (!nrow(TECH_MINERAL_INTENSITY_DT[Tech_Code == tid])) stop("Technology coefficients missing")
  }
  fwrite(comparison, file.path(out_root, "baseline_risk_comparison.csv"))
  cat("PASS: data structure, all 15 runs, loss application, scenario coverage and accounting\n")
  invisible(TRUE)
}

# Model execution and validation.
model_main <- function(mode, risk_data = NULL) {
  run_model_family(mode, risk_data)
  needed <- file.path(out_root, c("baseline", "climate_risk", "technology"), "scenario_summary.csv")
  if (all(file.exists(needed))) run_model_tests()
  cat(sprintf("Finished %s model.\n", mode))
  invisible(TRUE)
}

if ("--test-only" %in% commandArgs(trailingOnly = TRUE)) {
  run_model_tests()
} else {
  risk_data <- fread(file.path(project_root, "data", "climate_risk", "annual_mineral_losses.csv"))
  if (!all(c("Year", "Region", "Mineral_Code", "Loss_Rate") %in% names(risk_data))) {
    stop("Risk input columns are incomplete")
  }
  model_main("climate_risk", risk_data)
}
