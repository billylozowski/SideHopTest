# key updates from v4:
#
# - hop cycles added using FC -> next FC as the phase boundary
#   allows for COM to be visualised over the whole cycle
# - master outputs saved to "Data/Processed"
# - individual trial figures saved to "Tables & Figures"
# 
# may want to label hops as lateral or medial (depending on FP and hop leg)
# i.e., FP5 for SHR would be lateral, whereas FP would be medial

# ---------
# libraries
# ---------

library(tidyverse)
library(readr)
library(readxl)
library(patchwork)

# ----------------
# helper functions
# ----------------

# triangulated path length between consecutive samples
path_length <- function(x, y) {
   dx <- diff(x)
   dy <- diff(y)
   sum(sqrt(dx^2 + dy^2), na.rm = TRUE)
}

# safe min/max for vectors that may be all NA
safe_min <- function(x) {
   if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
   if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

# extract min/max for a defined set of variables within one data frame
extract_min_max <- function(df, vars, suffix) {
   df %>%
      summarise(across(all_of(vars),
                       list(min = safe_min, max = safe_max),
                       .names = paste0("{.col}_{.fn}_", suffix)))
}

# apply min/max extraction across a list of phase-specific data frames
summarise_phase_min_max <- function(phase_list, phase_info, vars, suffix) {
   map_dfr(seq_along(phase_list), function(i) {
      extract_min_max(phase_list[[i]], vars = vars, suffix = suffix) %>%
         mutate(Phase_No = phase_info$Phase_No[i],
                FP       = phase_info$FP[i],
                .before  = 1)
   })
}

# # create square limits using a supplied span
# square_limits <- function(min_val, max_val, span) {
#    mid <- (min_val + max_val) / 2
#    c(mid - span / 2, mid + span / 2)
# }

# ------------------------
# main processing function
# ------------------------

process_file <- function(path, height_lookup_path, fc_threshold = 20) {
   
   # -----------
   # read inputs
   # -----------
   
   d <- read_csv(path, col_names = FALSE, show_col_types = FALSE)
   
   trial_name <- tools::file_path_sans_ext(basename(path))
   participant_code <- substr(trial_name, 1, 6)
   
   height <- read_excel(height_lookup_path, sheet = "Study 2") %>%
      filter(Code == participant_code) %>%
      pull(Height_m)
   
   if (length(height) != 1 || is.na(height)) {
      stop(paste("Height not found or not unique for participant:", participant_code))
   }
   
   # --------------
   # clean raw data
   # --------------
   
   # build column names
   row1 <- trimws(as.character(d[1, ]))
   row2 <- trimws(as.character(d[2, ]))
   row3 <- trimws(as.character(d[3, ]))
   row4 <- trimws(as.character(d[4, ]))
   
   colnames(d) <- ifelse(row2 == "COFP",
                         paste(row1, row2, row4, sep = "_"),    # row1_COFP_row4
                         paste(row1, row4, sep = "_"))          # row1_row4
   
   # clean the .csv
   colnames(d) <- gsub("_+", "_", colnames(d))                  # collapse multiple underscores
   colnames(d) <- gsub("([^_A-Za-z0-9])", "", colnames(d))      # remove stray characters
   colnames(d) <- make.names(colnames(d), unique = TRUE)
   
   d <- d %>%
      rename(TIME   = TIME_0)
   
   d <- d[-c(1:4), ]                                            # remove the first 4 rows
   names(d) <- trimws(names(d))                                 # trim stray spaces in names
   
   d <- d %>%
      mutate(across(everything(), as.numeric))                  # remove the first column (either "NA." or "NA_NA")
   
   # ---------------------------------
   # check required force data exist
   # ---------------------------------
   
   required_force_cols <- c("FP5_Z", "FP6_Z",
                            "FP5_COFP_X", "FP5_COFP_Y",
                            "FP6_COFP_X", "FP6_COFP_Y")
   
   missing_force_cols <- setdiff(required_force_cols, names(d))
   
   if (length(missing_force_cols) > 0) {
      stop(paste0("Missing required force columns: ",
                  paste(missing_force_cols, collapse = ", ")))
   }
   
   # -----------------
   # global parameters
   # -----------------
   
   # identify sample rate from change in TIME
   sample_rate <- round(1 / median(diff(d$TIME), na.rm = TRUE))
   
   if (length(sample_rate) != 1 || is.na(sample_rate)) {
      stop(paste("Sample rate could not be determined for:", trial_name))
   }
   
   # all variables from COM X through to trunk tilt ML Y for generic min/max extraction
   vars_min_max <- d %>%
      select(CenterOfMass_X:last_col()) %>%
      names()
   
   # --------------------
   # detect FC / TO events
   # --------------------
   
   # logical vectors above fc_threshold
   FP5_above <- d$FP5_Z >= fc_threshold
   FP6_above <- d$FP6_Z >= fc_threshold
   
   # FP5 crossings
   FP5_FC_index <- which(diff(FP5_above) == 1)  + 1 # below -> above threshold
   FP5_TO_index <- which(diff(FP5_above) == -1) + 1 # above -> below threshold
   
   # FP6 crossings
   FP6_FC_index <- which(diff(FP6_above) == 1)  + 1 # below -> above threshold
   FP6_TO_index <- which(diff(FP6_above) == -1) + 1 # above -> below threshold
   
   # ------------------------------------------------
   # consolidate events across both force plates
   # ------------------------------------------------
   
   # foot contacts
   fc_all <- bind_rows(tibble(FP = "FP5", FC_index = FP5_FC_index),
                       tibble(FP = "FP6", FC_index = FP6_FC_index)) %>%
      arrange(FC_index)
   
   # take-offs
   to_all <- bind_rows(tibble(FP = "FP5", TO_index = FP5_TO_index),
                       tibble(FP = "FP6", TO_index = FP6_TO_index)) %>%
      arrange(TO_index)
   
   if (nrow(fc_all) < 20 || nrow(to_all) < 20) {
      stop(paste0("Insufficient FC/TO events. Found ",
                  nrow(fc_all), " FC events and ",
                  nrow(to_all), " TO events."))
      }
   
   # establish the middle 20 foot contacts
   start_fc <- floor((nrow(fc_all) - 20) / 2) + 1
   start_to <- floor((nrow(to_all) - 20) / 2) + 3
   
   # select just the middle 20 foot contacts
   fc_mid20 <- fc_all %>%
      slice(start_fc:(start_fc + 19))
   
   # select the 19 take-offs within the middle 20 contact cycles
   to_mid19 <- to_all %>%
      slice(start_to:(start_to + 18))
   
   # first and last foot contacts and take-offs
   first_fc <- min(fc_mid20$FC_index)
   first_to <- min(to_mid19$TO_index)
   last_fc  <- max(fc_mid20$FC_index)
   last_to  <- max(to_mid19$TO_index)
   
   # take-off after each kept foot contact on the same plate
   to_list <- to_all %>%
      group_by(FP) %>%
      summarise(TO_List = list(sort(TO_index)), .groups = "drop")
   
   # foot contact after each kept take-off on the same plate
   fc_list <- fc_all %>%
      group_by(FP) %>%
      summarise(FC_List = list(sort(FC_index)), .groups = "drop")
   
   # ------------------------------------------
   # pair events into contact, flight, and hop phases
   # ------------------------------------------
   
   # pair each foot contact with the take-off immediately following it (contact phase)
   pairs_contact <- fc_mid20 %>%
      left_join(to_list, by  = "FP") %>%
      mutate(TO_index        = map2_dbl(FC_index, TO_List, ~ .y[.y > .x][1])) %>%
      filter(!is.na(TO_index)) %>%
      mutate(Phase_No        = row_number(),
             Phase_ID        = sprintf("GC%02d", Phase_No),
             Height_m        = height,
             Contact_Time_ms = round((TO_index - FC_index) / sample_rate * 1000, 0)) %>%
      select(Phase_No, Phase_ID, Height_m, FP, FC_index, TO_index, Contact_Time_ms)
   
   # pair each take-off with the next foot contact (flight phase)
   pairs_flight <- to_mid19 %>%
      left_join(fc_list, by = "FP") %>%
      mutate(FC_index       = map2_dbl(TO_index, FC_List, ~ .y[.y > .x][1])) %>%
      filter(!is.na(FC_index)) %>%
      mutate(Phase_No       = row_number(),
             Phase_ID       = sprintf("FL%02d", Phase_No),
             Height_m       = height,
             Flight_Time_ms = round((FC_index - TO_index) / sample_rate * 1000, 0)) %>%
      select(Phase_No, Phase_ID, Height_m, FP, TO_index, FC_index, Flight_Time_ms)
   
   # pair each foot contact with the next foot contact (hop phase)
   pairs_hop <- fc_mid20 %>%
      mutate(FC_next_index = lead(FC_index)) %>%
      filter(!is.na(FC_next_index)) %>%
      mutate(Phase_No      = row_number(),
             Phase_ID      = sprintf("HC%02d", Phase_No),
             Height_m      = height,
             Start_FP      = FP,
             Hop_Time_ms   = round((FC_next_index - FC_index) / sample_rate * 1000, 0)) %>%
      select(Phase_No, Phase_ID, Height_m, Start_FP, FC_index, FC_next_index, Hop_Time_ms)
   
   # ------------------------------------
   # create phase-specific data frame lists
   # ------------------------------------
   
   # contact phases
   contact_list <- map(seq_len(nrow(pairs_contact)), function(i) {
      d %>%
         slice(pairs_contact$FC_index[i]:pairs_contact$TO_index[i]) %>%
         mutate(Phase_No          = pairs_contact$Phase_No[i],
                Sample_in_contact = row_number())
   })
   
   names(contact_list) <- sprintf("%02d", as.integer(pairs_contact$Phase_No))
   
   # flight phases
   flight_list <- map(seq_len(nrow(pairs_flight)), function(i) {
      d %>%
         slice(pairs_flight$TO_index[i]:pairs_flight$FC_index[i]) %>%
         mutate(Phase_No         = pairs_flight$Phase_No[i],
                Sample_in_flight = row_number())
   })
   
   names(flight_list) <- sprintf("%02d", as.integer(pairs_flight$Phase_No))
   
   # hop phases
   hop_list <- map(seq_len(nrow(pairs_hop)), function(i) {
      d %>%
         slice(pairs_hop$FC_index[i]:(pairs_hop$FC_next_index[i] - 1)) %>%
         mutate(Phase_No      = pairs_hop$Phase_No[i],
                Sample_in_hop = row_number())
   })
   
   names(hop_list) <- sprintf("%02d", as.integer(pairs_hop$Phase_No))
   
   # ----------------------------------------------
   # generic min/max extraction for contact / flight / hop
   # ----------------------------------------------
   
   min_max_contact <- summarise_phase_min_max(phase_list = contact_list,
                                              phase_info = pairs_contact,
                                              vars       = vars_min_max,
                                              suffix     = "contact")
   
   min_max_flight <- summarise_phase_min_max(phase_list = flight_list,
                                             phase_info = pairs_flight,
                                             vars       = vars_min_max,
                                             suffix     = "flight")
   
   min_max_hop <- summarise_phase_min_max(phase_list = hop_list,
                                          phase_info = pairs_hop %>% rename(FP = Start_FP),
                                          vars       = vars_min_max,
                                          suffix     = "hop") %>%
      rename(Start_FP = FP)
   
   # ---------------------------------
   # contact-specific derived variables
   # ---------------------------------
   
   # COP location (at peak GRF) within each contact phase
   cop_location <- tibble(Phase_No = as.integer(names(contact_list))) %>%
      mutate(res = map(contact_list, ~{
         
         i  <- which.max(pmax(.x$FP5_Z, .x$FP6_Z))                                  # index of peak vertical GRF
         fp <- if (.x$FP5_Z[i] >= .x$FP6_Z[i]) "FP5" else "FP6"                     # plate at peak GRF
         
         tibble(peak_grf = pmax(.x$FP5_Z[i], .x$FP6_Z[i]),                          # peak GRF magnitude
                FP       = fp,
                COPx     = if (fp == "FP5") .x$FP5_COFP_X[i] else .x$FP6_COFP_X[i], # COP x at peak
                COPy     = if (fp == "FP5") .x$FP5_COFP_Y[i] else .x$FP6_COFP_Y[i]) # COP y at peak
         
      })) %>%
      unnest(res)
   
   # compute a 95% confidence ellipse
   ellipse_area_contact <- bind_rows(
      cop_location %>%
         mutate(COPx_cm = COPx * 100,
                COPy_cm = COPy * 100) %>%
         group_by(FP) %>%
         summarise(`Contact_Area_cm^2` = pi * qchisq(0.95, 2) *
                      sqrt(det(cov(cbind(COPx_cm, COPy_cm)))),
                   .groups = "drop"),
      cop_location %>%
         mutate(COPx_cm = COPx * 100,
                COPy_cm = COPy * 100) %>%
         summarise(FP = "Total",
                   `Contact_Area_cm^2` = pi * qchisq(0.95, 2) *
                      sqrt(det(cov(cbind(COPx_cm, COPy_cm))))))
   
   # though "Total" is a combination of FP5 and FP6, it's being included in column
   # "FP" for ease of reading.
   
   # contact-specific COP and COM variables
   COP_COM_contact <- pairs_contact %>%
      left_join(min_max_contact, by = c("Phase_No", "FP")) %>%
      mutate(
         # COP path (xy) in centimetres - active plate only
         COP_Path_cm_xy_contact = map2_dbl(FP, Phase_No, ~{
            df <- contact_list[[sprintf("%02d", .y)]]
            if (.x == "FP5") {
               path_length(df$FP5_COFP_X, df$FP5_COFP_Y) * 100
            } else {
               path_length(df$FP6_COFP_X, df$FP6_COFP_Y) * 100
            }
         }),
         
         # COM path (xy) in centimetres
         COM_Path_cm_xy_contact = map_dbl(Phase_No, ~{
            df <- contact_list[[sprintf("%02d", .x)]]
            path_length(df$CenterOfMass_X, df$CenterOfMass_Y) * 100
         }),
         
         # COM path (xz) in centimetres
         COM_Path_cm_xz_contact = map_dbl(Phase_No, ~{
            df <- contact_list[[sprintf("%02d", .x)]]
            path_length(df$CenterOfMass_X, df$CenterOfMass_Z) * 100
         }),
         
         # COM path (yz) in centimetres
         COM_Path_cm_yz_contact = map_dbl(Phase_No, ~{
            df <- contact_list[[sprintf("%02d", .x)]]
            path_length(df$CenterOfMass_Y, df$CenterOfMass_Z) * 100
         }),
         
         # COP ML min/max (y) in centimetres - active plate only
         COP_ML_min_cm_contact = map2_dbl(FP, Phase_No, ~{
            df <- contact_list[[sprintf("%02d", .y)]]
            if (.x == "FP5") min(df$FP5_COFP_Y, na.rm = TRUE) * 100
            else             min(df$FP6_COFP_Y, na.rm = TRUE) * 100
         }),
         COP_ML_max_cm_contact = map2_dbl(FP, Phase_No, ~{
            df <- contact_list[[sprintf("%02d", .y)]]
            if (.x == "FP5") max(df$FP5_COFP_Y, na.rm = TRUE) * 100
            else             max(df$FP6_COFP_Y, na.rm = TRUE) * 100
         }),
         
         # COP AP min/max (x) in centimetres - active plate only
         COP_AP_min_cm_contact = map2_dbl(FP, Phase_No, ~{
            df <- contact_list[[sprintf("%02d", .y)]]
            if (.x == "FP5") min(df$FP5_COFP_X, na.rm = TRUE) * 100
            else             min(df$FP6_COFP_X, na.rm = TRUE) * 100
         }),
         COP_AP_max_cm_contact = map2_dbl(FP, Phase_No, ~{
            df <- contact_list[[sprintf("%02d", .y)]]
            if (.x == "FP5") max(df$FP5_COFP_X, na.rm = TRUE) * 100
            else             max(df$FP6_COFP_X, na.rm = TRUE) * 100
         })) %>%
      mutate(
         # centre of pressure
         COP_ML_range_cm_contact        = COP_ML_max_cm_contact - COP_ML_min_cm_contact,
         COP_ML_range_cm_contact_height = COP_ML_range_cm_contact / Height_m,
         COP_AP_range_cm_contact        = COP_AP_max_cm_contact - COP_AP_min_cm_contact,
         COP_AP_range_cm_contact_height = COP_AP_range_cm_contact / Height_m,
         
         # centre of mass (generic min/max columns are in metres, so convert to cm)
         COM_ML_range_cm_contact        = (CenterOfMass_Y_max_contact - CenterOfMass_Y_min_contact) * 100,
         COM_ML_range_cm_contact_height = COM_ML_range_cm_contact / Height_m,
         COM_AP_range_cm_contact        = (CenterOfMass_X_max_contact - CenterOfMass_X_min_contact) * 100,
         COM_AP_range_cm_contact_height = COM_AP_range_cm_contact / Height_m,
         COM_Z_range_cm_contact         = (CenterOfMass_Z_max_contact - CenterOfMass_Z_min_contact) * 100,
         COM_Z_range_cm_contact_height  = COM_Z_range_cm_contact / Height_m,
         
         # COM to COP ratio (raw)
         COM_COP_xy_Ratio_contact = COM_Path_cm_xy_contact / COP_Path_cm_xy_contact) %>%
      
      select(Phase_No, FP,
             COP_Path_cm_xy_contact,
             COM_Path_cm_xy_contact, COM_Path_cm_xz_contact, COM_Path_cm_yz_contact,
             COM_COP_xy_Ratio_contact,
             COP_ML_min_cm_contact, COP_ML_max_cm_contact, COP_ML_range_cm_contact,
             COP_ML_range_cm_contact_height,
             COP_AP_min_cm_contact, COP_AP_max_cm_contact, COP_AP_range_cm_contact,
             COP_AP_range_cm_contact_height,
             COM_ML_range_cm_contact, COM_ML_range_cm_contact_height,
             COM_AP_range_cm_contact, COM_AP_range_cm_contact_height,
             COM_Z_range_cm_contact, COM_Z_range_cm_contact_height)
   
   # --------------------------------
   # flight-specific derived variables
   # --------------------------------
   
   COM_flight <- min_max_flight %>%
      left_join(pairs_flight %>% select(Phase_No, FP, Height_m), by = c("Phase_No", "FP")) %>%
      mutate(
         # COM path (xy) in centimetres
         COM_Path_cm_xy_flight = map_dbl(Phase_No, ~{
            df <- flight_list[[sprintf("%02d", .x)]]
            path_length(df$CenterOfMass_X, df$CenterOfMass_Y) * 100
         }),
         
         # COM path (xz) in centimetres
         COM_Path_cm_xz_flight = map_dbl(Phase_No, ~{
            df <- flight_list[[sprintf("%02d", .x)]]
            path_length(df$CenterOfMass_X, df$CenterOfMass_Z) * 100
         }),
         
         # COM path (yz) in centimetres
         COM_Path_cm_yz_flight = map_dbl(Phase_No, ~{
            df <- flight_list[[sprintf("%02d", .x)]]
            path_length(df$CenterOfMass_Y, df$CenterOfMass_Z) * 100
         }),
         
         # absolute range (converted from metres to centimetres)
         COM_ML_range_cm_flight = (CenterOfMass_Y_max_flight - CenterOfMass_Y_min_flight) * 100,
         COM_AP_range_cm_flight = (CenterOfMass_X_max_flight - CenterOfMass_X_min_flight) * 100,
         COM_Z_range_cm_flight  = (CenterOfMass_Z_max_flight - CenterOfMass_Z_min_flight) * 100,
         
         # relative range (% height)
         COM_ML_range_cm_flight_height = COM_ML_range_cm_flight / Height_m,
         COM_AP_range_cm_flight_height = COM_AP_range_cm_flight / Height_m,
         COM_Z_range_cm_flight_height  = COM_Z_range_cm_flight / Height_m,
         
         # additional COM Z directional variables (% height)
         COM_Z_min_cm_flight_height = (CenterOfMass_Z_min_flight * 100) / Height_m,
         COM_Z_max_cm_flight_height = (CenterOfMass_Z_max_flight * 100) / Height_m) %>%
      
      select(Phase_No, FP, Height_m,
             COM_Path_cm_xy_flight, COM_Path_cm_xz_flight, COM_Path_cm_yz_flight,
             COM_ML_range_cm_flight, COM_ML_range_cm_flight_height,
             COM_AP_range_cm_flight, COM_AP_range_cm_flight_height,
             COM_Z_range_cm_flight,  COM_Z_range_cm_flight_height,
             COM_Z_min_cm_flight_height, COM_Z_max_cm_flight_height)
   
   # -----------------------------
   # hop phase-specific variables
   # -----------------------------
   
   COM_hop <- min_max_hop %>%
      left_join(pairs_hop %>% select(Phase_No, Height_m), by = "Phase_No") %>%
      mutate(
         # COM path (xy) in centimetres
         COM_Path_cm_xy_hop = map_dbl(Phase_No, ~{
            df <- hop_list[[sprintf("%02d", .x)]]
            path_length(df$CenterOfMass_X, df$CenterOfMass_Y) * 100
         }),
         
         # COM path (xz) in centimetres
         COM_Path_cm_xz_hop = map_dbl(Phase_No, ~{
            df <- hop_list[[sprintf("%02d", .x)]]
            path_length(df$CenterOfMass_X, df$CenterOfMass_Z) * 100
         }),
         
         # COM path (yz) in centimetres
         COM_Path_cm_yz_hop = map_dbl(Phase_No, ~{
            df <- hop_list[[sprintf("%02d", .x)]]
            path_length(df$CenterOfMass_Y, df$CenterOfMass_Z) * 100
         }),
         
         # absolute range (converted from metres to centimetres)
         COM_ML_range_cm_hop = (CenterOfMass_Y_max_hop - CenterOfMass_Y_min_hop) * 100,
         COM_AP_range_cm_hop = (CenterOfMass_X_max_hop - CenterOfMass_X_min_hop) * 100,
         COM_Z_range_cm_hop  = (CenterOfMass_Z_max_hop - CenterOfMass_Z_min_hop) * 100,
         
         # relative range (% height)
         COM_ML_range_cm_hop_height = COM_ML_range_cm_hop / Height_m,
         COM_AP_range_cm_hop_height = COM_AP_range_cm_hop / Height_m,
         COM_Z_range_cm_hop_height  = COM_Z_range_cm_hop / Height_m,
         
         # additional COM Z directional variables (% height)
         COM_Z_min_cm_hop_height = (CenterOfMass_Z_min_hop * 100) / Height_m,
         COM_Z_max_cm_hop_height = (CenterOfMass_Z_max_hop * 100) / Height_m) %>%
      
      select(Phase_No, Start_FP, Height_m,
             COM_Path_cm_xy_hop, COM_Path_cm_xz_hop, COM_Path_cm_yz_hop,
             COM_ML_range_cm_hop, COM_ML_range_cm_hop_height,
             COM_AP_range_cm_hop, COM_AP_range_cm_hop_height,
             COM_Z_range_cm_hop,  COM_Z_range_cm_hop_height,
             COM_Z_min_cm_hop_height, COM_Z_max_cm_hop_height)
   
   # --------------------
   # assemble final outputs
   # --------------------
   
   final_data_contact <- pairs_contact %>%
      left_join(cop_location, by = c("Phase_No", "FP")) %>%
      left_join(min_max_contact, by = c("Phase_No", "FP")) %>%
      left_join(COP_COM_contact, by = c("Phase_No", "FP")) %>%
      mutate(Trial_Name = trial_name, .before = 1)
   
   final_data_flight <- pairs_flight %>%
      left_join(min_max_flight, by = c("Phase_No", "FP")) %>%
      left_join(COM_flight, by = c("Phase_No", "FP", "Height_m")) %>%
      mutate(Trial_Name = trial_name, .before = 1)
   
   final_data_hop <- pairs_hop %>%
      left_join(min_max_hop, by = c("Phase_No", "Start_FP")) %>%
      left_join(COM_hop, by = c("Phase_No", "Start_FP", "Height_m")) %>%
      mutate(Trial_Name = trial_name, .before = 1)
   
   # --------------
   # visualisations
   # --------------
   
   # COP paths during each contact phase
   contact_long <- bind_rows(contact_list)
   
   contact_long_plot <- contact_long %>%
      select(Phase_No,
             FP5_COFP_X, FP5_COFP_Y,
             FP6_COFP_X, FP6_COFP_Y) %>%
      pivot_longer(cols          = -Phase_No,
                   names_to      = c("Plate", ".value"),
                   names_pattern = "(FP[56])_(.*)")
   

   # 2D COP locations (xy) during each contact phase
   plot_cop_location <- ggplot() +
      # COP locations at peak GRF
      geom_point(data = cop_location,
                 aes(x = COPy, y = COPx, colour = FP),
                 size = 2) +
      # ellipses per plate at FC
      stat_ellipse(data   = cop_location,
                   aes(x  = COPy, y = COPx, group = FP),
                   colour = "black", level = 0.95, linewidth = 0.5) +
      coord_equal() +
      scale_y_reverse() +
      scale_x_continuous(position = "top") +
      labs(title = "Centre of Pressure at Peak GRF",
           x     = "Mediolateral",
           y     = "Anteroposterior") +
      theme_classic() +
      theme(plot.title   = element_text(hjust = 0.5),
            aspect.ratio = 0.5)
   
   # 2D COP paths (xy) during each contact phase
   plot_cop_path <- ggplot() +
      geom_path(data  = contact_long_plot,
                aes(x = COFP_Y, y = COFP_X, group  = interaction(Phase_No, Plate), colour = Plate),
                alpha = 0.7) +
      # COP points at FC (first row of each contact phase)
      geom_point(data = contact_long %>%
                    group_by(Phase_No) %>%
                    slice(1) %>%
                    ungroup() %>%
                    left_join(pairs_contact %>% select(Phase_No, FP), by = "Phase_No") %>%
                    mutate(COPx  = if_else(FP == "FP5", FP5_COFP_X, FP6_COFP_X),
                           COPy  = if_else(FP == "FP5", FP5_COFP_Y, FP6_COFP_Y)),
                 aes(x = COPy, y = COPx, colour = FP),
                 size = 2, shape = 4) +
      coord_equal() +
      scale_y_reverse() +
      scale_x_continuous(position = "top") +
      labs(title = "COP Paths During Each Contact Phase",
           x     = "Mediolateral",
           y     = "Anteroposterior") +
      theme_classic() +
      theme(aspect.ratio = 0.5,
            plot.title   = element_text(hjust = 0.5))
   
   # COM paths during each hop phase
   hop_long <- bind_rows(hop_list)
   
   # fixed square limits
   
   yz_x_limits <- c(0.3, 1.0)
   yz_y_limits <- c(0.7, 1.4)
   
   xz_x_limits <- c(0.0, 0.6)
   xz_y_limits <- c(0.7, 1.4)
   
   xy_x_limits <- c(0.4, 0.9)
   xy_y_limits <- c(0.0, 0.5)
     
   # frontal COM paths (yz) during each hop phase
   plot_com_yz_hop <- ggplot(hop_long,
                             aes(x = CenterOfMass_Y,
                                 y = CenterOfMass_Z,
                                 group  = Phase_No,
                                 colour = CenterOfMass_Z)) +
      geom_path(alpha = 0.7) +
      scale_x_continuous(position = "bottom",
                         limits   = yz_x_limits,
                         breaks   = seq(yz_x_limits[1], yz_x_limits[2], by = 0.1),
                         labels   = \(x) sprintf("%.1f", x)) +
      scale_y_continuous(limits   = yz_y_limits,
                         breaks   = seq(yz_y_limits[1], yz_y_limits[2], by = 0.1),
                         labels   = \(x) sprintf("%.1f", x)) +
      coord_equal() +
      labs(title = "YZ COM Path Across Hop Cycle",
           x     = "Mediolateral",
           y     = "Vertical") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 8)))
   
   # sagittal COM paths (xz) during each hop phase
   plot_com_xz_hop <- ggplot(hop_long,
                             aes(x = CenterOfMass_X,
                                 y = CenterOfMass_Z,
                                 group  = Phase_No,
                                 colour = CenterOfMass_Z)) +
      geom_path(alpha = 0.7) +
      scale_x_continuous(position = "bottom",
                         limits   = xz_x_limits,
                         breaks   = seq(xz_x_limits[1], xz_x_limits[2], by = 0.1),
                         labels = \(x) sprintf("%.1f", x)) +
      scale_y_continuous(limits = xz_y_limits,
                         breaks = seq(xz_y_limits[1], xz_y_limits[2], by = 0.1),
                         labels = \(x) sprintf("%.1f", x)) +
      coord_equal() +
      labs(title = "XZ COM Path Across Hop Cycle",
           x     = "Anteroposterior",
           y     = "Vertical") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 8)))
   
   # transverse COM paths (xy) during each hop phase
   plot_com_xy_hop <- ggplot(hop_long,
                             aes(x = CenterOfMass_Y,
                                 y = CenterOfMass_X,
                                 group  = Phase_No,
                                 colour = CenterOfMass_X)) +
      geom_path(alpha = 0.7) +
      scale_x_continuous(position = "top",
                         limits   = xy_x_limits,
                         breaks   = seq(xy_x_limits[1], xy_x_limits[2], by = 0.1),
                         labels   = \(x) sprintf("%.1f", x)) +
      scale_y_continuous(limits   = xy_y_limits,
                         breaks   = seq(xy_y_limits[1], xy_y_limits[2], by = 0.1),
                         labels   = \(x) sprintf("%.1f", x)) +
      coord_equal() +
      scale_y_reverse() +
      labs(title = "XY COM Path Across Hop Cycle",
           x     = "Mediolateral",
           y     = "Anteroposterior") +
      theme_classic() +
      theme(plot.title   = element_text(hjust = 0.5, margin = margin(b = 8)),
            aspect.ratio = 0.5)
   
   # composite COP figure
   plot_composite <- plot_cop_location / plot_cop_path +
      plot_layout(heights = c(1, 1)) +
      plot_annotation(tag_levels = "A")
   
   # composite COM figure
   plot_composite_com <- ((plot_com_yz_hop + theme(legend.position = "none")) |
                          (plot_com_xz_hop + theme(legend.position = "none"))) / 
                          (plot_com_xy_hop + theme(legend.position = "none")) +
      plot_layout(widths = c(1, 1.1), heights = c(1, 1)) +
      plot_annotation(tag_levels = "A") &
      theme(plot.margin = margin(5, 5, 5, 5))

   # ------
   # return
   # ------
   
   list(trial_name           = trial_name,
        participant_code     = participant_code,
        final_data_contact   = final_data_contact,
        final_data_flight    = final_data_flight,
        final_data_hop       = final_data_hop,
        ellipse_area_contact = ellipse_area_contact,
        plot_cop_location    = plot_cop_location,
        plot_cop_path        = plot_cop_path,
        plot_com_yz_hop      = plot_com_yz_hop,
        plot_com_xz_hop      = plot_com_xz_hop,
        plot_com_xy_hop      = plot_com_xy_hop,
        plot_composite       = plot_composite,
        plot_composite_com   = plot_composite_com)
   }

# ----------------
# output folders
# ----------------

dir.create("Data/Processed", showWarnings = FALSE)
dir.create("Tables & Figures", showWarnings = FALSE)

# ---------------
# process all files
# ---------------

csv.files <- list.files("Data/Split Data", pattern = "\\.csv$", full.names = TRUE)

results <- list()
skipped_trials <- tibble(Trial_Name = character(),
                         File_Path  = character(),
                         Reason     = character())

for (csv.file in csv.files) {
   
   trial_name <- tools::file_path_sans_ext(basename(csv.file))
   
   res <- tryCatch(
      process_file(path               = csv.file,
                   height_lookup_path = "Data/SHT Data_biomech.xlsx",
                   fc_threshold       = 20),
      error = function(e) {
         skipped_trials <<- bind_rows(skipped_trials,
                                      tibble(Trial_Name = trial_name,
                                             File_Path  = csv.file,
                                             Reason     = conditionMessage(e)))
         NULL
         }
      )
   
   if (is.null(res)) next
   
   results[[res$trial_name]] <- res
   
   #save the COP plot
   ggsave(filename = file.path("Tables & Figures",
                           paste0(res$trial_name, " - COP.png")),
          plot = res$plot_composite, width = 9, height = 9, dpi = 600)
   
   # save the COM plot
   ggsave(filename = file.path("Tables & Figures",
                               paste0(res$trial_name, " - COM.png")),
          plot = res$plot_composite_com, width = 9, height = 9, dpi = 600)
}

# --------------------------
# combined master data files
# --------------------------

contact_all <- map_dfr(results, "final_data_contact")
flight_all  <- map_dfr(results, "final_data_flight")
hop_all     <- map_dfr(results, "final_data_hop")

write_csv(contact_all, "Data/Processed/Contact Phase.csv")
write_csv(flight_all,  "Data/Processed/Flight Phase.csv")
write_csv(hop_all,     "Data/Processed/Hop Cycle.csv")

# ----------------
# skipped trial log
# ----------------

write_csv(skipped_trials, "Data/Processed/Skipped Trials.csv")