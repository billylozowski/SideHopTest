# key updates from v3:
#
# - reshuffle of code chunks
# - generic min/max extraction added for all columns from

# ---------
# libraries
# ---------

library(tidyverse)
library(readr)
library(readxl)

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
      summarise(
         across(all_of(vars),
                list(min = safe_min,
                     max = safe_max),
                .names = paste0("{.col}_{.fn}_", suffix))
      )
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

# -----------
# read inputs
# -----------

path <- "Data/Split Data/MBB001_SHL_001.csv"
d    <- read_csv(path, col_names = FALSE)

participant_code <- substr(path, 17, 22)

height <- read_excel("Data/SHT Data_biomech.xlsx", sheet = "Study 2") %>%
   filter(Code == participant_code) %>%
   pull(Height_m)

# --------------
# clean raw data
# --------------

# build column names
row1 <- trimws(as.character(d[1, ]))
row2 <- trimws(as.character(d[2, ]))
row3 <- trimws(as.character(d[3, ]))
row4 <- trimws(as.character(d[4, ]))

colnames(d) <- ifelse(row2 == "COFP",
                      paste(row1, row2, row4, sep = "_"), # row1_COFP_row4
                      paste(row1, row4, sep = "_"))        # row1_row4

# clean the .csv
colnames(d) <- gsub("_+", "_", colnames(d))                  # collapse multiple underscores
colnames(d) <- gsub("([^_A-Za-z0-9])", "", colnames(d))      # remove stray characters
colnames(d) <- make.names(colnames(d), unique = TRUE)

d <- d %>%
   rename(FRAMES = FRAMES_0,
          TIME   = TIME_0)

d <- d[-c(1:4), ]                                            # remove the first 4 rows
names(d) <- trimws(names(d))                                 # trim stray spaces in names

d <- d %>%
   mutate(across(everything(), as.numeric)) %>%
   select(-1)                                                # remove the first column (either "NA." or "NA_NA")

# -----------------
# global parameters
# -----------------

# identify sample rate from TIME and FRAMES columns (FRAME no. -1 when TIME == 1)
sample_rate <- d %>%
   filter(TIME == 1) %>%
   pull(FRAMES) - 1

# force threshold - set as desired
fc_threshold <- 20

# all variables from COM X through to trunk tilt ML Y for generic min/max extraction
vars_min_max <- d %>%
   select(CenterOfMass_X:trunk_tilt_ml_Y) %>%
   names()

# ------------------
# detect FC / TO events
# ------------------

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
# pair events into contact and flight phases
# ------------------------------------------

# pair each foot contact with the take-off immediately following it (contact phase)
pairs_contact <- fc_mid20 %>%
   left_join(to_list, by = "FP") %>%
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

# ----------------------------------------------
# generic min/max extraction for contact / flight
# ----------------------------------------------

min_max_contact <- summarise_phase_min_max(
   phase_list = contact_list,
   phase_info = pairs_contact,
   vars       = vars_min_max,
   suffix     = "contact"
)

min_max_flight <- summarise_phase_min_max(
   phase_list = flight_list,
   phase_info = pairs_flight,
   vars       = vars_min_max,
   suffix     = "flight"
)

# ---------------------------------
# contact-specific derived variables
# ---------------------------------

# COP location (at peak GRF) within each contact phase
cop_location <- tibble(Phase_No = as.integer(names(contact_list))) %>%
   mutate(res = map(contact_list, ~{
      
      i  <- which.max(pmax(.x$FP5_Z, .x$FP6_Z))                              # index of peak vertical GRF
      fp <- if (.x$FP5_Z[i] >= .x$FP6_Z[i]) "FP5" else "FP6"                 # plate at peak GRF
      
      tibble(
         peak_grf = pmax(.x$FP5_Z[i], .x$FP6_Z[i]),                          # peak GRF magnitude
         FP       = fp,
         COPx     = if (fp == "FP5") .x$FP5_COFP_X[i] else .x$FP6_COFP_X[i], # COP x at peak
         COPy     = if (fp == "FP5") .x$FP5_COFP_Y[i] else .x$FP6_COFP_Y[i]  # COP y at peak
      )
      
   })) %>%
   unnest(res)

# create a starting contact output
final_data_contact <- pairs_contact %>%
   left_join(cop_location, by = c("Phase_No", "FP"))

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
                   sqrt(det(cov(cbind(COPx_cm, COPy_cm)))))
)

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
      })
   ) %>%
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
      COM_Z_max_cm_flight_height = (CenterOfMass_Z_max_flight * 100) / Height_m
   ) %>%
   select(Phase_No, FP, Height_m,
          COM_Path_cm_xy_flight, COM_Path_cm_xz_flight, COM_Path_cm_yz_flight,
          COM_ML_range_cm_flight, COM_ML_range_cm_flight_height,
          COM_AP_range_cm_flight, COM_AP_range_cm_flight_height,
          COM_Z_range_cm_flight,  COM_Z_range_cm_flight_height,
          COM_Z_min_cm_flight_height, COM_Z_max_cm_flight_height)

# --------------------
# assemble final outputs
# --------------------

final_data_contact <- final_data_contact %>%
   left_join(min_max_contact,   by = c("Phase_No", "FP")) %>%
   left_join(COP_COM_contact,   by = c("Phase_No", "FP")) %>%
   select(-Height_m)

final_data_flight <- pairs_flight %>%
   left_join(min_max_flight, by = c("Phase_No", "FP")) %>%
   left_join(COM_flight,     by = c("Phase_No", "FP", "Height_m")) %>%
   select(-Height_m)

# --------------
# visualisations
# --------------

# COP locations during each contact phase
ggplot(cop_location,
       aes(x = COPy,
           y = COPx,
           colour = FP)) +
   geom_point() +
   coord_equal() +
   scale_y_reverse() +
   scale_x_continuous(position = "top") +
   labs(title = "Centre of Pressure at Peak GRF",
        x     = "Mediolateral",
        y     = "Anteroposterior") +
   theme_classic() +
   theme(plot.title   = element_text(hjust = 0.5),
         aspect.ratio = 0.5)

# COP paths during each contact phase
contact_long <- bind_rows(contact_list)

contact_long_plot <- contact_long %>%
   select(Phase_No,
          FP5_COFP_X, FP5_COFP_Y,
          FP6_COFP_X, FP6_COFP_Y) %>%
   pivot_longer(cols          = -Phase_No,
                names_to      = c("Plate", ".value"),
                names_pattern = "(FP[56])_(.*)")

ggplot(contact_long_plot,
       aes(x = COFP_Y,
           y = COFP_X,
           group  = interaction(Phase_No, Plate),
           colour = Plate)) +
   geom_path(alpha = 0.7) +
   coord_equal() +
   scale_y_reverse() +
   scale_x_continuous(position = "top") +
   labs(title = "COP Path",
        x     = "Mediolateral",
        y     = "Anteroposterior") +
   theme_classic() +
   theme(aspect.ratio = 0.5,
         plot.title = element_text(hjust = 0.5))

# COM paths during each flight phase
flight_long <- bind_rows(flight_list)

# axis limits
max_x <- round(max(flight_long$CenterOfMass_Y, na.rm = TRUE) * 1.1, 1)
min_x <- round(min(flight_long$CenterOfMass_Y, na.rm = TRUE) * 0.9, 1)
max_y <- round(max(flight_long$CenterOfMass_Z, na.rm = TRUE) * 1.1, 1)
min_y <- round(min(flight_long$CenterOfMass_Z, na.rm = TRUE) * 0.9, 1)

ggplot(flight_long,
       aes(x = CenterOfMass_Y,
           y = CenterOfMass_Z,
           group  = Phase_No,
           colour = CenterOfMass_Z)) +
   geom_path(alpha = 0.7) +
   scale_x_continuous(position = "bottom",
                      limits   = c(min_x, max_x),
                      breaks   = seq(min_x, max_x, by = 0.1)) +
   scale_y_continuous(limits = c(min_y, max_y),
                      breaks = seq(min_y, max_y, by = 0.1)) +
   coord_equal() +
   labs(title = "YZ COM Path",
        x     = "Mediolateral",
        y     = "Vertical") +
   theme_classic() +
   theme(plot.title = element_text(hjust = 0.5))

# ----------------------
# future derived variables
# ----------------------
# now, I can calculate:
# - COM path (min to min? TO to TO? FC to FC?)
# - min/max joint angles & segment orientations
# - COM variables as a % of height