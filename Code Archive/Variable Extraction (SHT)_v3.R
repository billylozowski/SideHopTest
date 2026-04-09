# key updates from v2:
# 
# - flight phases calculated (19 amongst middle 20 foot contacts)
# - phase (contact) appended to end of variable names

library(tidyverse)
library(readr)
library(readxl)

# ------------------------
# read and clean the data
# ------------------------

path <- "Data/Split Data/MBB001_SHL_001.csv"
d      <- read_csv(path, col_names = FALSE)

participant_code <- substr(path, 17,22) 

height <- read_excel("Data/SHT Data_biomech.xlsx", sheet = "Study 2") %>%
   filter(Code == participant_code) %>%
   pull(Height_m)
      

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


# ----------------------
# set parameters
# ----------------------
# identify sample rate from TIME and FRAMES columns (FRAME no. -1 when TIME == 1)
sample_rate <- d %>%
   filter(TIME == 1) %>% # just the row at 1s
   pull(FRAMES) - 1      # get the value in FRAMES and subtract 1

# force threshold - set as desired
fc_threshold <- 20

# --------------------------------------------
# foot contacts and take-offs for FP5 and FP6
# --------------------------------------------
# logical vector above fc_threshold
FP5_above <- d$FP5_Z >= fc_threshold
FP6_above <- d$FP6_Z >= fc_threshold

# FP5 crossings
FP5_FC_index <- which(diff(FP5_above) == 1)  + 1 # all instances where FP5_Z moves below → above fc_threshold
FP5_TO_index <- which(diff(FP5_above) == -1) + 1 # all instances where FP5_Z moves above → below fc_threshold

# FP6 crossings
FP6_FC_index <- which(diff(FP6_above) == 1)  + 1 # all instances where FP5_Z moves below → above fc_threshold
FP6_TO_index <- which(diff(FP6_above) == -1) + 1 # all instances where FP5_Z moves above → below fc_threshold

# diff() checks changes between consecutive logical samples (0, 0, 0, *1, 1, ... , 1, 1, *0, 0, 0, ...)
# 0  →  1 transition returns +1 (foot contact)
# 1  →  0 transition returns -1 (take-off)
# which() returns indices where transitions occur
# "+1" shifts to the first new-phase sample

# -------------------------------------------------------
# combine foot contacts and takeoffs across force plates
# -------------------------------------------------------
# foot contacts
fc_all <- bind_rows(tibble(FP = "FP5", FC_index = FP5_FC_index),
                    tibble(FP = "FP6", FC_index = FP6_FC_index)) %>%
   arrange(FC_index)

# takeoffs
to_all <- bind_rows(tibble(FP = "FP5", TO_index = FP5_TO_index),
                    tibble(FP = "FP6", TO_index = FP6_TO_index)) %>%
   arrange(TO_index)

# establish the middle 20 foot contacts
start_fc <- floor((nrow(fc_all) - 20) / 2) + 1
start_to <- floor((nrow(to_all) - 20) / 2) + 3 

# select just the middle 20 foot contacts
fc_mid20 <- fc_all %>%
   slice(start_fc:(start_fc + 19))

# select the 19 takeoffs within the middle 20 contact cycles
to_mid19 <- to_all %>%
   slice(start_to:(start_to + 18))

# first and last foot contacts an takeoffs
first_fc <- min(fc_mid20$FC_index)
first_to <- min(to_mid19$TO_index)
last_fc  <- max(fc_mid20$FC_index)
last_to  <- max(to_mid19$TO_index)

# takeoff after each kept foot contact on the same plate
to_list <- to_all %>%
   group_by(FP) %>%
   summarise(TO_List = list(sort(TO_index)), .groups = "drop")

# foot contact after each kept foot contact on the same plate
fc_list <- fc_all %>%
   group_by(FP) %>%
   summarise(FC_List = list(sort(FC_index)), .groups = "drop")

# ----------------------------------------------------
# group each foot contact with its respective takeoff
# ----------------------------------------------------
# pair each foot contact with the takeoff immediately following it (contact cycle)
pairs_contact <- fc_mid20 %>%
   left_join(to_list, by = "FP") %>%
   mutate(TO_index        = map2_dbl(FC_index, TO_List, ~ .y[.y > .x][1])) %>%
   filter(!is.na(TO_index)) %>%
   mutate(Phase_No        = row_number(),                  # numeric row ID
          Phase_ID        = sprintf("GC%02d", Phase_No),   # zero-padded row ID (i.e. 01, 02, ...)
          Height_m        = height,                        # create a height column
          Contact_Time_ms = round((TO_index - FC_index) / sample_rate * 1000, 0)) %>%
   select(Phase_No, Phase_ID, Height_m, FP, FC_index, TO_index, Contact_Time_ms)

# do the same for takeoff and foot contact (flight cycle) - only 19 though
pairs_flight <- to_mid19 %>%
   left_join(fc_list, by = "FP") %>%
   mutate(FC_index       = map2_dbl(TO_index, FC_List, ~ .y[.y > .x][1])) %>%
   filter(!is.na(FC_index)) %>%
   mutate(Phase_No       = row_number(),                  # numeric row ID
          Phase_ID       = sprintf("FL%02d", Phase_No),   # zero-padded row ID (i.e. 01, 02, ...)
          Height_m       = height,                        # create a height column
          Flight_Time_ms = round((FC_index - TO_index) / sample_rate * 1000, 0)) %>%
   select(Phase_No, Phase_ID, Height_m, FP, TO_index, FC_index, Flight_Time_ms)

# -------------------------------------------------------------
# create a subset data frame for each contact and flight cycle
# -------------------------------------------------------------
# contacts phases
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
      mutate(Phase_No          = pairs_flight$Phase_No[i],
             Sample_in_flights = row_number())
})

names(flight_list) <- sprintf("%02d", as.integer(pairs_flight$Phase_No))

# -----------------------------------------------------
# COP location (at peak GRF) within each contact cycle
# -----------------------------------------------------

cop_location <- tibble(Phase_No = as.integer(names(contact_list))) %>%       # initialise output with contact labels
   mutate(res = map(contact_list, ~{                                         # iterate over each contact cycle
      
      i  <- which.max(pmax(.x$FP5_Z, .x$FP6_Z))                              # index of peak vertical GRF (either plate)
      fp <- if (.x$FP5_Z[i] >= .x$FP6_Z[i]) "FP5" else "FP6"                 # determine which plate produced peak
      
      tibble(                                                                # return one-row tibble per contact
         peak_grf = pmax(.x$FP5_Z[i], .x$FP6_Z[i]),                          # peak GRF magnitude
         FP       = fp,                                                      # plate at peak GRF
         COPx     = if (fp == "FP5") .x$FP5_COFP_X[i] else .x$FP6_COFP_X[i], # COP x at peak
         COPy     = if (fp == "FP5") .x$FP5_COFP_Y[i] else .x$FP6_COFP_Y[i]) # COP y at peak
      
   })) %>%
   unnest(res)                                                               # expand list-column into standard columns

# ------------------------------------------------
# create a data frame for all extracted variables
# ------------------------------------------------

final_data_contact <- pairs_contact %>%
   left_join(cop_location, by = c("Phase_No", "FP"))

# ---------------------------------
# compute a 95% confidence ellipse 
# ---------------------------------

ellipse_area_contact <- bind_rows(cop_location %>%
      # convert COP values to centimetres
      mutate(COPx_cm = COPx * 100,
             COPy_cm = COPy * 100) %>%
      group_by(FP) %>%
      # contact area per plate
      summarise(`Contact_Area_cm^2` = pi * qchisq(0.95, 2) *
                   sqrt(det(cov(cbind(COPx_cm, COPy_cm)))),
                .groups = "drop"),
   cop_location %>%
      mutate(COPx_cm = COPx * 100,
             COPy_cm = COPy * 100) %>%
      # total contact area
      summarise(FP = "Total",
                `Contact_Area_cm^2` = pi * qchisq(0.95, 2) *
                   sqrt(det(cov(cbind(COPx_cm, COPy_cm))))))

# though "Total" is a combination of FP5 and FP6, it's being included in column 
# "FP" for ease of reading.

# --------------------------------
# visualise the COP location data
# --------------------------------
# COP locations during each contact cycle
ggplot(cop_location, 
       aes(x = COPy, 
           y = COPx, 
           color = FP)) +
   geom_point() +
   coord_equal() +
   scale_y_reverse() +                    # invert y axis
   scale_x_continuous(position = "top") + # move x axis to top of plot
   labs(title = "Centre of Pressure at Peak GRF",
        x     = "Mediolateral",
        y     = "Anteroposterior") +
   theme_classic() +
   theme(plot.title = element_text(hjust = 0.5),
         aspect.ratio = 0.5)              # 2:1 width vs. height ratio

# COP paths during each contact cycle
contact_long <- bind_rows(contact_list, .id = "Phase_No")

contact_long_plot <- contact_long %>%
   select(Phase_No,
          FP5_COFP_X, FP5_COFP_Y,
          FP6_COFP_X, FP6_COFP_Y) %>%
   pivot_longer(cols          = -Phase_No,
                names_to      = c("Plate", ".value"),
                names_pattern = "(FP[56])_(.*)")

# build the plot
ggplot(contact_long_plot,
       aes(x = COFP_Y,                    # ml direction
           y = COFP_X,                    # ap direction
           group  = interaction(Phase_No, Plate),
           colour = Plate)) +
   geom_path(alpha = 0.7) +               # ml direction
   # I want to add a marker for the first row of data (x and y) to see where path starts
   # geom_point(x = COFP_Y[1],                    
   #            y = COFP_X[1]) +
   coord_equal() +
   scale_y_reverse() +                    # invert y axis
   scale_x_continuous(position = "top") + # move x axis to top of plot
   labs(title = "COP Path",
        x     = "Mediolateral",
        y     = "Anteroposterior") +
   theme_classic() +
   theme(aspect.ratio = 0.5)              # 2:1 width vs. height ratio

# --------------------------------------------------
# calculate COP & COM paths during each phase cycle
# --------------------------------------------------
# create a function for computing the triangulated distance between rows
path_length <- function(x, y) {
   dx <- diff(x)
   dy <- diff(y)
   sum(sqrt(dx^2 + dy^2), na.rm = TRUE)
}

# ---------------------------------
# contact phase-specific variables
# ---------------------------------

COP_COM_contact <- pairs_contact %>%
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
      # COP ML min (y) in centimetres
      COP_ML_min_cm_contact = map2_dbl(FP, Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .y)]]
         if (.x == "FP5") min(df$FP5_COFP_Y, na.rm = TRUE) * 100
         else             min(df$FP6_COFP_Y, na.rm = TRUE) * 100
      }),
      # COP ML max (y) in centimetres
      COP_ML_max_cm_contact = map2_dbl(FP, Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .y)]]
         if (.x == "FP5") max(df$FP5_COFP_Y, na.rm = TRUE) * 100
         else             max(df$FP6_COFP_Y, na.rm = TRUE) * 100
      }),
      # COP AP min (x) in centimetres
      COP_AP_min_cm_contact = map2_dbl(FP, Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .y)]]
         if (.x == "FP5") min(df$FP5_COFP_X, na.rm = TRUE) * 100
         else             min(df$FP6_COFP_X, na.rm = TRUE) * 100
      }),
      # COP AP max (x) in centimetres
      COP_AP_max_cm_contact = map2_dbl(FP, Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .y)]]
         if (.x == "FP5") max(df$FP5_COFP_X, na.rm = TRUE) * 100
         else             max(df$FP6_COFP_X, na.rm = TRUE) * 100
      }),
      # COM ML min (y) in centimetres
      COM_ML_min_cm_contact = map_dbl(Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .x)]]
         min(df$CenterOfMass_Y, na.rm = TRUE) * 100
      }),
      # COM ML max (y) in centimetres
      COM_ML_max_cm_contact = map_dbl(Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .x)]]
         max(df$CenterOfMass_Y, na.rm = TRUE) * 100
      }),
      # COM AP min (x) in centimetres
      COM_AP_min_cm_contact = map_dbl(Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .x)]]
         min(df$CenterOfMass_X, na.rm = TRUE) * 100
      }),
      # COM AP max (x) in centimetres
      COM_AP_max_cm_contact = map_dbl(Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .x)]]
         max(df$CenterOfMass_X, na.rm = TRUE) * 100
      }),
      # COM vertical min (z) in centimetres
      COM_Z_min_cm_contact = map_dbl(Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .x)]]
         min(df$CenterOfMass_Z, na.rm = TRUE) * 100
      }),
      # COM vertical max (z) in centimetres
      COM_Z_max_cm_contact = map_dbl(Phase_No, ~{
         df <- contact_list[[sprintf("%02d", .x)]]
         max(df$CenterOfMass_Z, na.rm = TRUE) * 100
      })) %>%
   mutate(
      # centre of pressure
      COP_ML_range_cm_contact    = COP_ML_max_cm_contact - COP_ML_min_cm_contact, # COP ML range
      COP_ML_range_cm_contact_height = COP_ML_range_cm_contact / Height_m,         # as a % of height   
      COP_AP_range_cm_contact    = COP_AP_max_cm_contact - COP_AP_min_cm_contact, # COP AP range
      COP_AP_range_cm_contact_height = COP_AP_range_cm_contact / Height_m,         # as a % of height  
      # centre of mass
      COM_ML_range_cm_contact    = COM_ML_max_cm_contact - COM_ML_min_cm_contact, # COM ML range
      COM_ML_range_cm_contact_height = COM_ML_range_cm_contact / Height_m,         # as a % of height  
      COM_AP_range_cm_contact    = COM_AP_max_cm_contact - COM_AP_min_cm_contact, # COM AP range
      COM_AP_range_cm_contact_height = COM_AP_range_cm_contact / Height_m,         # as a % of height  
      COM_Z_range_cm_contact     = COM_Z_max_cm_contact  - COM_Z_min_cm_contact,  # COM  Z range
      COM_Z_range_cm_contact_height = COM_Z_range_cm_contact / Height_m,           # as a % of height  
      # COM to COP ratio (raw)
      COM_COP_xy_Ratio_contact = COM_Path_cm_xy_contact / COP_Path_cm_xy_contact,
      ) %>% 
   select(Phase_No, FP, 
          COP_Path_cm_xy_contact, COM_Path_cm_xy_contact, COM_COP_xy_Ratio_contact,
          COM_Path_cm_xz_contact, COM_Path_cm_yz_contact,
          COP_ML_min_cm_contact,  COP_ML_max_cm_contact,  COP_ML_range_cm_contact, # COP ML
          COP_AP_min_cm_contact,  COP_AP_max_cm_contact,  COP_AP_range_cm_contact, # COP AP
          COM_ML_min_cm_contact,  COM_ML_max_cm_contact,  COM_ML_range_cm_contact, # COM ML
          COM_AP_min_cm_contact,  COM_AP_max_cm_contact,  COM_AP_range_cm_contact, # COM AP
          COM_Z_min_cm_contact,   COM_Z_max_cm_contact,   COM_Z_range_cm_contact)  # COM Z

# join COP/COM data to final_data_contact
final_data_contact <- final_data_contact %>%
   left_join(COP_COM_contact, by = c("Phase_No", "FP"))

# --------------------------------
# flight phase-specific variables
# --------------------------------

COM_flight <- pairs_flight %>%
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
      # COM ML (y) min in centimetres
      COM_ML_min_cm_flight = map_dbl(Phase_No, ~{
         df <- flight_list[[sprintf("%02d", .x)]]
         min(df$CenterOfMass_Y, na.rm = TRUE) * 100
      }),
      # COM ML (y) max in centimetres
      COM_ML_max_cm_flight = map_dbl(Phase_No, ~{
         df <- flight_list[[sprintf("%02d", .x)]]
         max(df$CenterOfMass_Y, na.rm = TRUE) * 100
      }),
      # COM AP (x) min in centimetres
      COM_AP_min_cm_flight = map_dbl(Phase_No, ~{
         df <- flight_list[[sprintf("%02d", .x)]]
         min(df$CenterOfMass_X, na.rm = TRUE) * 100
      }),
      # COM AP (x) max in centimetres
      COM_AP_max_cm_flight = map_dbl(Phase_No, ~{
         df <- flight_list[[sprintf("%02d", .x)]]
         max(df$CenterOfMass_X, na.rm = TRUE) * 100
      }),
      # COM vertical (z) min in centimetres
      COM_Z_min_cm_flight = map_dbl(Phase_No, ~{
         df <- flight_list[[sprintf("%02d", .x)]]
         min(df$CenterOfMass_Z, na.rm = TRUE) * 100
      }),
      # COM vertical (z) max in centimetres
      COM_Z_max_cm_flight = map_dbl(Phase_No, ~{
         df <- flight_list[[sprintf("%02d", .x)]]
         max(df$CenterOfMass_Z, na.rm = TRUE) * 100
      })
      ) %>%
   mutate(
      # absolute range
      COM_ML_range_cm_flight = COM_ML_max_cm_flight - COM_ML_min_cm_flight,
      COM_AP_range_cm_flight = COM_AP_max_cm_flight - COM_AP_min_cm_flight,
      COM_Z_range_cm_flight  = COM_Z_max_cm_flight  - COM_Z_min_cm_flight,
      # relative range (% height)
      COM_ML_range_cm_flight_height = COM_ML_range_cm_flight / Height_m,
      COM_AP_range_cm_flight_height = COM_AP_range_cm_flight / Height_m,
      COM_Z_range_cm_flight_height  = COM_Z_range_cm_flight  / Height_m,
      # additional COM Z directional variables (% height)
      COM_Z_min_cm_flight_height = COM_Z_min_cm_flight / Height_m,
      COM_Z_max_cm_flight_height = COM_Z_max_cm_flight / Height_m) %>%
   
   select(Phase_No, FP, Height_m,
          COM_Path_cm_xy_flight, COM_Path_cm_xz_flight, COM_Path_cm_yz_flight,  # COM paths
          COM_ML_min_cm_flight,  COM_ML_max_cm_flight,  COM_ML_range_cm_flight, # COM ML
          COM_ML_range_cm_flight_height,
          COM_AP_min_cm_flight,  COM_AP_max_cm_flight,  COM_AP_range_cm_flight, # COM AP
          COM_AP_range_cm_flight_height,
          COM_Z_min_cm_flight,   COM_Z_max_cm_flight,   COM_Z_range_cm_flight,  # COM Z
          COM_Z_range_cm_flight_height,
          COM_Z_min_cm_flight_height, COM_Z_max_cm_flight_height)
          
# join COP/COM data to final_data_flight
final.data.flight <- pairs_flight %>%
   left_join(COM_flight, by = c("Phase_No", "FP", "Height_m"))

# --------------------
# visualise COM paths
# --------------------
# COM paths during each flight cycle
flight_long <- bind_rows(flight_list, .id = "Phase_No")

# axis limits
max_x <- round(max(flight_long$CenterOfMass_Y) * 1.1, 1)
min_x <- round(min(flight_long$CenterOfMass_Y) * 0.9, 1)
max_y <- round(max(flight_long$CenterOfMass_Z) * 1.1, 1)
min_y <- round(min(flight_long$CenterOfMass_Z) * 0.9, 1)

# plot
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

# ---------------------------------------------------------
# now, I can calculate:
# - COM path (min to min? TO to TO? FC to FC?)
# - min/max joint angles & segment orientations
# - COM variables as a % of height
# ---------------------------------------------------------