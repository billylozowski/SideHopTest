
library(tidyverse)
library(readr)

# ------------------------
# read and clean the data
# ------------------------

d <- read_csv("Data/Split Data/WBB012_SHR_001.csv", col_names = FALSE)

# build column names
row.1 <- trimws(as.character(d[1, ]))
row.2 <- trimws(as.character(d[2, ]))
row.3 <- trimws(as.character(d[3, ]))
row.4 <- trimws(as.character(d[4, ]))

colnames(d) <- ifelse(row.2 == "COFP",
                      paste(row.1, row.2, row.4, sep = "_"), # row1_COFP_row4
                      paste(row.1, row.4, sep = "_"))        # row1_row4

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
   select(-1) # remove the first column (either "NA." or "NA_NA")

# ----------------------
# determine sample rate
# ----------------------

sample.rate <- d %>%
   filter(TIME == 1) %>% # just the row at 1s
   pull(FRAMES) - 1      # get the value in FRAMES and subtract 1

# --------------------------------------------
# foot contacts and take-offs for FP5 and FP6
# --------------------------------------------
# force threshold - set as desired
fc.threshold <- 20

# logical vector above fc_threshold
FP5.above <- d$FP5_Z >= fc.threshold
FP6.above <- d$FP6_Z >= fc.threshold

# FP5 crossings
FP5.fc.index <- which(diff(FP5.above) == 1)  + 1 # all instances where FP5_Z moves below → above fc.threshold
FP5.to.index <- which(diff(FP5.above) == -1) + 1 # all instances where FP5_Z moves above → below fc.threshold

# FP6 crossings
FP6.fc.index <- which(diff(FP6.above) == 1)  + 1 # all instances where FP5_Z moves below → above fc.threshold
FP6.to.index <- which(diff(FP6.above) == -1) + 1 # all instances where FP5_Z moves above → below fc.threshold

# diff() checks changes between consecutive logical samples (0, 0, 0, *1, 1, ... , 1, 1, *0, 0, 0, ...)
# 0  →  1 transition returns +1 (foot contact)
# 1  →  0 transition returns -1 (take-off)
# which() returns indices where transitions occur
# "+1" shifts to the first new-phase sample

# -------------------------------------------------------
# combine foot contacts and takeoffs across force plates
# -------------------------------------------------------
# foot contacts
fc.all <- bind_rows(tibble(FP = "FP5", FC_Index = FP5.fc.index),
                    tibble(FP = "FP6", FC_Index = FP6.fc.index)) %>%
   arrange(FC_Index)

# takeoffs
to.all <- bind_rows(tibble(FP = "FP5", TO_Index = FP5.to.index),
                    tibble(FP = "FP6", TO_Index = FP6.to.index)) %>%
   arrange(TO_Index)

# establish the middle 20 foot contacts
start <- floor((nrow(fc.all) - 20) / 2) + 1

# select just the middle 20 foot contacts
fc.mid20 <- fc.all %>%
   slice(start:(start + 19))

# first and last foot contacts
first_fc <- min(fc.mid20$FC_Index)
last_fc  <- max(fc.mid20$FC_Index)

# takeoff after each kept foot contact on the same plate
to.list <- to.all %>%
   group_by(FP) %>%
   summarise(TO_List = list(sort(TO_Index)), .groups = "drop")

# ----------------------------------------------------
# group each foot contact with its respective takeoff
# ----------------------------------------------------
# pair each foot contact with the takeoff immediately following it (contact cycle)
pairs <- fc.mid20 %>%
   left_join(to.list, by = "FP") %>%
   mutate(TO_Index = map2_dbl(FC_Index, TO_List, ~ .y[.y > .x][1])) %>%
   filter(!is.na(TO_Index)) %>%
   mutate(`Contact No.` = row_number(),                  # numeric row ID
          `Contact ID` = sprintf("%02d", `Contact No.`), # zero-padded row ID (i.e. 01, 02, ...)
          `Contact Time (ms)` = round((TO_Index - FC_Index) / sample.rate * 1000, 0)) %>%
   select(FP, `Contact No.`, `Contact ID`, FC_Index, TO_Index, `Contact Time (ms)`)

# --------------------------------------------------
# create a subset data frame for each contact cycle
# --------------------------------------------------

contact.list <- map(seq_len(nrow(pairs)), function(i) {
   d %>%
      slice(pairs$FC_Index[i]:pairs$TO_Index[i]) %>%
      mutate(Contact_No = pairs$`Contact No.`[i],
             Sample_in_contact = row_number())
   })

names(contact.list) <- sprintf("%02d", as.integer(pairs$`Contact No.`))

# -----------------------------------------------------
# COP location (at peak GRF) within each contact cycle
# -----------------------------------------------------

cop.location <- tibble(`Contact No.` = as.integer(names(contact.list))) %>%    # initialise output with contact labels
   mutate(res = map(contact.list, ~{                                           # iterate over each contact cycle
      
      i  <- which.max(pmax(.x$FP5_Z, .x$FP6_Z))                                # index of peak vertical GRF (either plate)
      fp <- if (.x$FP5_Z[i] >= .x$FP6_Z[i]) "FP5" else "FP6"                   # determine which plate produced peak
      
      tibble(                                                                  # return one-row tibble per contact
         peak.grf   = pmax(.x$FP5_Z[i], .x$FP6_Z[i]),                          # peak GRF magnitude
         peak.plate = fp,                                                      # plate at peak GRF
         COPx       = if (fp == "FP5") .x$FP5_COFP_X[i] else .x$FP6_COFP_X[i], # COP x at peak
         COPy       = if (fp == "FP5") .x$FP5_COFP_Y[i] else .x$FP6_COFP_Y[i]) # COP y at peak
      
   })) %>%
   unnest(res)                                                                # expand list-column into standard columns

# ------------------------------------------------
# create a data frame for all extracted variables
# ------------------------------------------------

final.data <- pairs %>%
   left_join(cop.location, by = "Contact No.") %>%
   select(-peak.plate)

# ---------------------------------
# compute a 95% confidence ellipse 
# ---------------------------------

ellipse.area <- bind_rows(cop.location %>%
      # convert COP values to centimetres
      mutate(COPx_cm = COPx * 100,
             COPy_cm = COPy * 100) %>%
      group_by(peak.plate) %>%
      # contact area per plate
      summarise(`Contact Area (cm^2)` = pi * qchisq(0.95, 2) *
                   sqrt(det(cov(cbind(COPx_cm, COPy_cm)))),
                .groups = "drop"),
   cop.location %>%
      mutate(COPx_cm = COPx * 100,
             COPy_cm = COPy * 100) %>%
      # total contact area
      summarise(peak.plate = "Total",
                `Contact Area (cm^2)` = pi * qchisq(0.95, 2) *
                   sqrt(det(cov(cbind(COPx_cm, COPy_cm)))))
   ) %>%
   rename(FP = peak.plate)

# though "Total" is a combination of FP5 and FP6, it's being included in column 
# "FP" for ease of reading.

# --------------------------------
# visualise the COP location data
# --------------------------------
# COP locations during each contact cycle
ggplot(cop.location, 
       aes(x = COPy, 
           y = COPx, 
           color = peak.plate)) +
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
contact_long <- bind_rows(contact.list, .id = "Contact No.")

contact_long_plot <- contact_long %>%
   select(`Contact No.`,
          FP5_COFP_X, FP5_COFP_Y,
          FP6_COFP_X, FP6_COFP_Y) %>%
   pivot_longer(cols          = -`Contact No.`,
                names_to      = c("Plate", ".value"),
                names_pattern = "(FP[56])_(.*)")

# build the plot
ggplot(contact_long_plot,
       aes(x = COFP_Y,                    # ml direction
           y = COFP_X,                    # ap direction
           group  = interaction(`Contact No.`, Plate),
           colour = Plate)) +
   geom_path(alpha = 0.7) +
   coord_equal() +
   scale_y_reverse() +                    # invert y axis
   scale_x_continuous(position = "top") + # move x axis to top of plot
   labs(title = "COP Path",
        x     = "Mediolateral",
        y     = "Anteroposterior") +
   theme_classic() +
   theme(aspect.ratio = 0.5)              # 2:1 width vs. height ratio

# -----------------------------------------------------
# calculate COP displacement during each contact cycle
# -----------------------------------------------------
# create a function for computing the triangulated distance between rows
path_length <- function(x, y) {
   dx <- diff(x)
   dy <- diff(y)
   sum(sqrt(dx^2 + dy^2), na.rm = TRUE)
}

# calculate COP and COM variables during each contact cycle
cop.displacement <- pairs %>%
   mutate(
      # COP path (xy) in centimetres - active plate only
      `COP Path (cm)` = map2_dbl(FP, `Contact No.`, ~{
         df <- contact.list[[sprintf("%02d", .y)]]
         if (.x == "FP5") 
            path_length(df$FP5_COFP_X, df$FP5_COFP_Y) * 100
         else
            path_length(df$FP6_COFP_X, df$FP6_COFP_Y) * 100
         }),
      
      # COM path (xy) in centimetres
      `COM Path (cm)` = map2_dbl(FP, `Contact No.`, ~{       
         df <- contact.list[[sprintf("%02d", .y)]]
         if (.x == "FP5") 
            path_length(df$CenterOfMass_X, df$CenterOfMass_Y) * 100
         else
            path_length(df$CenterOfMass_X, df$CenterOfMass_Y) * 100
         }),
      
      # COP ML range (y) in centimetres
      `COP_ML (cm)` = map2_dbl(FP, `Contact No.`, ~{       
         df <- contact.list[[sprintf("%02d", .y)]]
         if (.x == "FP5") 
            (max(df$FP5_COFP_Y, na.rm = TRUE) - min(df$FP5_COFP_Y, na.rm = TRUE)) * 100
         else
            (max(df$FP6_COFP_Y, na.rm = TRUE) - min(df$FP6_COFP_Y, na.rm = TRUE)) * 100
      }),
      
      # COP AP range (x) in centimetres
      `COP_AP (cm)` = map2_dbl(FP, `Contact No.`, ~{       
         df <- contact.list[[sprintf("%02d", .y)]]
         if (.x == "FP5") 
            (max(df$FP5_COFP_X, na.rm = TRUE) - min(df$FP5_COFP_X, na.rm = TRUE)) * 100 
         else
            (max(df$FP6_COFP_X, na.rm = TRUE) - min(df$FP6_COFP_X, na.rm = TRUE)) * 100
      }),
      
      # COM ML range (y) in centimetres
      `COM_ML (cm)` = map2_dbl(FP, `Contact No.`, ~{       
         df <- contact.list[[sprintf("%02d", .y)]]
         if (.x == "FP5") 
            (max(df$CenterOfMass_Y, na.rm = TRUE) - min(df$CenterOfMass_Y, na.rm = TRUE)) * 100
         else
            (max(df$CenterOfMass_Y, na.rm = TRUE) - min(df$CenterOfMass_Y, na.rm = TRUE)) * 100
      }),
      
      # COM AP range (x) in centimetres
      `COM_AP (cm)` = map2_dbl(FP, `Contact No.`, ~{       
         df <- contact.list[[sprintf("%02d", .y)]]
         if (.x == "FP5") 
            (max(df$CenterOfMass_X, na.rm = TRUE) - min(df$CenterOfMass_X, na.rm = TRUE)) * 100 
         else
            (max(df$CenterOfMass_X, na.rm = TRUE) - min(df$CenterOfMass_X, na.rm = TRUE)) * 100
      }),
      
      # COM vertical range (z) in centimetres
      `COM_Z (cm)` = map2_dbl(FP, `Contact No.`, ~{       
         df <- contact.list[[sprintf("%02d", .y)]]
         if (.x == "FP5") 
            (max(df$CenterOfMass_Z, na.rm = TRUE) - min(df$CenterOfMass_Z, na.rm = TRUE)) * 100 
         else
            (max(df$CenterOfMass_Z, na.rm = TRUE) - min(df$CenterOfMass_Z, na.rm = TRUE)) * 100
      })
      ) %>%
   select(`Contact No.`, FP, `COP Path (cm)`, `COM Path (cm)`,
          `COP_ML (cm)`, `COP_AP (cm)`, `COM_ML (cm)`, `COM_AP (cm)`, `COM_Z (cm)`)

# join COP/COM data to final.data
final.data <- final.data %>%
   left_join(cop.displacement, by = "Contact No.") %>%
   rename(FP = FP.x) %>%
   select(-FP.y) %>%
   # calculate ratio of COM path to COP path
   mutate(`COM:COP Path Ratio` = `COM Path (cm)` / `COP Path (cm)`)

# ---------------------------------------------------------
# now, I can calculate:
# - COP displacement during each contact phase
# - COM displacement during each contact phase
# - COM-to-COP ratio
# ---------------------------------------------------------