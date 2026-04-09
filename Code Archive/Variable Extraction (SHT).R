
library(tidyverse)
library(readr)

# ------------------------
# read and clean the data
# ------------------------

d <- read_csv("Data/Split Data/MBB001_SHL_001.csv", col_names = FALSE)

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
   select(-NA.)

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

fc.mid20 <- fc.all %>%
   slice(start:(start + 19)) %>%
   mutate(`Contact No.` = row_number())

# first and last foot contacts
first_fc <- min(fc.mid20$FC_Index)
last_fc  <- max(fc.mid20$FC_Index)

# takeoff after each kept foot contact on the same plate
to.list <- to.all %>%
   group_by(FP) %>%
   summarise(TO_List = list(sort(TO_Index)), .groups = "drop")

pairs <- fc.mid20 %>%
   left_join(to.list, by = "FP") %>%
   mutate(TO_Index = map2_dbl(FC_Index, TO_List, ~ .y[.y > .x][1])) %>%
   select(FP, `Contact No.`, FC_Index, TO_Index) %>%
   filter(!is.na(TO_Index))

# build the 40-event table (starts FC_01, ends TO_20)
events <- bind_rows(pairs %>%
                       transmute(FP, Index = FC_Index, Event_ID = sprintf("FC_%02d", `Contact No.`)),
                    pairs %>%
                       transmute(FP, Index = TO_Index, Event_ID = sprintf("TO_%02d", `Contact No.`))) %>%
   arrange(Index)

# -------------------------------------------
# pull rows at each foot contact and takeoff
# -------------------------------------------

# event.data <- events %>%
#    mutate(Row = Index) %>%
#    left_join(d %>%
#                 mutate(Row = row_number()),
#              by = "Row") %>%
#    select(-Row, -HEIGHT_X, -MASS_X) %>%
#    rename(Time = TIME_0)

# ------------------------
# contact durations (×20)
# ------------------------

contact.times <- pairs %>%
   mutate(`Contact Time (ms)` = round((TO_Index - FC_Index) / sample.rate * 1000, 0),
          `Contact No.` = sprintf("FC%02d", `Contact No.`)) %>%
   select(`Contact No.`, `Contact Time (ms)`)


# --------------------------------------------------
# create a subset data frame for each contact cycle
# --------------------------------------------------

contact.list <- map(seq_len(nrow(pairs)), function(i) {
   d %>%
      slice(pairs$FC_Index[i]:pairs$TO_Index[i]) %>%
      mutate(Contact_No = sprintf("FC%02d", pairs$`Contact No.`[i]),
             Sample_in_contact = row_number())
   })

names(contact.list) <- sprintf("FC%02d", pairs$`Contact No.`)

# -----------------------------------------------------
# COP location (at peak GRF) within each contact cycle
# -----------------------------------------------------

cop.location <- tibble(`Contact No.` = names(contact.list),
   peak.grf   = map_dbl(contact.list, ~ max(pmax(.x$FP5_Z, .x$FP6_Z), na.rm = TRUE)),
   peak.plate = map_chr(contact.list, ~{
      i <- which.max(pmax(.x$FP5_Z, .x$FP6_Z))
      if (.x$FP5_Z[i] >= .x$FP6_Z[i]) "FP5" else "FP6"
   }),
   COPx = map_dbl(contact.list, ~{
      i <- which.max(pmax(.x$FP5_Z, .x$FP6_Z))
      if (.x$FP5_Z[i] >= .x$FP6_Z[i]) .x$FP5_COFP_X[i] else .x$FP6_COFP_X[i]
   }),
   COPy = map_dbl(contact.list, ~{
      i <- which.max(pmax(.x$FP5_Z, .x$FP6_Z))
      if (.x$FP5_Z[i] >= .x$FP6_Z[i]) .x$FP5_COFP_Y[i] else .x$FP6_COFP_Y[i]
   })
   )

# ---------------------------------
# compute a 95% confidence ellipse 
# ---------------------------------

ellipse.area <- bind_rows(
   cop.location %>%
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
   )

# plot COP locations at peak force during each contact cycle
ggplot(cop.location, aes(x = COPx, y = COPy, color = peak.plate)) +
   geom_point() +
   coord_flip() +
   labs(title = "Centre of Pressure at Peak GRF",
        x = "Anteroposterior",
        y = "Mediolateral") +
   theme_classic() +
   theme(plot.title = element_text(hjust = 0.5))

# plot COP paths during each contact cycle
contact_long <- bind_rows(contact.list, .id = "Contact_No")

contact_long_plot <- contact_long %>%
   select(Contact_No,
          FP5_COFP_X, FP5_COFP_Y,
          FP6_COFP_X, FP6_COFP_Y) %>%
   pivot_longer(cols = -Contact_No,
                names_to = c("Plate", ".value"),
                names_pattern = "(FP[56])_(.*)")

ggplot(contact_long_plot,
       aes(x = COFP_X,
           y = COFP_Y,
           group = interaction(Contact_No, Plate),
           colour = Plate)) +
   geom_path(alpha = 0.7) +
   coord_equal() +
   coord_flip() +
   labs(title = "COP Path",
        x = "Anteroposterior",
        y = "Mediolateral") +
   theme_classic() +
   theme(plot.title = element_text(hjust = 0.5))

# -------------------------------------
# 
# -------------------------------------

# ---------------------------------------------------------
# now, I can calculate:
# - COP displacement during each contact phase
# - COM displacement during each contact phase
# - COM-to-COP ratio
# ---------------------------------------------------------