
# this script identifies candidate predictor variables for side-hop test (SHT) performance
# to do this we use the following workflow:
# 
# 1. With SHT performance as the outcome (mean of right & left), use "glmulti" to determine 
# the best linear model for the data.
#  - "glmulti" creates every possible model to ensure the best fit for the data. I doesn't
#    over- or underfit the model, and its diagnostics allow us to see how many models (as 
#    a % each variable is included in).
#    
# 2. With a streamlined set of predictor variables identified, run the same glmulti model,
# but this time with interactions.
#  - mainly to see if any of the predictor relationships are dependent on CAI classifiation.
# 
# 3. Summarise the final model(s) to establish which variables, if any, best predict SHT
# performance.
# 
# BENEFITS of predicting SHT performance?
# - SHT amongst one of the best FMS tests for CAI populations.
# - Could provide links between FMS and other screening tests, allowing practitioners to 
#   catch increases in injury risk earlier (or at least confirm them more confidently).
#   

# ------------------------------------------------------------------------------
# core packages
# ------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(performance)
library(DHARMa)

# -------------------------
# load data
# -------------------------

data.path <- "Data/SHT Data_biomech.xlsx"

data.study2 <- read_excel(data.path, sheet = "Study 2") %>%
   mutate(`Study no.` = 2)

# collapse spaces in column names
names(data.study2) <- gsub(" ", "_", names(data.study2))

# raw phase data
contact   <- read_csv("Data/Processed/SHT - Contact Phase.csv") %>% mutate(Code = substr(Trial_Name, 1, 6))
flight    <- read_csv("Data/Processed/SHT - Flight Phase.csv")  %>% mutate(Code = substr(Trial_Name, 1, 6))
hop.cycle <- read_csv("Data/Processed/SHT - Hop Cycle.csv")     %>% mutate(Code = substr(Trial_Name, 1, 6))

# -------------------------------------------
# combine base data with kinematic variables
# -------------------------------------------
# Contact Phase:
#    - Landing variability (COP location area)?

# select variables from demographics
data.study2.select <- data.study2 %>%
   select(Code, Height_m, Mass_kg, Age, SHT_R, SHT_L, CAI)

# select variables from contact phase
contact.select <- contact %>%
   group_by(Code, Trial_Name) %>%
   summarise(L_Ankle_DF_range_contact     = max(l_ankle_df_X_max_contact, na.rm = TRUE)    - min(l_ankle_df_X_min_contact, na.rm = TRUE),
             R_Ankle_DF_range_contact     = max(r_ankle_df_X_max_contact, na.rm = TRUE)    - min(r_ankle_df_X_min_contact, na.rm = TRUE),
             L_Knee_Flexion_range_contact = max(l_knee_fle_X_max_contact, na.rm = TRUE)    - min(l_knee_fle_X_min_contact, na.rm = TRUE),
             R_Knee_Flexion_range_contact = max(r_kne_fle_X_max_contact, na.rm = TRUE)     - min(r_kne_fle_X_min_contact, na.rm = TRUE),
             L_Hip_Flexion_range_contact  = max(l_hip_fle_X_max_contact, na.rm = TRUE)     - min(l_hip_fle_X_min_contact, na.rm = TRUE),
             R_Hip_Flexion_range_contact  = max(r_hip_fle_X_max_contact, na.rm = TRUE)     - min(r_hip_fle_X_min_contact, na.rm = TRUE),
             L_Hip_Abd_range_contact      = max(l_hip_abd_Y_max_contact, na.rm = TRUE)     - min(l_hip_abd_Y_min_contact, na.rm = TRUE),
             R_Hip_Abd_range_contact      = max(r_hip_abd_Y_max_contact, na.rm = TRUE)     - min(r_hip_abd_Y_min_contact, na.rm = TRUE),
             Trunk_Tilt_AP_contact        = max(trunk_tilt_ap_X_max_contact, na.rm = TRUE) - min(trunk_tilt_ap_X_min_contact, na.rm = TRUE),
             Trunk_Tilt_ML_contact        = max(trunk_tilt_ml_Y_max_contact, na.rm = TRUE) - min(trunk_tilt_ml_Y_min_contact, na.rm = TRUE),
             Contact_Time_ms              = mean(Contact_Time_ms, na.rm = TRUE),
             COP_Depth                    = max(COPx, na.rm = TRUE) - min(COPx, na.rm = TRUE),
             COP_Width                    = max(COPy, na.rm = TRUE) - min(COPy, na.rm = TRUE),
             peak_grf                     = mean(peak_grf, na.rm = TRUE),
             .groups = "drop") %>%
   select(Trial_Name, Code, Contact_Time_ms,
          COP_Depth, COP_Width, peak_grf, 
          L_Ankle_DF_range_contact, R_Ankle_DF_range_contact, 
          L_Knee_Flexion_range_contact, R_Knee_Flexion_range_contact,
          L_Hip_Flexion_range_contact, R_Hip_Flexion_range_contact,
          L_Hip_Abd_range_contact, R_Hip_Abd_range_contact,
          Trunk_Tilt_AP_contact, Trunk_Tilt_ML_contact)

# select variables from flight phase
flight.select <- flight %>%
   group_by(Code, Trial_Name) %>%
   summarise(Flight_Time_ms = mean(Flight_Time_ms, na.rm = TRUE),
             .groups = "drop")

# select variables across entire hop cycle
hop.cycle.select <- hop.cycle %>%
   group_by(Code, Trial_Name) %>%
   summarise(COM_Z_Excursion_hop      = max(CenterOfMass_Z_max_hop, na.rm = TRUE)  - min(CenterOfMass_Z_min_hop, na.rm = TRUE),
             L_Ankle_DF_range_hop     = max(l_ankle_df_X_max_hop, na.rm = TRUE)    - min(l_ankle_df_X_min_hop, na.rm = TRUE),
             R_Ankle_DF_range_hop     = max(r_ankle_df_X_max_hop, na.rm = TRUE)    - min(r_ankle_df_X_min_hop, na.rm = TRUE),
             L_Knee_Flexion_range_hop = max(l_knee_fle_X_max_hop, na.rm = TRUE)    - min(l_knee_fle_X_min_hop, na.rm = TRUE),
             R_Knee_Flexion_range_hop = max(r_kne_fle_X_max_hop, na.rm = TRUE)     - min(r_kne_fle_X_min_hop, na.rm = TRUE),
             L_Hip_Flexion_range_hop  = max(l_hip_fle_X_max_hop, na.rm = TRUE)     - min(l_hip_fle_X_min_hop, na.rm = TRUE),
             R_Hip_Flexion_range_hop  = max(r_hip_fle_X_max_hop, na.rm = TRUE)     - min(r_hip_fle_X_min_hop, na.rm = TRUE),
             L_Hip_Abd_range_hop      = max(l_hip_abd_Y_max_hop, na.rm = TRUE)     - min(l_hip_abd_Y_min_hop, na.rm = TRUE),
             R_Hip_Abd_range_hop      = max(r_hip_abd_Y_max_hop, na.rm = TRUE)     - min(r_hip_abd_Y_min_hop, na.rm = TRUE),
             Trunk_Tilt_AP_hop        = max(trunk_tilt_ap_X_max_hop, na.rm = TRUE) - min(trunk_tilt_ap_X_min_hop, na.rm = TRUE),
             Trunk_Tilt_ML_hop        = max(trunk_tilt_ml_Y_max_hop, na.rm = TRUE) - min(trunk_tilt_ml_Y_min_hop, na.rm = TRUE),
             .groups = "drop")

# ----------------------------------------------------
# join select data frames together and save as a .csv
# ----------------------------------------------------
sht.select <- contact.select %>%
   left_join(flight.select,      by = c("Code", "Trial_Name")) %>%
   left_join(hop.cycle.select,   by = c("Code", "Trial_Name")) %>%
   left_join(data.study2.select, by = "Code") %>%
   select(Trial_Name, Code, Height_m, Mass_kg, Age, SHT_R, SHT_L, CAI, everything()) 

   write.csv(file = "Data/Processed/SHT - Study 2 - Kinematic Predictors (both trials).csv", row.names = FALSE)

# ----------------------------------------------
# person-level summary (one row per individual)
# ----------------------------------------------
trial.summary.data <- sht.select %>%
   mutate(Side = case_when(str_detect(Trial_Name, "SHL") ~ "L",
                           str_detect(Trial_Name, "SHR") ~ "R"),
          
          SHT_time = case_when(Side == "L" ~ SHT_L,
                               Side == "R" ~ SHT_R),
          
          Ankle_DF_range_contact = case_when(Side == "L" ~ L_Ankle_DF_range_contact,
                                             Side == "R" ~ R_Ankle_DF_range_contact),
          
          Knee_Flexion_range_contact = case_when(Side == "L" ~ L_Knee_Flexion_range_contact,
                                                 Side == "R" ~ R_Knee_Flexion_range_contact),
          
          Hip_Flexion_range_contact = case_when(Side == "L" ~ L_Hip_Flexion_range_contact,
                                                Side == "R" ~ R_Hip_Flexion_range_contact),
          Hip_Abd_range_contact = case_when(Side == "L" ~ L_Hip_Abd_range_contact,
                                            Side == "R" ~ R_Hip_Abd_range_contact),
          
          Ankle_DF_range_hop = case_when(Side == "L" ~ L_Ankle_DF_range_hop,
                                         Side == "R" ~ R_Ankle_DF_range_hop),
          
          Knee_Flexion_range_hop = case_when(Side == "L" ~ L_Knee_Flexion_range_hop,
                                             Side == "R" ~ R_Knee_Flexion_range_hop),
          
          Hip_Flexion_range_hop = case_when(Side == "L" ~ L_Hip_Flexion_range_hop,
                                            Side == "R" ~ R_Hip_Flexion_range_hop),
          
          Hip_Abd_range_hop = case_when(Side == "L" ~ L_Hip_Abd_range_hop,
                                        Side == "R" ~ R_Hip_Abd_range_hop)) %>%
   
   select(Code, CAI,
          Height_m, Mass_kg, Age,
          SHT_time,
          Contact_Time_ms, Flight_Time_ms,
          COM_Z_Excursion_hop,
          COP_Depth, COP_Width, peak_grf,
          Ankle_DF_range_contact, Knee_Flexion_range_contact,
          Hip_Flexion_range_contact, Hip_Abd_range_contact,
          Ankle_DF_range_hop, Knee_Flexion_range_hop,
          Hip_Flexion_range_hop, Hip_Abd_range_hop,
          Trunk_Tilt_AP_contact, Trunk_Tilt_ML_contact,
          Trunk_Tilt_AP_hop, Trunk_Tilt_ML_hop) %>%
   group_by(Code) %>%
   summarise(CAI = first(CAI),
             n_trials   = n(),
             across(Height_m:Trunk_Tilt_ML_hop, mean, na.rm = TRUE),
             .groups = "drop")

write_csv(trial.summary.data, file = "Data/Processed/SHT - Study 2 - Kinematic Predictors (summarised across trials).csv")

# --------------------------
# screen and model the data (use "model.data")
# --------------------------

model.data <- read_csv("Data/Processed/SHT - Study 2 - Kinematic Predictors (summarised across trials).csv") %>%
   mutate(CAI = as.factor(CAI))

# ------------------------------------------------------------------------------
# random forest model - SHT performance
# ------------------------------------------------------------------------------
library(party)
library(flextable)

# construct a random forest model
model.data.rf <- model.data %>%
   select(Height_m:last_col()) 

model.data.rf <- cforest(SHT_time ~ ., data = model.data.rf)

# extract Variable Importance
model.data.rf.vi <- as.data.frame(varimp(model.data.rf)) %>%
   rownames_to_column("Variable") %>%
   rename("V.I." = "varimp(model.data.rf)") %>%
   arrange(desc(`V.I.`)) %>% # ensure highest ranked variables are on top
   mutate(Variable = factor(Variable, levels = Variable))

# variable importance table
vi.table <- model.data.rf.vi %>%
   flextable() %>%
   autofit()

save_as_docx(vi.table, 
             path = "Tables & Figures/SHT - Study 2 - Kinematic Predictors of SHT Performance/SHT - Ranked Variable Importance for SHT Performance.docx")

# top 10 predictors
model.data.rf.vi <- head(model.data.rf.vi, 10)

# visualise V.I. for clearer picture
p.sht.rf.vi <- ggplot(model.data.rf.vi, aes(Variable, V.I., fill = V.I.)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(x     = NULL,
        y     = "V.I.",
        title = "Ranked Variable Importance") +
   theme(axis.text.x  = element_blank(),
         axis.ticks.x = element_blank(),
         plot.title   = element_text(hjust = 0.5)) +
   scale_y_continuous(expand = expansion(mult = c(0, 0.25)))

ggsave(p.sht.rf.vi, 
       file = "Tables & Figures/SHT - Study 2 - Kinematic Predictors of SHT Performance/SHT - Ranked Variable Importance for SHT Performance.png",
       dpi = 600)

# 5 variables are clearly much more closely associated with SHT performance (time to completion)

#                      Variable        V.I.
# 1              Flight_Time_ms 0.231237749
# 2             Contact_Time_ms 0.181362765
# 3         COM_Z_Excursion_hop 0.155177011
# 4          Ankle_DF_range_hop 0.144420336
# 5      Ankle_DF_range_contact 0.135569813

# Consequently, these will be passed through the glmulti portion of the analysis!

# -----------------
# glmulti analysis
# -----------------

# I've had issues loading "rJava", so this is a hard-coded file path the the files required
Sys.setenv(JAVA_HOME = "C:/Program Files/Eclipse Adoptium/jdk-25.0.2.10-hotspot")

library(rJava)
library(glmulti)

# identify the top n performing variables from the random forest model
model.variables <- model.data.rf.vi$Variable[1:5] # modify this as required

# subset model.data to only include top performing variables, plus "Height_m" and "Mass_kg"
# height and mass included given their modulating effect from study 1
model.data.rf.subset <- model.data %>%
   select(SHT_time, Height_m, Mass_kg, any_of(model.variables))

# predictor names only
predictor.names <- setdiff(names(model.data.rf.subset), "SHT_time")

# safer formula construction
model.formula <- reformulate(predictor.names, response = "SHT_time")

glmulti(model.formula,
        data        = model.data.rf.subset,
        crit        = bic,   # model fit criterion - good for both large and small sample sizes
        level       = 2,     # 1 without interactions, 2 with interactions
        method      = "d",   # computes number of all possible models (doesn't execute)
        # other options include "h" (exhaustive) or "g" (genetic)
        family      = gaussian,
        fitfunction = glm,   # specify the model type (lm or glm)
        confsetsize = 100)   # keep the 100 best models (confidence set)

# without interactions, the candidate set contains 128 models.
# with interactions, candidate set contains 268435456 models!!!
# for SHT performance, we will exclude any interactions.

model.sht.glmulti.g <- glmulti(model.formula,
                               data        = model.data.rf.subset,
                               crit        = bic,   # model fit criterion - good for both large and small sample sizes
                               level       = 1,     # 1 without interactions, 2 with interactions
                               method      = "g",   
                               family      = gaussian,
                               fitfunction = glm,   # specify the model type (lm or glm)
                               confsetsize = 100)   # keep the 100 best models (confidence set)

# inspect candidate terms
attr(terms(model.formula), "term.labels")

# model results
print(model.sht.glmulti.g)

# plot the best 100 models base on information criterion (IC)
plot.ic <- plot(model.sht.glmulti.g) # red line indicates 2 IC units from the best model

ggsave(plot.ic,
       file = "Tables & Figures/SHT - Study 2 - Kinematic Predictors of SHT Performance/SHT - Study 2 - IC Profile for SHT Performance.png",
       dpi = 600)


# model-averaged importance terms
# --------------------------------
# specify output parameters
png(filename = "Tables & Figures/SHT - Study 2 - Kinematic Predictors of SHT Performance/SHT - Study 2 - Model Importance Terms for SHT Performance.png",
    width = 2000, height = 1200, res = 150)

# plot margins
par(mar = c(5, 12, 1, 2) + 0.1)

# plot the model
plot(model.sht.glmulti.g, type = "s")

# close the .png window
dev.off()


# summarise the best 3 models
# ----------------------------

glmulti.table <- weightable(model.sht.glmulti.g)[1:3,] %>% # select the best models
   regulartable() %>%
   autofit()

save_as_docx(glmulti.table, 
             path = "Tables & Figures/SHT - Study 2 - Kinematic Predictors of SHT Performance/SHT - Top 3 Models for SHT Performance.docx")


# based on the model-averaged importance terms, the top 2 variables are:
# 
# - "Flight_Time_ms"
# - "COM_Z_Excursion_hop"
# 
# there is then a cluster of 3 additional variable following these:
# 
# - "Ankle_DF_range_contact"
# - "Contact_Time_ms" 
# - "Ankle_DF_range_hop" 
# 
# "Height_m" and "Mass_kg" are the two worst performing predictors in this subset!

# --------------------
# build linear models
# --------------------

# linear model (without interaction terms for height and mass)
m1 <- lm(SHT_time ~ 
            # predictors
            Flight_Time_ms + COM_Z_Excursion_hop + Ankle_DF_range_contact,
         data = model.data)

# linear model (with interaction terms for height and mass)
m2 <- lm(SHT_time ~ 
            # predictors
            Flight_Time_ms + COM_Z_Excursion_hop +
            # interaction terms
            Flight_Time_ms*Height_m + COM_Z_Excursion_hop*Height_m +
            Flight_Time_ms*Mass_kg  + COM_Z_Excursion_hop*Mass_kg,
         data = model.data)

anova(m1, m2)

# model m1a (with ankle df)

# create a table of the final model outputs, and save as a word .doc
# -------------------------------------------------------------------

library(sjPlot)
tab_model(m1,
          show.df   = FALSE,
          show.se   = TRUE, 
          string.se = "SE",
          digits    = 3,
#          dv.labels = c("No Interactions", "Interactions"),
          file = "Tables & Figures/SHT - Study 2 - Kinematic Predictors of SHT Performance/SHT - Study 2 - Model Summaries.doc")


# summarise the strongest model
# ------------------------------

library(car)
library(effectsize)

summary(m1)           # estimates and coefficients
model_performance(m1) # fit statistics

anova_m1 <- car::Anova(m1, type = 2)                             # model evaluation (for effect sizes)
eta_squared(anova_m1, partial = TRUE, alternative = "two.sided") # effect sizes

