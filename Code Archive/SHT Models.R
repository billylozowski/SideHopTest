
# this script identifies candidate predictor variables for side-hop test (SHT) performance
# to do this we use the following workflow:
# 
# 1. With SHT performance as the outcome (composite, then right & left separately), run 
# a random-forest (RF) model and select the most appropriate variables based on variable 
# importance (VI).
# 
# 2. Use "glmulti" with SHT performance as the outcome to determine the best linear model
# for the data.
#  - "glmulti" creates every possible model to ensure the best fit for the data. I doesn't
#    over- or underfit the model, and its diagnostics allow us to see how many models (as 
#    a % each variable is included in).
# 
# 3. Summarise the final model(s) to establish whic variables, if any, best predict SHT
# performance.
# 
# BENEFITS of predicting SHT performance?
# - SHT amongst one of the best FMS tests for CAI populations.
# - Could provide links between FMS and other screening tests, allowing practitioners to 
#   catch increases in injury risk earlier (or at least confirm them more confidently).

# ------------------------------------------------------------------------------
# core packages
# ------------------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(performance)

# ------------------------------------------------------------------------------
# load the data
# ------------------------------------------------------------------------------
d.sht <- read_excel("SHT Data.xlsx", sheet = "Nm")

d.sht <- d.sht %>%
   mutate(`Side Hop Agg.` = rowMeans(cbind(`Side Hop R`, `Side Hop L`), na.rm = TRUE)) %>%
   relocate(`Side Hop Agg.`, .after = `Side Hop L`) %>%
   rename(Height = `Height (m)`)

# collapse spaces in column names
names(d.sht) <- gsub(" ", "_", names(d.sht))

# ------------------------------------------------------------------------------
# remove any outliers
# ------------------------------------------------------------------------------
# loop over each numeric column in a data frame
remove_outliers <- function(data) {
   data[] <- lapply(data, function(column) {
      if (is.numeric(column)) {
         
         # calculate Q1, Q3, and IQR
         Q1 <- quantile(column, 0.25, na.rm = TRUE)
         Q3 <- quantile(column, 0.75, na.rm = TRUE)
         IQR_val <- Q3 - Q1
         
         # set thresholds for outliers
         lower_limit <- Q1 - 2 * IQR_val
         upper_limit <- Q3 + 2 * IQR_val
         
         # replace outliers with NA (or you can choose to filter them out)
         column[column < lower_limit | column > upper_limit] <- NA
      }
      return(column)
   })
   return(data)
}

d.sht.clean <- remove_outliers(d.sht) 

# ------------------------------------------------------------------------------
# random forest model - aggregated SHT performance
# ------------------------------------------------------------------------------
library(party)

# remove R or L columns for the aggregate analysis
d.sht.numeric <- d.sht.clean %>%
   select(Height:ncol(d.sht.clean), 
          -Side_Hop_R, -Side_Hop_L)

# construct the rf model for the aggregated SHT score
rf.sht.agg <- cforest(Side_Hop_Agg. ~ ., data = d.sht.numeric)

# extract Variable Importance
variable.importance.sht.agg <- as.data.frame(varimp(rf.sht.agg)) %>%
   rownames_to_column("Variable") %>%
   rename("V.I." = "varimp(rf.sht.agg)") %>%
   arrange(desc(`V.I.`)) %>% # ensure highest ranked variables are on top
   mutate(Variable = factor(Variable, levels = Variable))

# top 10 predictors
sht.agg.VI <- head(variable.importance.sht.agg, 10)

# visualise V.I. for clearer picture
p.agg <- ggplot(sht.agg.VI, aes(Variable, V.I., fill = V.I.)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(x = NULL,
        y = "Variable Importance") +
   theme(axis.text.x = element_blank(),
         x.ticks = NULL) +
   scale_y_continuous(expand = expansion(mult = c(0, 0.25)))

ggsave(p.agg, file = "Variable Importance for Aggregated SHT Performance.png",
       dpi = 600)

# it appears that only 3-5 variables are necessary in the aggregated model

# ------------------------------------------------------------------------------
# random forest model - SHT R performance
# ------------------------------------------------------------------------------
# remove any non-numeric columns & retain only right limb
d.sht.numeric.R <- d.sht.clean %>%
   select(Height:ncol(d.sht.clean),
          -Side_Hop_Agg.,
          -matches("LExt|LFle")) # remove left limb data

# construct the rf model for the right SHT score
rf.sht.r <- cforest(Side_Hop_R ~ ., data = d.sht.numeric.R)

# extract Variable Importance
variable.importance.sht.r <- as.data.frame(varimp(rf.sht.r)) %>%
   rownames_to_column("Variable") %>%
   rename("V.I." = "varimp(rf.sht.r)") %>%
   arrange(desc(`V.I.`)) %>% # ensure highest ranked variables are on top
   mutate(Variable = factor(Variable, levels = Variable))

# top 10 predictors
sht.r.VI <- head(variable.importance.sht.r, 10)

# return the variable names
names.right <- sht.r.VI$Variable 

# visualise V.I. for clearer picture
p.r <- ggplot(sht.r.VI, aes(Variable, V.I., fill = V.I.)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(x = NULL,
        y = "Variable Importance") +
   theme(axis.text.x = element_blank(),
         x.ticks = NULL) +
   scale_y_continuous(expand = expansion(mult = c(0, 0.25)))

ggsave(p.r, file = "Variable Importance for Right SHT Performance.png",
       dpi = 600)

# the right model only appears to require 4 variables

# ------------------------------------------------------------------------------
# random forest model - SHT L performance
# ------------------------------------------------------------------------------
# remove any non-numeric columns & retain only left limb
d.sht.numeric.L <- d.sht.clean %>%
   select(Height:ncol(d.sht.clean),
          -Side_Hop_Agg.,
          -matches("RExt|RFle")) # remove right limb data

# construct the rf model for the left SHT score
rf.sht.l <- cforest(Side_Hop_L ~ ., data = d.sht.numeric.L)

# extract Variable Importance
variable.importance.sht.l <- as.data.frame(varimp(rf.sht.l)) %>%
   rownames_to_column("Variable") %>%
   rename("V.I." = "varimp(rf.sht.l)") %>%
   arrange(desc(`V.I.`)) %>% # ensure highest ranked variables are on top
   mutate(Variable = factor(Variable, levels = Variable))

# top 10 predictors
sht.l.VI <- head(variable.importance.sht.l, 10)

# return the variable names
names.left <- sht.l.VI$Variable 

# visualise V.I. for clearer picture
p.l <- ggplot(sht.l.VI, aes(Variable, V.I., fill = V.I.)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(x = NULL,
        y = "Variable Importance") +
   theme(axis.text.x = element_blank(),
         x.ticks = NULL) +
   scale_y_continuous(expand = expansion(mult = c(0, 0.25)))

ggsave(p.l, file = "Variable Importance for Left SHT Performance.png",
       dpi = 600)

# likely that the top 5 variables should be sufficient for the glmulti analysis

# ------------------------------------------------------------------------------
# save top 10 VI for each random forest
# ------------------------------------------------------------------------------
library(writexl)

write_xlsx(list(
   `Aggregated VI` = sht.agg.VI,
   `Right VI`      = sht.r.VI,
   `Left VI`       = sht.l.VI),
   path = "Variable Importance for 3 RF Models.xlsx")

# ------------------------------------------------------------------------------
# "glmulti" models with top 10 VI variables from each Random Forest
# ------------------------------------------------------------------------------
library(glmulti)
library(rJava)

################################################################################
# aggregated model
################################################################################
vars.to.keep <- sht.agg.VI$Variable[1:5]

aggregate.subset <- d.sht.clean %>%
   select(Side_Hop_Agg., any_of(vars.to.keep))

# return the variable names
names.agg <- colnames(aggregate.subset)

glmulti(as.formula(paste("Side_Hop_Agg. ~", 
                         # dynamically input top 10 variable names
                         paste(names.agg, collapse = " + "))), 
        data = aggregate.subset,
        crit = bic,          # model fit criterion - good for both large and small sample sizes
        level = 1,           # 1 without interactions, 2 with interactions
        method = "d",        # computes number of all possible models (doesn't execute)
                             # other options include "h" (exhaustive) or "g" (genetic)
        family = gaussian,
        fitfunction = glm,   # specify the model type (lm or glm)
        confsetsize = 100)   # keep the 100 best models (confidence set)

# without interactions, the candidate set contains 64 models.
# with interactions, candidate set contains 2097152 models!!!
# for aggregated SHT performance, we will exclude any interactions.

sht.agg.model.g <- glmulti(as.formula(paste("Side_Hop_Agg. ~", 
                                    # dynamically input variable names from subset 2
                                    paste(names.agg, collapse = " + "))),
                   data = aggregate.subset,
                   crit = bic,         # model fit criterion
                   level = 1,          # 1 without interactions, 2 with interactions
                   method = "g",       # exhaustive screening algorithm (better with fewer predictors)
                   family = gaussian,
                   fitfunction = glm,  # specify the model type (lm or glm)
                   confsetsize = 100)  # keep the 100 best models (confidence set)

# model results
print(sht.agg.model.g)

# plot the best 100 models base on information criterion (IC)
plot(sht.agg.model.g) # red line indicates 2 IC units from the best model

################################################################################
# Model-Averaged Importance Terms
################################################################################
# specify output parameters
png(filename = "Model Importance Terms (aggregated SHT).png",
    width = 2000, height = 1200, res = 150)

# plot margins
par(mar = c(5, 12, 1, 2) + 0.1)

# plot the model
plot(sht.agg.model.g, type = "s")

# close the .png window
dev.off()


################################################################################
# summarise the best 3 models
################################################################################
library(flextable)

weightable(sht.agg.model.g)[1:3,] %>% # select the best models
   flextable::regulartable() %>%
   flextable::autofit()

# final model
m0 <- lm(Side_Hop_Agg. ~
            # predictors
            HU_workPRCV_LFle_180 +
            HU_timePeakTorqueHeld_RFle_180 +
            HU_jointAnglePeakTorque_LFle_60 +
            HU_workPR_RFle_60,
         data = d.sht.clean)

# create a table of the final model outputs, and save as a word .doc
library(sjPlot)
tab_model(m0, 
          show.std  = TRUE,
          show.df   = FALSE,
          show.se   = TRUE, 
          string.se = "SE",
          dv.labels = c("Aggregated SHT Model"),
          file = "Aggregated SHT Models.doc")

# summarise the strongest model
summary(m0)
model_performance(m0)

################################################################################
# right model
################################################################################
vars.to.keep <- sht.r.VI$Variable[1:4]

right.subset <- d.sht.clean %>%
   select(Side_Hop_R, any_of(vars.to.keep))

names.right <- colnames(right.subset)

# #####################################
# this is running more than 32 models!!
# #####################################
glmulti(as.formula(paste("Side_Hop_R ~", 
                         # dynamically input top 10 variable names
                         paste(names.right, collapse = " + "))), 
        data = right.subset,
        crit = bic,         # model fit criterion - good for both large and small sample sizes
        level = 1,          # 1 without interactions, 2 with interactions
        method = "d",       # computes number of all possible models (doesn't execute)
        # other options include "h" (exhaustive) or "g" (genetic)
        family = gaussian,
        fitfunction = glm,  # specify the model type (lm or glm)
        confsetsize = 100)  # keep the 100 best models (confidence set)

# without interactions, the candidate set contains 32 models.
# with interactions, candidate set contains 32768 models!
# for aggregated SHT performance, we will exclude any interactions.

sht.r.model.g <- glmulti(as.formula(paste("Side_Hop_R ~", 
                                            # dynamically input variable names from subset 2
                                            paste(names.right, collapse = " + "))),
                           data = right.subset,
                           crit = bic,         # model fit criterion
                           level = 1,          # 1 without interactions, 2 with interactions
                           method = "g",       # genetic screening algorithm (better with fewer predictors)
                           family = gaussian,
                           fitfunction = glm,  # specify the model type (lm or glm)
                           confsetsize = 100)  # keep the 100 best models (confidence set)

# converged after 460 generations

# model results
print(sht.r.model.g)

# plot the best 100 models base on information criterion (IC)
plot(sht.r.model.g) # red line indicates 2 IC units from the best model

################################################################################
# Model-Averaged Importance Terms
################################################################################
# specify output parameters
png(filename = "Model Importance Terms (Right SHT).png",
    width = 2000, height = 1200, res = 150)

# plot margins
par(mar = c(5, 12, 1, 2) + 0.1)

# plot the model
plot(sht.r.model.g, type = "s")

# close the .png window
dev.off()

################################################################################
# summarise the best 3 models
################################################################################
weightable(sht.r.model.g)[1:3,] %>% # select the best models
   flextable::regulartable() %>%
   flextable::autofit()

# final model
m1 <- lm(Side_Hop_R ~
            # predictors
            Side_Hop_L + 
            HU_timePeakTorqueHeld_LFle_180 + 
            HU_timePeakTorque_LFle_60 + 
            HU_reciprocalDelayCV_RExt_60 + 
            HU_delayTimeCV_RExt_60 + 
            HU_jointAnglePeakTorqueCV_RExt_60 + 
            HU_workPRCV_LFle_180,
         data = d.sht.clean)


# create a table of the final model outputs, and save as a word .doc
library(sjPlot)
tab_model(m1, 
          show.std  = TRUE,
          show.df   = FALSE,
          show.se   = TRUE, 
          string.se = "SE",
          dv.labels = c("Final Model"),
          file = "Right SHT Models.doc")

# summarise the strongest model
summary(m1)
model_performance(m1)

################################################################################
# left model
################################################################################
vars.to.keep <- sht.l.VI$Variable[1:4]

left.subset <- d.sht.clean %>%
   select(Side_Hop_L, any_of(vars.to.keep))

glmulti(as.formula(paste("Side_Hop_L ~", 
                         # dynamically input top 10 variable names
                         paste(names.left, collapse = " + "))), 
        data = left.subset,
        crit = bic,         # model fit criterion - good for both large and small sample sizes
        level = 1,          # 1 without interactions, 2 with interactions
        method = "d",       # computes number of all possible models (doesn't execute)
                            # other options include "h" (exhaustive) or "g" (genetic)
        family = gaussian,
        fitfunction = glm,  # specify the model type (lm or glm)
        confsetsize = 100)  # keep the 100 best models (confidence set)

# without interactions, the candidate set contains 1024 models.
# with interactions, candidate set contains more than 1 billion (1e9) models!!!
# for aggregated SHT performance, we will exclude any interactions.

sht.l.model.g <- glmulti(as.formula(paste("Side_Hop_L ~", 
                                          # dynamically input variable names from subset 2
                                          paste(names.left, collapse = " + "))),
                         data = left.subset,
                         crit = bic,         # model fit criterion
                         level = 1,          # 1 without interactions, 2 with interactions
                         method = "g",       # genetic screening algorithm (better with fewer predictors)
                         family = gaussian,
                         fitfunction = glm,  # specify the model type (lm or glm)
                         confsetsize = 100)  # keep the 100 best models (confidence set)

# converged after 410 generations

# model results
print(sht.l.model.g)

# plot the best 100 models base on information criterion (IC)
plot(sht.l.model.g) # red line indicates 2 IC units from the best model

################################################################################
# Model-Averaged Importance Terms
################################################################################
# specify output parameters
png(filename = "Model Importance Terms (Left SHT).png",
    width = 2000, height = 1200, res = 150)

# plot margins
par(mar = c(5, 12, 1, 2) + 0.1)

# plot the model
plot(sht.l.model.g, type = "s")

# close the .png window
dev.off()

################################################################################
# summarise the best 3 models
################################################################################
weightable(sht.l.model.g)[1:3,] %>% # select the best models
   flextable::regulartable() %>%
   flextable::autofit()

# final model
m2 <- lm(Side_Hop_L ~
            # predictors
            Side_Hop_R + 
            HU_workPRCV_LFle_180 + 
            HU_timePeakTorque_LFle_60 + 
            HU_reciprocalDelayCV_LExt_180 + 
            HU_delayTime_LFle_60,
         data = d.sht.clean)


# create a table of the final model outputs, and save as a word .doc
tab_model(m2, 
          show.std  = TRUE,
          show.df   = FALSE,
          show.se   = TRUE, 
          string.se = "SE",
          dv.labels = c("Final Model"),
          file = "Left SHT Models.doc")

# summarise the strongest model
summary(m2)
model_performance(m2)
effect