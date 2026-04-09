
# this script runs statistical models for 3 side hop test (SHT) studies.
# 
# study 1: how do anthropometric constraints contribute/modulate SHT performance in collegiate athletes.
# study 2: how do DI athletes achieve faster SHT times?
# study 3: how does SHT strategy differ between CAI and non-CAI athletes?

library(tidyverse)
library(DHARMa)
library(readxl)

# -------------------------
# load data
# -------------------------
data.path <- "SHT Data_biomech.xlsx"

data.study1 <- read_excel(data.path, sheet = "Study 1") %>%
   mutate(`Study no.` = 1)
data.study2 <- read_excel(data.path, sheet = "Study 2") %>%
   mutate(`Study no.` = 2)
data.study3 <- read_excel(data.path, sheet = "Study 3") %>%
   mutate(`Study no.` = 3)

# -------------------------
# study 1 model(s)
# -------------------------
library(sjPlot)

# "Sport" did not add meaningful explanatory power to m0, so it was removed.
# It will subsequently be excluded from all further models

# -------------------------
# core anthropometrics
# -------------------------
m0 <- data.study1 %>%
   select(-R_Ankle_DF, -L_Ankle_DF) %>%       # remove ankle variables
   lm(SHT_agg ~ Height_m + Mass_kg, data = .) # use reduced dataset

# interaction
m0a <- data.study1 %>%
   select(-R_Ankle_DF, -L_Ankle_DF) %>%       # remove ankle variables
   lm(SHT_agg ~ Height_m + Mass_kg, data = .) # use reduced dataset

tab_model(m0, 
          show.df   = TRUE,
          show.se   = TRUE, string.se = "SE",
          dv.labels = "",
          title     = "Table 1. Side hop test (SHT) performance predicted by height and mass.",
          file      = "Tables & Figures/Table 1.doc")

# the core anthropometrics model suggests:
#     - SHT time increase with height (0.58s with every 10cm)
#     - SHT time decrease with body mass (0.04s with every +1kg)
#     - these two variables only explain ~21% of the variance in SHT time though.
# 
# taller athletes tend to be slower.
#     - taller athletes have longer levers and potentially greater lateral displacement
# heavier athletes (for a given height) tend to be faster.
#     - mass (likely reflecting strength) aids force production.

# ------------------------------------
# substituting leg length for height
# ------------------------------------
m1 <- data.study1 %>%
   select(-R_Ankle_DF, -L_Ankle_DF) %>%   # remove ankle variables
   lm(SHT_agg ~ Leg_m_agg + Mass_kg, 
      data = .)                           # use reduced dataset

tab_model(m0, m1, 
          show.df   = TRUE,
          show.se   = TRUE, string.se = "SE",
          dv.labels = c("Height", "Leg Length"),
          title     = "Table 2. Side hop test (SHT) performance model comparison: height and mass vs leg length and mass.",
          file      = "Tables & Figures/Table 2.doc")

# ------------------------------------
# subset anthropometrics and ankle DF
# ------------------------------------
# core anthropometrics
m2_base <- data.study1 %>%
   filter(Sport != "MBB") %>%     # remove MBB
            lm(SHT_agg ~ Height_m + Mass_kg,
               data = .)

# including ankle DF ROM
m2_full <- data.study1 %>%
   filter(Sport != "MBB") %>%     # remove MBB
   lm(SHT_agg ~ Height_m + Mass_kg +
         R_Ankle_DF + L_Ankle_DF,
      data = .)

tab_model(m2_base, m2_full, 
          show.df   = TRUE,
          show.se   = TRUE, string.se = "SE",
          dv.labels = c("Reduced", "Full"),
          title     = "Table 3. Model comparisons for predicting side hop test (SHT) performance in a subset sample with the inclusion of ankle dorsiflexion range of motion (ROM).",
          file      = "Tables & Figures/Table 3.doc")

# there is no clear evidence that ankle DF ROM meaningfully improved prediction of 
# SHT time in this subset sample of DI collegiate athletes.


# substituting mean leg length for height yielded similar results, indicating that 
# overall body size, rather than specific segment length, underpins the association.

# summary of findings...
#     - taller athlete (& longer leg length) ≈ slower SHT
#     - greater mass ≈ slightly faster SHT
#     - ankle DF ROM ≈ minimal additional explanatory value
#     - anthropometrics explain ~15–20% of variance (adj. R²)

m3 <- data.study1 %>%
   filter(Sport != "WVB") %>%
   lm(SHT_agg ~ CAI + Height_m + Mass_kg,
      data = .)

m4 <- data.study1 %>%
   filter(Sport != "WVB") %>%
   lm(SHT_agg ~ CAI,
      data = .)

tab_model(m4, m3, 
          show.df   = TRUE,
          show.se   = TRUE, string.se = "SE",
          dv.labels = c("Reduced", "Full"),
          title     = "Table 4. Side hop test (SHT) performance CAI.",
          file      = "Tables & Figures/Table 4.doc")

# Although CAI participants demonstrated a mean difference of approximately 0.4s 
# slower SHT performance, the confidence interval was wide and included trivial 
# effects. Therefore, the present data do not provide clear evidence that SHT 
# performance differs by CAI status when accounting for anthropometric 
# characteristics.
