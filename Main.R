#################################
#Install and load packages      #
#################################

#Packages to be used
packages<-c("readxl","here","tidyverse","ggplot2","fmsb","knitr","multcompView",
            "logistf","MASS", "car", "lmm","lme4","here", "haven","tidyr","dplyr",
            "lsmeans", "mgcv","optimx","flexsurv","survminer","lmerTest","splines", 
            "coxme", "gridExtra","BRugs","coda","rjags","rgl","survival","reghelper")



# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


#Creating directories
dir.create(here("Figures"),showWarnings = F)
dir.create(here("Outputs"),showWarnings = F)


#Halo diameter
source(here("Scripts","routine_halo.R"))

#MIC
source(here("Scripts","routine_MIC.R"))
