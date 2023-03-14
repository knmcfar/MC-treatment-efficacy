library(wgeesel)
library(mice)
library(VIM)
library(tidyverse)
library(geepack)
library(lme4)
library(broom.mixed)

#load in data and reshape-------------------------------------------------------
mcf <- read.csv("mucosal_candidiasis_final.csv")
View(mcf)

mcf_wide <- reshape(mcf, 
                    timevar = "months", 
                    idvar = "personid",
                    v.names = c("posculture"),
                    direction = "wide")
View(mcf_wide)

mcf_long <- reshape2::melt(mcf_wide, id = c(1:9))
mcf_long <- mcf_long[order(mcf_long$personid), ] #order the data
#now we have rows for each measurement, including missing values
View(mcf_long)

#exploratory analyses ----------------------------------------------------------

#investigate missingness
md.pattern(mcf_wide, plot = FALSE)

par(mfrow = c(2, 2))
spineMiss(mcf_wide[, c("trt", "posculture.3")])
spineMiss(mcf_wide[, c("trt", "posculture.6")])
spineMiss(mcf_wide[, c("trt", "posculture.9")])
spineMiss(mcf_wide[, c("trt", "posculture.12")])

aggr_plot <- aggr(mcf_wide[, 10:14],
                  numbers = TRUE,
                  sortVars = TRUE,
                  labels = names(mcf_wide[, 10:14]),
                  cex.axis = 0.7,
                  gap = 3,
                  ylab = c("Histogram of missing data","Pattern"))

#descriptive stats
mcf %>% 
  group_by(trt, months) %>% 
  summarise(n = n(), 
            prop.posculture = mean(posculture), 
            mean_cd4 = mean(cd4bl),
            mean_arusebl = mean(arusebl))

mcf_wide %>% 
  filter(mcf_wide$canobl == 1) %>% 
  summarise(n = n())

mcf_wide %>%
  filter(mcf_wide$canvbl == 1) %>% 
  summarise(n = n())

mcf_wide %>% 
  filter(mcf_wide$canobl == 1 & mcf_wide$canvbl == 1) %>%
  summarise(n = n())

#baseline characteristics
mcf %>% 
  group_by(trt) %>% 
  filter(months == 0) %>% 
  summarise(n = n(), 
            age = mean(age), 
            race = mean(white), 
            cd4 = mean(cd4bl), 
            prior = mean(podbl), 
            arv = mean(arusebl), 
            vc = mean(canvbl), 
            oc = mean(canobl), 
            presence = mean(posculture))

mcf %>% 
  group_by(trt) %>% 
  filter(months == 0) %>% 
  summarise(n = n(), age = sd(age), cd4 = sd(cd4bl))

#confirmatory analysis ---------------------------------------------------------

#long format imputation
mcf_long <- rename(mcf_long,
                   months = "variable",
                   posculture = "value")
mcf_long$months <- as.character(mcf_long$months)

mcf_long <- mcf_long %>% 
              mutate(months = replace(months, months == "posculture.0", 0))
mcf_long <- mcf_long %>% 
              mutate(months = replace(months, months == "posculture.3", 3))
mcf_long <- mcf_long %>% 
              mutate(months = replace(months, months == "posculture.6", 6))
mcf_long <- mcf_long %>%
              mutate(months = replace(months, months == "posculture.9", 9))
mcf_long <- mcf_long %>% 
              mutate(months = replace(months, months == "posculture.12", 12))

mcf_imp <- data.frame(mcf_long[c(1:11)])
View(mcf_imp)

pred <- make.predictorMatrix(mcf_imp)*2 
pred["posculture", "personid"] <- -2 # indicate clustering by id
imp_long <- mice(mcf_imp, 
                 method = "2l.bin", 
                 predictorMatrix = pred, 
                 seed = 20, 
                 maxit = 1, 
                 m = 20, 
                 print = FALSE)

densityplot(imp_long, xlab = "Imputations of Candidiasis Fungus")

fit_imp_longLMM <- with(imp_long, 
                        glmer(posculture ~ trt*as.numeric(months) + (1+as.numeric(months)|personid), 
                              family = binomial(link = "logit"), 
                              nAGQ = 1))
summary(pool(fit_imp_longLMM))

#model w/ no imputations
glmm_fit <- glmer(posculture ~ trt*as.numeric(months) + (1+as.numeric(months)|personid), 
                  data = mcf, 
                  family = binomial(link="logit"), 
                  nAGQ = 1)
summary(glmm_fit)

#diagnostics: best and worse case models----------------------------------------

#best case: treatment A
mcfl_worst <- mcf_long
mcfl_worst[mcfl_worst$trt == 1 & is.na(mcfl_worst$posculture), "posculture"] <- 1
mcfl_worst[mcfl_worst$trt == 0 & is.na(mcfl_worst$posculture), "posculture"] <- 0

mcfl_worst %>% 
  group_by(trt, months) %>% 
  summarise(n = n(), proportion = mean(posculture))

glmm_fit2 <- glmer(posculture ~ trt*as.numeric(months) + (1+as.numeric(months)|personid), 
                   data = mcfl_worst, 
                   family = binomial(link = "logit"), 
                   nAGQ = 1)
summary(glmm_fit2)

#best case: treatment B
mcfl_best <- mcf_long
mcfl_best[mcfl_best$trt == 1 & is.na(mcfl_best$posculture), "posculture"] <- 0
mcfl_best[mcfl_best$trt == 0 & is.na(mcfl_best$posculture), "posculture"] <- 1

mcfl_best %>% 
  group_by(trt, months) %>% 
  summarise(n = n(), proportion = mean(posculture))

glmm_fit3 <- glmer(posculture ~ trt*as.numeric(months) + (1+as.numeric(months)|personid), 
                   data = mcfl_best, 
                   family = binomial(link = "logit"), 
                   nAGQ = 1, 
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
summary(glmm_fit3)
