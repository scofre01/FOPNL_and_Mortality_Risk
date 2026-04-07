# Author: Sebastian Cofre
# Date: 2025-12-03
# Description: Data processing and analysis procedure to assess the causal effects of FOPNL on all-cause mortality in the 2016-2017 Chilean National Health Survey.

#libraries
library(dplyr)
library(mice)
library(naniar)
library(AIPW)
library(tmle)
library(SuperLearner)
library(ggplot2)

#Load complete database 
load("./input_data/bbdd.RData")

#=================================================
#============== Data pre-processing ==============
#================================================

#### We select the variables of interest for the analysis and create the outcome variable (death) and the time variable (time to death or censoring) and the intervention variable (warning labels based on the consideration of warning labels).
df <-bbdd |> 
  select(Event,FOPLN,Sex,Age,Education,DM,BMI,waist,pa,ethnicity,alcohol,HTN)

### We exclude all people with FOPNL = 3 (consideration of warning labels) to compare only those with FOPNL = 1 (no consideration of warning labels) and FOPNL = 2 (consideration of warning labels)
df <- df %>%
  filter(FOPLN != "3") ### filtramos todas personas con sellos 3 


#=================================================
#============== Estimation of causal effects (ATE) ====
#================================================

### We created an object that include the covariates 
cov = c("Sex","Age","Education","DM","BMI","waist","pa","ethnicity","alcohol","HTN")

##1. Estimation based on AIPWL
AIPW_SL <- AIPW$new(Y= df$Event,
                    A= df$FOPNL,
                    W= subset(df,select=cov), 
                    Q.SL.library = c("SL.mean","SL.glm"),
                    g.SL.library = c("SL.mean","SL.glm"),
                    k_split = 10,
                    verbose=FALSE)$
  fit()$
  #Default truncation
  summary(g.bound = c(0.025,0.975))$
  plot.p_score()$
  plot.ip_weights()

AIPW_SL$summary(g.bound = 0.025) 
AIPW_SL$result ### results
suppressWarnings({
  AIPW_SL$stratified_fit()$summary()
}) # stratify by treatment
AIPW_SL$g.plot##

##2. Estimation based on TMLE

tmle_fit <- tmle(Y=df$Event,
                 A=df$FOPNL,
                 W=subset(df,select=cov), 
                 Q.SL.library=c("SL.mean","SL.glm"),
                 g.SL.library=c("SL.mean","SL.glm"),
                 family="binomial",
                 cvQinit = TRUE)

cat("\nEstimates from TMLE\n")
unlist(tmle_fit$estimates$ATE)
unlist(tmle_fit$estimates$RR)
unlist(tmle_fit$estimates$OR)

cat("\nEstimates from TMLE\n")
a_tmle <- AIPW_tmle$
  new(A=df$FOPNL,Y=df$Event,tmle_fit = tmle_fit,verbose = TRUE)$
  summary(g.bound=0.025)

a_tmle$result ### see results


##3. Estimation based on IPTW
# Calculate  PS with WeightIt
library(WeightIt)
W.out <- weightit(FOPNL ~ Sex + Age +Education + DM+ Waist + pa + ethnicity+ alcohol+ HTN,
                  data = df,
                  method = "ps", # means propensity score
                  estimand = "ATE",
                  stabilize = TRUE) # Average Treatment Effect

# Add weights to the original data frame 
bd$iptw_weights <- W.out$weights

# We estimated the effect of the intervention (FOPNL) on the outcome (death) using a weighted logistic regression model.
model_iptw <- glm(Event ~ FOPLN,
                   data = df,
                   family = quasibinomial(),
                   weights = iptw_weights)

summary(model_iptw)

#  Extract the Odds Ratio
or_iptw <- exp(coef(model_iptw)["FOPLNSi"])
ci_iptw <- exp(confint(model_iptw)["FOPLNSi", ]) #with confidence intervals

###Results
cat("--- Results of Analysis based on IPTW ---\n")
cat(sprintf("Odds Ratio (OR) estimated: %.3f\n", or_iptw))
cat(sprintf("95% CI: [%.3f, %.3f]\n\n", ci_iptw[1], ci_iptw[2]))

###################################
# Save the results in an RData file 
save(df1, file = "df1.RData")

