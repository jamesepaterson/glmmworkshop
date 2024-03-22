## ---------------------------
##
## Script name: 02_glmm_workshop_code.R 
##
## Purpose of script: Code for introduction to Generalized Linear Mixed-Effects Models
##
## Author: James E Paterson
##
## Date Updated: 2024-03-21
##
## 
##
## ---------------------------

# A workflow for fitting GLMMs
# 1. Fit model(s)
# 2. Diagnostics and assumptions tests with base R or DHARMa
# 3. Make inferences (test effects and effect sizes)
# 4. Make predictions
# 5. Plot
# ...Iterate, as needed!

# Outline:
# 1. Gaussian Example Workflow
# 2. Count data Example Workflow
# 3. Zero-inflated Example Workflow

# Load libraries
library(dplyr)
library(lme4)
library(glmmTMB)
library(performance)
library(car)
library(DHARMa)
library(ggplot2)

##### 1. Gaussian Example Workflow ----------------------------------------------------


## Fit model(s) -------------------------------------------------

# Use example
sleepstudy <- data.frame(sleepstudy)
str(sleepstudy)

# Random slopes & intercepts model
# Using "lmer" implies family = "gaussian"
sleep.lmer <- lmer(Reaction~Days + (Days|Subject), 
                   data = sleepstudy)

# Random intercepts model
sleep.lmer2 <- lmer(Reaction~Days + (1|Subject), 
                   data = sleepstudy)

##### Diagnostics -----------------------------------------------------------

# For diagnostics, examine normality of residuals (left) and patterns in residuals (right)

# Diagnostic plots with base R
par(mfrow = c(1, 2))
qqnorm(resid(sleep.lmer))
qqline(resid(sleep.lmer))
plot(resid(sleep.lmer)~fitted(sleep.lmer))
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(sleep.lmer), residuals(sleep.lmer)))

# Diagnostic plots with DHARMa
# See the vignette for more details: 
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
simulationOutput <- simulateResiduals(fittedModel = sleep.lmer, 
                                      plot = T)

##### Inference ---------------------------------------------------------------

# Output of model with main effect estimates
summary(sleep.lmer)

# Random effect estimates
ranef(sleep.lmer)

# Quick plot of random effects
# install.packages("lattice")
lattice::dotplot(ranef(sleep.lmer, condVar = TRUE))

# Testing the main effects
# type = 2 when no interactions, type 3 when there are interactions
Anova(sleep.lmer, test = "F", type = 2)

# Testing the effect of a random effect
# Drop random effect and compare two models with a LRT
# Requires models with with ML instead of REML
anova(sleep.lmer, sleep.lmer2)
# Strong evidence that there are different slopes (Reaction by Days) for each subject

##### Predictions -------------------------------------------------------------

# Function for Confidence intervals
# From: https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html#prediction
# only for linear mixed-effects models
easyPredCI <- function(model,newdata = NULL, alpha = 0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form = NA,newdata = newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model, fixed.only = TRUE)[-2],newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link function
  linkinv <- family(model)$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low = pred0-crit*pred.se,
                conf.high = pred0+crit*pred.se))
}

# Predict without random effects (on Population Level)
new.sleep.df <- data.frame(Days = 0:9) %>%
  # re.form = NA sets random effects to zero (the population mean)
  mutate(pred = predict(sleep.lmer, re.form = NA, 
                        type = "response", 
                        newdata = .)) %>%
  cbind(., 
        easyPredCI(sleep.lmer, .))

# Examine dataframe
new.sleep.df

# Predictions for random effects (adding to main dataframe)
sleepstudy <- sleepstudy %>%
  mutate(pred_individuals = predict(sleep.lmer, sleepstudy)) # predict based on fixed and random effects

##### Plotting ----------------------------------------------------------------

# Presenting two ways to plot data: showing population effect (no random effects) and showing random effects
# Choose approach that best suits your questions or interest and study objectives.

# Plot raw data + population prediction (with confidence intervals)
ggplot() +
  # raw data points
  geom_point(data = sleepstudy, aes(x = Days, y = Reaction), pch = 21) +  
  # Predicted population mean values
  stat_smooth(data = new.sleep.df, aes(x = Days, y = pred),
              formula = 'y~x', col = "black", method = "lm", se = FALSE) +
  # Confidence interval of predicted population mean values
  geom_ribbon(data = new.sleep.df,
              aes(x = Days, ymin = conf.low, ymax = conf.high), alpha = 0.3) +  
  labs(x = "Days", y = "Reaction (s)") +
  theme_classic()

# Plot raw data + individual (random) effects
ggplot() +
  # raw data points
  geom_point(data = sleepstudy, aes(x = Days, y = Reaction, col = Subject), pch = 21) +  
  # Predicted regression for each individual
  stat_smooth(data = sleepstudy, aes(x = Days, y = pred_individuals, col = Subject),
              method = "lm", formula = 'y~x') +  
  labs(x = "Days", y = "Reaction (s)") +
  theme_classic()

##### Example 2 (Count data) ----------------------------------------------------

# Same workflow:
# Fit
# Diagnostics
# Inference
# Predict 
# Plot

# This example deliberately takes iterations because that's often how analyses happen

# Example dataset in richdf.RData (object = rich.df)
load(file = "richdf.RData")

# Examine structure
str(rich.df)

# Basic plot
ggplot(data = rich.df, aes(x = dist.disturb, y = richness, col = as.factor(site))) +
  geom_point(pch = 21) +
  theme_classic()
         
# Try Gaussian model
rich.lmer <- lmer(richness~dist.disturb + (1|site), data = rich.df)
# singular fit = standard deviation of site ~= 0

# Diagnostics (DHARMa)
simulationOutput <- simulateResiduals(fittedModel = rich.lmer, 
                                      plot = T)
# Model doesn't fit data well

# Fit with a Poisson model (a good starting point for count data)
rich.glmer.poisson <- glmer(richness~dist.disturb + (1|site), 
                      data = rich.df,
                      family = "poisson") 
# model convergence/ scale warning

# A few potential solutions: rescale predictor, change optimizer, increase the maximum number of iterations to converge
rich.glmer.poisson.opt2 <- glmer(richness~dist.disturb + (1|site), 
                           data = rich.df,
                           family = "poisson",
                           control=glmerControl(optimizer = "bobyqa",
                                                optCtrl = list(maxfun = 10^8)))
# model convergence/ scale warning

# That didn't work so we'll scale the predictor variable (mean = 0, sd = 1). 
# See help(scale) for details
rich.df <- rich.df %>%
  mutate(scaled.dist.disturb = scale(dist.disturb))

# Fit Poisson model on scaled predictor data
rich.glmer.poisson.scale <- glmer(richness~scaled.dist.disturb + (1|site), 
                            data = rich.df,
                            family = "poisson")

# Diagnostics (DHARMa) on scaled data model
simulationOutput <- simulateResiduals(fittedModel = rich.glmer.poisson.scale, 
                                      plot = T)
# Model fit still has issues in the residual patterns

# Is there overdispersion? (Poisson assumes variance = mean, overdispersed when the variance >> mean)
check_overdispersion(rich.glmer.poisson.scale)
# Overdispersion detected

# We're running out of options in lme4. 
# To deal with overdispersion, we can use a negative binomial family (variance > mean, by estimating dispersion parameter), 
# but need to switch to glmmTMB package
rich.negbin <- glmmTMB(richness~scaled.dist.disturb + (1|site), 
                       data = rich.df,
                       family = nbinom1) # There are two formulations of the nbinom distribution available

# Diagnostics (DHARMa) on scaled data model
simulationOutput <- simulateResiduals(fittedModel = rich.negbin, 
                                      plot = T)
# Fit is much improved

# Is there overdispersion? (For count data, is the variance >> mean)
check_overdispersion(rich.negbin)
# No overdispersion detected (i.e. negative binomial family worked to deal with that issue)

# Summary
summary(rich.negbin)

# Testing main effects with LRT (drop main effect, compare with anova)
rich.negbin.null <- glmmTMB(richness~(1|site), 
                            data = rich.df,
                            family = nbinom1)
anova(rich.negbin, rich.negbin.null)
# Strong evidence richness is related to scaled.dist.disturb

# Predict on new data
rich.new.data <- data.frame(scaled.dist.disturb = scale(0:500),
                            site = 1) %>%
  mutate(pred = predict(rich.negbin, ., 
                        re.form = NA, type = "response"),
         pred_se = predict(rich.negbin, ., 
                           re.form = NA, type = "response", se.fit = TRUE)$se.fit,
         lower_conf = pred-1.96*pred_se,
         upper_conf = pred+1.96*pred_se)

# ggplot
ggplot(data = rich.new.data, aes(x = scaled.dist.disturb, y = pred)) +
  geom_point(data = rich.df, aes(x = scaled.dist.disturb, y = richness), pch = 21) +
  stat_smooth(method = "loess", formula = 'y~x', se = FALSE, col = "black") + 
  geom_ribbon(aes(ymin = lower_conf, ymax = upper_conf), alpha = 0.3) +
  labs(x = "Scaled distance to disturbance", y = "Species richness") +
  theme_classic()

##### More complicated example: zero-inflation --------------------------------

# Use Salamander data set, 
# see more in glmmTMB vignette: https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf
# I include zero-inflated models below, but glmmTMB also can fit Hurdle models
data("Salamanders")

# Fit zero-inflated poisson or negative binomial
salamanders_nb1 <- glmmTMB(count ~ spp + (1|site),
                    family = nbinom1(), 
                    zi = ~1, # zero inflation formula for probability of a structural zero (logit link)
                    data = Salamanders)

salamanders_nb1_alt <- glmmTMB(count ~ spp + (1|site),
                           family = nbinom1(), 
                           zi = ~spp, # zero inflation formula for probability of a structural zero (logit link)
                           data = Salamanders)

# Diagnostics (DHARMa) on scaled data model
simulationOutput <- simulateResiduals(fittedModel = salamanders_nb1, 
                                      plot = T)

simulationOutput <- simulateResiduals(fittedModel = salamanders_nb1_alt, 
                                      plot = T)


# Inference
AIC(salamanders_nb1, salamanders_nb1_alt)

# Model summary
summary(salamanders_nb1)

# Predict


# Plot


