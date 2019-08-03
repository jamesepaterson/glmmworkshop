# Introduction to Generalized Linear Mixed-Effects Models
# 2019/02/06
# James Paterson


# Outline:
# 1. Fitting GLMMs
# 2. Diagnostics
# 3. Inference & interpretation
# 4. Predictions
# 5. Graphing
# 6. Non-gaussian example


##### Fitting GLMMs with lme4 -------------------------------------------------

library(lme4)

# Use example
sleepstudy <- data.frame(sleepstudy)
str(sleepstudy)

# Random slopes & intercepts model
# Using "lmer" implies family = "gaussian"
sleep.lmer <- lmer(Reaction~Days + (Days|Subject), 
                   data = sleepstudy) # data

# Random intercepts model
# ________
sleep.lmer2 <- lmer(Reaction~Days + (1|Subject), 
                   data = sleepstudy) # data




##### Diagnostics -----------------------------------------------------------

# Plots
par(mfrow = c(2, 2))
hist(resid(sleep.lmer),breaks = 50)
qqnorm(resid(sleep.lmer))
qqline(resid(sleep.lmer))
plot(resid(sleep.lmer)~fitted(sleep.lmer))
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(sleep.lmer), residuals(sleep.lmer)))
hist(sleepstudy$Reaction, breaks = 50)

# Leverage plot. Maybe useful, but doesn't show Cook's Distance
library(ggplot2)
ggplot(data.frame(lev=hatvalues(sleep.lmer),pearson=residuals(sleep.lmer,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw()


##### Inference ---------------------------------------------------------------

# Output of model with main effect estimates
summary(sleep.lmer)

# Random effect estimates
ranef(sleep.lmer)

# Quick plot of random effects
# install.packages("lattice")
lattice::dotplot(ranef(sleep.lmer, condVar = TRUE))

# Testing the main effects
# install.packages("car")
library(car)
Anova(sleep.lmer, test = "F")

# Testing the effect of a random effect
# Drop random effect and compare two models with a LRT
sleep.lmer.null <- lmer(Reaction~Days + (1|Subject), sleepstudy)
anova(sleep.lmer.null, sleep.lmer)

##### Predictions -------------------------------------------------------------

# Predict without random effects (on Population Level)
new.sleep.df <- data.frame(Days = 0:9)
# re.form = NA sets random effects to zero 
new.sleep.df$pred <- predict(sleep.lmer, re.form = NA, 
                                 type = "response", newdata = new.sleep.df)

# Confidence intervals
# From: https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html#prediction
easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link function
  linkinv <- family(model)$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low=pred0-crit*pred.se,
                conf.high=pred0+crit*pred.se))
}

sleep.pred.CI <- easyPredCI(sleep.lmer, new.sleep.df)

# Put in same dataframe
new.sleep.df <- cbind(new.sleep.df, sleep.pred.CI)

# Examine dataframe
new.sleep.df

##### Graphing ----------------------------------------------------------------

library(ggplot2)

# Plot raw data + prediction (with confidence intervals)
ggplot(data = new.sleep.df, aes(x = Days, y = pred)) +
  geom_point(data = sleepstudy, aes(x = Days, y = Reaction)) +  # raw data points
  stat_smooth(col = "black", method = "lm", se = FALSE) +  # Predicted regression
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +  # Confidence interval
  labs(x = "Days", y = "Reaction (s)") 
 stat_smooth(data = sleepstudy, aes(x = Days, y = Reaction), col = "red", method = "lm")

# Plotting individual level effects
sleepstudy2 <- sleepstudy
sleepstudy2$pred <- predict(sleep.lmer, sleepstudy) # predict based on fixed and random effects

ggplot(data = sleepstudy2, aes(x = Days, y = pred, col = Subject, group = Subject)) +
  stat_smooth(method = "lm") +  # Predicted regression for each individual
  geom_point(data = sleepstudy, aes(x = Days, y = Reaction, col = Subject)) +
  labs(x = "Days", y = "Reaction (s)")  # Axis titles

##### Non-gaussian example ----------------------------------------------------

# Example dataset in richdf.RData (object = rich.df)
load(file = "richdf.RData")

# Examine structure
str(rich.df)

# Basic plot
plot(richness~dist.disturb, rich.df, col = as.factor(site))

# Try linear model
rich.lmer <- lmer(richness~dist.disturb + (1|site), data = rich.df)
# singular fit

# Diagnostics
par(mfrow = c(2, 2))
hist(resid(rich.lmer),breaks = 50)
qqnorm(resid(rich.lmer))
qqline(resid(rich.lmer))
plot(resid(rich.lmer)~fitted(rich.lmer))
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(rich.lmer), residuals(rich.lmer)))
hist(rich.df$richness, breaks = 50)
# Model doesn't fit data well

rich.glmer.p <- glmer(richness~dist.disturb + (1|site), 
                      data = rich.df,
                      family = "poisson") # mode convergence / scale warning

rich.glmer.p.opt2 <- glmer(richness~dist.disturb + (1|site), 
                           data = rich.df,
                           family = "poisson",
                           control=glmerControl(optimizer="bobyqa")) # mode convergence / scale warning

# Scale data first
rich.df$scaled.dist.disturb <- scale(rich.df$dist.disturb)
rich.glmer.p.scale <- glmer(richness~scaled.dist.disturb + (1|site), 
                            data = rich.df,
                            family = "poisson")

# Diagnostics on scale error
par(mfrow = c(2, 2))
hist(resid(rich.glmer.p.scale),breaks = 50)
qqnorm(resid(rich.glmer.p.scale))
qqline(resid(rich.glmer.p.scale))
plot(resid(rich.glmer.p.scale)~fitted(rich.glmer.p.scale))
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(rich.glmer.p.scale), residuals(rich.glmer.p.scale)))
hist(rich.df$richness, breaks = 50)

# Summary
summary(rich.glmer.p.scale)

# Testing main effects with LRT (drop main effect, compare with anova)
rich.glmer.p.scale.null <- glmer(richness~(1|site), 
                            data = rich.df,
                            family = "poisson")
anova(rich.glmer.p.scale, rich.glmer.p.scale.null)

# Predict on new data
rich.new.data <- data.frame(scaled.dist.disturb = scale(0:500))
rich.new.data$pred <- predict(rich.glmer.p.scale, rich.new.data, re.form = NA, type = "response")
rich.new.data <- cbind(rich.new.data,
                       data.frame(easyPredCI(rich.glmer.p.scale, rich.new.data)))


# ggplot
ggplot(data = rich.new.data, aes(x = scaled.dist.disturb, y = pred)) +
  geom_point(data = rich.df, aes(x = scaled.dist.disturb, y = richness)) +
  stat_smooth(method = "loess", se = FALSE) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  labs(x = "Scaled distance to disturbance", y = "Species richness")




# sjPlot for graphs ------------------------------------------------------------------


# install.packages("sjPlot")
library(sjPlot)

plot_model(sleep.lmer)

plot_model(sleep.lmer, type = "diag")

plot_model(sleep.lmer, type = "slope")

plot_model(rich.glmer.p.scale, type = "diag")

plot_model(rich.glmer.p.scale, type = "pred")

### Zero-inflated

library(glmmTMB)

zi.model <- glmmTMB(response~fixed + (1|random),
                    ziformula =~fixed + (1|random),
                    family = "nbinom1",
                    data = mydata)
