# GLMM workshop extras

library(glmmTMB)

# Based on examples in vignette here:
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf

# Data set Owls included in package
Owls <- data.frame(Owls)
str(Owls)

Owls <- transform(Owls,
                  Nest=reorder(Nest,NegPerChick),
                  NCalls=SiblingNegotiation,
                  FT=FoodTreatment)

# Zero-inflated model (ziformula just fits an intercept)
fit_zipoisson <- glmmTMB(NCalls~(FT+ArrivalTime)*SexParent+(1|Nest),
                         data=Owls,
                         ziformula=~1,
                         family=poisson)

summary(fit_zipoisson)

# Zero-inflation model estimates for probability response is 0
# Conditional model estimates for number of calls (when positive)


# Plotting with sjPlot

# install.packages("sjPlot")
library(sjPlot)

plot_model(sleep.lmer, type = "slope")  # model from other code
plot_model(sleep.lmer, type = "diag")  # see help(plot_model) for other options


