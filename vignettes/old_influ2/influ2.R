## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE, message=FALSE-------------------------------------------------
library(knitr)
library(tidyverse)
library(reshape2)
library(brms)
library(influ2)
library(bayesplot)

theme_set(theme_bw())

options(mc.cores = 2)

## ----sim-data, echo=TRUE, message=FALSE---------------------------------------
data(lobsters_per_pot)
nrow(lobsters_per_pot)
head(lobsters_per_pot)

## ----plot-bubble-1, echo=TRUE, fig.height=6, fig.width=6, message=FALSE, fig.cap = "Distribution of data by year and month."----
plot_bubble(df = lobsters_per_pot, group = c("year", "month")) +
  labs(x = "Month", y = "Year")

## ----plot-bubble-2, echo=TRUE, fig.height=6, fig.width=6, message=FALSE, fig.cap = "Distribution of data by year and month, also coloured by month."----
plot_bubble(df = lobsters_per_pot, group = c("year", "month"), 
            fill = "month") +
  labs(x = "Month", y = "Year") +
  theme(legend.position = "none")

## ----fit-model----------------------------------------------------------------
fit1 <- brm(lobsters ~ year + month + poly(soak, 3), 
            data = lobsters_per_pot, family = poisson, 
            chains = 2, iter = 1500, refresh = 0, seed = 42, 
            file = "fit1", file_refit = "never")

fit2 <- brm(lobsters ~ year + (1 | month) + s(depth, k = 3), 
            data = lobsters_per_pot, family = negbinomial(), 
            control = list(adapt_delta = 0.99),
            chains = 2, iter = 3000, refresh = 0, seed = 1, 
            file = "fit2", file_refit = "never")

## ----add-criterion------------------------------------------------------------
criterion <- c("loo", "waic", "bayes_R2")
fit1 <- add_criterion(x = fit1, criterion = criterion, file = "fit1")
fit2 <- add_criterion(x = fit2, criterion = criterion, file = "fit2")

## ----do-glm-------------------------------------------------------------------
fit_glm <- glm(lobsters ~ year + month + poly(soak, 3), 
               data = lobsters_per_pot, family = poisson)

## ----influ1-setup, echo=TRUE, fig.height=6, fig.width=6, message=FALSE, fig.cap = "The original CDI plot from the influ package."----
myInfl <- Influence$new(model = fit_glm)
myInfl$init()
myInfl$calc()
myInfl$cdiPlot(term = "month")

## ----plot-cdi-fe, echo=TRUE, fig.height=6, fig.width=6, message=FALSE, fig.cap = "The new Bayesian CDI plot from the influ2 package."----
cdi_month <- plot_bayesian_cdi(fit = fit1, xfocus = "month", yfocus = "year", 
                               xlab = "Month", ylab = "Year")
cdi_month

## ----influ1-step, echo=TRUE, fig.height=6, fig.width=6, message=FALSE, fig.cap = "The original step-plot from the influ package."----
myInfl$stepPlot()

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE, fig.cap = "The new Bayesian step-plot from the influ2 package. Note that the new step-plot requires that all models/steps be run in brms before the function can be used."----
fits <- list(fit1, fit2)
plot_step(fits = fits)

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE----------------------
plot_index(fit = fit2)

## ----influ1-plot, echo=TRUE, fig.height=6, fig.width=6, message=FALSE---------
myInfl$stanPlot()

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE----------------------
plot_influ(fit = fit1)

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE----------------------
myInfl$influPlot()

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE, fig.cap = "Comparison of influences produced by the original influ package (points) and the new influ2 calculation (lines)."----
i1 <- myInfl$influences %>%
  pivot_longer(cols = -level, names_to = "variable") %>%
  rename(year = level)
i_month <- get_influ(fit = fit1, group = c("year", "month")) %>%
  group_by(year) %>%
  summarise(value = mean(delta)) %>%
  mutate(variable = "month")
# i_soak <- get_influ(fit = fit1, group = c("year", "soak")) %>%
#   group_by(year) %>%
#   summarise(value = mean(delta)) %>%
#   mutate(variable = "poly(soak, 3)")
# i2 <- bind_rows(i_month, i_soak)
# 
# ggplot(data = i2, aes(x = year, y = exp(value))) +
#   geom_point(data = i1) +
#   geom_line(aes(colour = variable, group = variable)) +
#   facet_wrap(variable ~ ., ncol = 1) +
#   labs(x = "Year", y = "Influence") +
#   theme(legend.position = "none")

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE----------------------
fit1$criteria$loo
fit1$criteria$waic

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE----------------------
loo_compare(fit1, fit2, criterion = "loo") %>% kable(digits = 1)

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE----------------------
loo_compare(fit1, fit2, criterion = "waic") %>% kable(digits = 1)

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE----------------------
# get_bayes_R2(fits)

## ----echo=TRUE, fig.height=6, fig.width=6, message=FALSE, fig.cap = "CDI plot for a fixed effect."----
# cdi_month

