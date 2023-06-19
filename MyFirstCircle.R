### Title: Baysian GLM with a circular distrobution
### Authors: Skylar Lynch and Brett Melbourne
### Date: 11/11/2022


#install.packages("brms")
#install.packages("remotes")
#install.packages("circglmbayes")
#install.packages("circular")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("rstanarm")


#Load Libraries
library(brms)
library(remotes)
library(circglmbayes)
library(circular)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstanarm) #nb overrides default ggplot theme
options(mc.cores=parallel::detectCores())

# hpdi.R wasn't working seperately so it is below

hpdi <- function (samp, prob = 0.95) {
  vals <- sort(samp)
  nsamp <- length(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- which.min(vals[init + gap,drop=FALSE] - vals[init, drop=FALSE])
  ans <- cbind(lower=vals[inds], upper=vals[inds + gap])
  return(ans)
}

# Set theme to black and white because that's the only color scheme I seem to be able to exist with in life

theme_set(theme_bw())

# Read in data

setwd("C:/Users/skly5321/Downloads/ChapterOne/ChapterOne/Data")

plant <- read.csv("Amazon.csv")


# Plot a blob of nonsence, because why not. Bassically shows flowering happens all year, and has consistently over time.

plant %>% 
  ggplot(mapping=aes(x=Year, y=Day.of.the.year)) +
  geom_point()

# Plot in a circle

plant %>% 
  mutate(doy_bin=cut_width(Day.of.the.year, width=365/12, boundary=0, labels=FALSE)) %>% 
  group_by(doy_bin) %>% 
  summarize(histo=n()) %>%
  mutate(doy_mid=doy_bin*365/12 - 365/24) %>% 
  ggplot() +
  geom_col(aes(x=doy_mid, y=histo)) +
  coord_polar() +
  scale_x_continuous(breaks=seq(365/24, 365, 365/12),
                     labels = c("J","F","M","A","M","J","J","A","S","O","N","D"),
                     limits = c(0,365)) +
  labs(x=NULL, y="Number of records")



plant %>% 
  mutate(doy_bin=cut_width(Day.of.the.year, width=365/12, boundary=0, labels=FALSE)) %>% 
  mutate(Period=cut_number(Year, n=6, dig.lab=4)) %>% 
  group_by(doy_bin, Period) %>% 
  summarize(histo=n()) %>%
  mutate(doy_mid=doy_bin*365/12 - 365/24) %>% 
  ggplot() +
  geom_col(aes(x=doy_mid, y=histo)) +
  facet_wrap(vars(Period)) +
  coord_polar() +
  scale_x_continuous(breaks=seq(365/24, 365, 365/12),
                     labels = c("J","F","M","A","M","J","J","A","S","O","N","D"),
                     limits = c(0,365)) +
  labs(main= "INPA Reserves", x=NULL, y="Number of records")




## Fit lm based on circular package

#Convert Julian dates to degrees

plant$degrees <- (plant$Day.of.the.year/366)*360

# Define x and y and make sure date is circular
x <- plant$Year
y <- circular(plant$degrees, units = c("degrees"))

# Fit a circular-circular linear regression
circ.lm <- lm.circular(y, x, order=1)

# Make data frame
data <- data.frame(x,y)

formula <- y ~ x + y:x



# Usage of circGLM

glm <- circGLM(
  formula=y ~ x + y:x,
  data=data.frame(y=y,x=x),
  X = if (missing(th)) {
    model.matrix(formula, data)[, -1, drop = FALSE]
  } else {
    
    matrix(nrow = length(th), ncol = 0)
  },
  conj_prior = rep(0, 3),
  bt_prior_musd = c(mu = 0, sd = 1),
  starting_values = c(0, 1, rep(0, 4)),
  bwb = rep(0.05, 4),
  Q = 1000,
  burnin = 5,
  thin = 1,
  kappaModeEstBandwith = 0.1,
  CIsize = 0.95,
  r = 2,
  returnPostSample = TRUE,
  returnLLEachPar = FALSE,
  output = "list",
  SDDBFDensEstMethod = "density",
  reparametrize = TRUE,
  groupMeanComparisons = TRUE,
  skipDichSplit = FALSE,
  centerOnly = FALSE
)


# Look at results
print(glm)
print(glm, type = "mcmc")
print(glm, type = "coef")
print(glm, type = "all")
summary(glm, digits=4)

# Plot
plot(glm)

par(mar=c(1, 1, 1, 1))

# Traceplot stack
plot(glm, type = "tracestack")

# Prediction plot
plot(glm, type = "predict")

# Mean comparisons
plot(glm, type = "meancompare")
plot(glm, type = "meanboxplot")


# hpdi values

hpdi(samples$alpha[,1], prob=0.95)
