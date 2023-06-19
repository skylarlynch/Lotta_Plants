

library(circglmbayes)
library(coda) #for mcmc tools to use with circglmbayes
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)

par(mar = c(1, 1, 1, 1))

setwd("C:/Users/skly5321/Downloads/ChapterOne/ChapterOne/Data")


# The link function for a circular GLM is the tan half function:
#   tan_half(x) = tan(x/2)
# So the inverse link function is the inverse tan half

inv_tan_half <- function(x) {
  return(2 * atan(x))
}

# Converting day to circular

day2circ <- function(day_of_year) {
  return(2 * pi * day_of_year / 365 - pi)
}

circ2day <- function(circ) {
  return(365 * (circ + pi) / (2 * pi))
}

day2circ-leap <- function(day_of_year) {
  return(2 * pi * day_of_year / 366 - pi)
}

#For a circular response y on -pi:pi and a single
# continuous predictor x, the model is:

# $$
# y = \beta_0 + inv_tan_half(\beta_1 * x)
# $$ 




v <- c("Amazon", "Bangladesh", "Bolivia", "Caatinga", "Cameroon", "CochaCashu(new)", "FrenchPolynesia", "Ghana", "Guinnea", "JatunSacha", "LasCruces", "MadhyaPradesh", "Miombo", "Myanmar", "Udzungwa", "WesternGhatsIndia")

traceplots_list <- list()
finalplot_list <- list()
coef_list <- list()

for (i in 1:16){
  df <- read.csv(paste(v[i],".csv", sep = ""))
  df$species <- factor(df$species)
  spec <- levels(df$species)
  for(j in length(spec)){
    plant2 <- na.omit(df[df$species == spec[j],])
    year <- plant2$Year
    year_s <- scale(year) #year scaled and centered
    year_center <- attr(year_s, "scaled:center")
    year_scale <- attr(year_s, "scaled:scale")
   
    
     for (m in 1:year)
      lp <- leap_year(year)
      if(lp == FALSE){
        circ <- day2circ(plant2$Day.of.the.year)
    } else {
        circ <- day2circ-leap(plant2$Day.of.the.year)
    }
    
   
     #circ <- day2circ(plant2$Day.of.the.year)
    sdat <- data.frame(year, year_s, circ)
    nchains <- 4
    chains <- list()
    for (k in 1:nchains ) {
      fit <- circGLM(circ ~ year_s, data=sdat)
      chains[[k]] <- fit$all_chains
    }
    chains <- mcmc.list(chains)
    traceplots_list[[spec[j]]] <- traceplot(chains)
    autocorr.plot(chains, ask=FALSE)
    gelman.diag(chains[,"b0_chain"])     
    gelman.diag(chains[,"kp_chain"])     
    gelman.diag(chains[,"bt_chain"])  
    fit <- circGLM(circ ~ year_s, data=sdat, burnin=200, thin=30, Q=2500)
    fit
    samples <- fit$all_chains
    samplesdf <- data.frame(samples)
    names(samplesdf) <- names(samplesdf) %>% 
      gsub("_chain", "", .)
    samplesdf %>% 
      select("b0","kp","bt") %>% 
      pivot_longer(cols=everything(), names_to="parameter", values_to="sample_value") %>% 
      ggplot() +
      geom_histogram(mapping=aes(x=sample_value, y=after_stat(density)), bins=75) +
      facet_wrap(vars(parameter), scales="free")
    
    mean(samplesdf$b0)
    mean(samplesdf$kp)
    mean(samplesdf$bt)
    
    hpdi <- function (samp, prob = 0.95) {
      vals <- sort(samp)
      nsamp <- length(vals)
      gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
      init <- 1:(nsamp - gap)
      inds <- which.min(vals[init + gap,drop=FALSE] - vals[init, drop=FALSE])
      ans <- cbind(lower=vals[inds], upper=vals[inds + gap])
      return(ans)
    }
    
    hpdi(samplesdf$b0, prob=0.95)
    hpdi(samplesdf$kp, prob=0.95)
    hpdi(samplesdf$bt, prob=0.95)
    
    quantile(samplesdf$b0, prob=c(0.025,0.975))
    quantile(samplesdf$kp, prob=c(0.025,0.975))
    quantile(samplesdf$bt, prob=c(0.025,0.975))
    
    
    year <- seq(from=min(sdat$year), to=max(sdat$year), by=1) 
    year_s <- scale(year, center=year_center, scale=year_scale)
    n <- length(year)
    results <- matrix(NA, nrow=n, ncol=5) 
    colnames(results) <- c("mnmu","mulo95","muhi95","ppdlo95","ppdhi95")
    
    for ( i in 1:n ) {
      
      mu <- samplesdf$b0 + inv_tan_half(samplesdf$bt * year_s[i])
    
      ppd <- rvon_mises(n=length(mu), mu=mu, kappa=samplesdf$kp)
      
      results[i,1] <- mean(mu)
      #results[i,2:3] <- hpdi(mu, prob=0.95)
      results[i,2:3] <- quantile(mu, prob=c(0.025,0.975)) #CPI
      #results[i,4:5] <- hpdi(ppd, prob=0.95)
      results[i,4:5] <- quantile(ppd, prob=c(0.025,0.975)) #CPI
    }
    results <- circ2day(results) #transform to day of year
    preds <- data.frame(year, year_s, results)
    rm(year, year_s, n, results, mu, ppd) #clean up
    
    coef_list[[spec[j]]] <- coef(fit)
    
    finalplot_list[[spec[j]]] <- preds %>%
      ggplot() +
      geom_ribbon(mapping=aes(x=year, ymin=mulo95, ymax=muhi95), alpha=0.2) +
      geom_point(data=sdat, 
                 mapping=aes(x=year, y=plant2$Day.of.the.year)) +
      geom_line(mapping=aes(x=year, y=mnmu)) +
      geom_line(mapping=aes(x=year, y=ppdlo95), lty=2) +
      geom_line(mapping=aes(x=year, y=ppdhi95), lty=2) + labs(title= spec,
                                                            x ="Year", y = "Julian Date")
  }
}



