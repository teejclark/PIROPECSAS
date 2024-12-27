#### Testing flat and varied detection functions for old PIROP data

library(tidyverse)
library(ggplot2)
library(Distance)
library(lubridate)
library(sf)
library(here)

#### GANNETS ####

# load data
# NOTE: use "flatfile format
the.data <- readRDS("./Data/RData/all_PIROP_data.RDS")
data <- the.data$flatfiledata # flatfile format for analyses (like ECSAS test data)
mcp_sf_gannet <- st_read("mcp_sf_gannet.shp") %>% st_set_crs("+proj=longlat +datum=WGS84 +no_defs")

data$year <- year(data$Date)
data$month <- month(data$Date)
data <- arrange(data, Date)

# clip spatial data
data2 <- st_as_sf(data, coords = c("LongStart", "LatStart"), 
                  crs = "+proj=longlat +datum=WGS84 +no_defs")
data3 <- st_intersection(data2, mcp_sf_gannet)
data3 <- st_drop_geometry(data3)
# 43789 rows

# for now, let's only look at the original, unlimited width PIROP data
PIROP <- data3 %>% filter(ProgramText == "PIROP (original - unlimited width)") %>%
  #filter(WhatCount == "1") %>%
  filter(year >= "1969" & year <= "1983") %>% # %>% filter(nWhatCount != "2") %>%
  rename(size = Count)

# let's extract only gannet data
gannet_PIROP <- PIROP %>%
  filter(English == "Northern Gannet") # 2,377 rows, 2193 rows on 29/01/2024

# Get transects inside the MCP
transects <- st_as_sf(the.data$transects, coords = c("LongStart", "LatStart"), 
                      crs = "+proj=longlat +datum=WGS84 +no_defs") %>% 
  st_intersection(mcp_sf_gannet) %>%
  st_drop_geometry() %>%
  mutate(Region.Label = as.factor(year(.$Date))) %>% 
  select(Sample.Label, Region.Label, Effort)

# Get PIROP linear density estimates by year
# Get the conversion factor
load(here("Data/Rdata/conversion_factor.rda"))

PIROP_abund <- gannet_PIROP %>% 
  group_by(year) %>% 
  summarize(tot_birds = sum(size)) %>% 
  mutate(year = as.factor(year))

PIROP_yearly_effort <- transects %>% 
  group_by(Region.Label) %>% 
  summarize(Effort = sum(Effort))

# Get ADJUSTED PIROP estimates in units of birds/km^2 on same scale as ECSAS
PIROP_estimates <- left_join(PIROP_abund, PIROP_yearly_effort, by = join_by(year == Region.Label)) %>% 
  mutate(mean_unadj = (tot_birds/Effort),
         mean = (tot_birds/Effort)/conversion.factor) %>%
  # If using gam to predict equivalent ECSAS linear density for PIROP lin dens uncomment next line
  # mean = predict(m1, newdata = data.frame(PIROP_lindens = tot_birds/Effort))) %>%
  select(year, mean, tot_birds, Effort, mean_unadj)

##################################################################################
# Method to get cv on PIROP linear densities. Use dht2 to get yearly sq densities
# using a dummy.ddf and 3 different transect widths to ensure it makes no difference
# to the cv on yearly densities. Therefore, the cv for each yearly square density estimate
# is a measure of encounter rate and group size variance. Then we use this same CV
# with the linear densities to generate LCI and UCI for the plot.
gannet_PIROP$distance <- 0 # create column with 0s because we have no distance information

# need unique ID for each observation
gannet_PIROP$object <- seq(1, nrow(gannet_PIROP), 1)

# for flatfile, need "Region.Label" and "Area" for the strata
gannet_PIROP$Region.Label <- as.factor(year(gannet_PIROP$Date))

# 0 => just density estimates
gannet_PIROP$Area <- 0

# unit conversion, since distances are in m, but Effort is in km
cu <- convert_units("metre","kilometre","square kilometer")

# observations
uf_gannet <- unflatten(gannet_PIROP)

# fit dummy detection functions for various 300m strip width
gannet.dummy.300 <- dummy_ddf(gannet_PIROP, width = 300) 

# get cv's
gannet_PIROP_300 <- dht2(ddf = gannet.dummy.300, flatfile = gannet_PIROP,
                         observations = uf_gannet$obs.table,
                         stratification = "geographical",
                         strat_formula = ~Region.Label,
                         convert_units = cu)
gannet_PIROP_300

# Get uncertainty in PIROP estimates
PIROP_estimates <- PIROP_estimates %>%
  mutate(cv = head(gannet_PIROP_300$Abundance_CV, - 1),
    se = mean * cv,
    LCI = lognormal.ci(mean, cv)$lcl,
    UCI = lognormal.ci(mean, cv)$ucl
    )

# plot estimates
PIROP_estimates %>%
  ggplot(aes(year, mean, group = 1)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(year, ymin = LCI, ymax = UCI), linetype = 2, fill = "black", alpha = 0.4)

# gannet_dht2
load("gannet_ECSAS_results_280224_fixedsize.Rdata")

# Make sure transect width is correct
gannet_dht2b$Area/gannet_dht2b$Effort

ECSAS_estimates <- data.frame(year = as.numeric(as.character(gannet_dht2b$Region.Label[1:length(gannet_dht2b$Region.Label)-1])),
                              mean = gannet_dht2b$Abundance[1:length(gannet_dht2b$Abundance)-1],
                              LCI = gannet_dht2b$LCI[1:length(gannet_dht2b$Abundance)-1],
                              UCI = gannet_dht2b$UCI[1:length(gannet_dht2b$Abundance)-1])

overall_estimates <- data.frame(year = seq(1969,2021,1),
                                PIROP_mean = c(PIROP_estimates$mean, rep(NA,length(seq(1984,2021)))),
                                PIROP_LCI = c(PIROP_estimates$LCI, rep(NA,length(seq(1984,2021)))),
                                PIROP_UCI = c(PIROP_estimates$UCI, rep(NA,length(seq(1984,2021)))),
                                ECSAS_mean = c(rep(NA,length(seq(1969,2005))), ECSAS_estimates$mean),
                                ECSAS_LCI = c(rep(NA,length(seq(1969,2005))), ECSAS_estimates$LCI),
                                ECSAS_UCI = c(rep(NA,length(seq(1969,2005))), ECSAS_estimates$UCI))

# max value for y axis
ymax <- max(overall_estimates$ECSAS_UCI, na.rm = TRUE) * 1.05

overall_estimates %>%
  ggplot()+
  geom_line(aes(year,PIROP_mean), color = "black", linewidth = 1.5)+
  geom_ribbon(aes(year, ymin = PIROP_LCI,
                  ymax = PIROP_UCI), linetype = 2, fill = "black", alpha = 0.4)+
  geom_line(aes(year, ECSAS_mean), color = "red", linewidth = 1.5)+
  geom_ribbon(aes(year, ymin = ECSAS_LCI, ymax = ECSAS_UCI), linetype = 2, fill = "red", alpha = 0.4)+
  ylim(0,ymax)+xlab("Year")+ylab("Density") +
  ggtitle("ADJUSTED")

#### Calculate Rate of Increase

# first fit and graph linear model
fit <- data.frame(year = overall_estimates$year,
                  mean = c(overall_estimates$PIROP_mean[1:37],
                           overall_estimates$ECSAS_mean[38:53]),
                  Program = c(rep("PIROP", 37), rep("ECSAS", 16)))

summary(lm(fit$mean ~ fit$year))

ggplot(fit, aes(year, mean))+
  geom_point()+
  geom_smooth(method = "lm")

# generate predicted values and calculate increase off of this?
fitlm <- lm(fit$mean ~ fit$year)
predlm <- predict.lm(fitlm, newdata = data.frame(fit$year),
                     type = "response")

(predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)]
mean((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)])
quantile((predlm[2:length(predlm)]-predlm[1:(length(predlm)-1)])/predlm[1:(length(predlm)-1)], c(0.05,0.95))

# Same as above - just simpler to understand
numerator <- rep(coef(fitlm)[2], length(predlm) - 1)
denominator <- head(predlm, length(predlm) - 1)
growth.rates <- numerator/denominator
growth.rates
hist(growth.rates)
(med <- round(median(growth.rates), 4))
mean(growth.rates)
(quant <- round(quantile(growth.rates, c(0.05,0.95)), 4))


a <- overall_estimates %>%
  ggplot()+
  geom_smooth(aes(fit$year, fit$mean), method = "lm", color = "black")+
  geom_line(aes(year,PIROP_mean/conversion.factor,color = "PIROP"), linewidth = 1.5)+
  geom_ribbon(aes(year, ymin = PIROP_LCI/conversion.factor, ymax = PIROP_UCI/conversion.factor), linetype = 2, fill = "#3d348b", alpha = 0.4)+
  geom_line(aes(year, ECSAS_mean,color = "ECSAS"), linewidth = 1.5)+
  geom_ribbon(aes(year, ymin = ECSAS_LCI, ymax = ECSAS_UCI), linetype = 2, fill = "#f35b04", alpha = 0.4)+
  scale_color_manual(name = "Program",
                     values = c("PIROP" = "#3d348b", "ECSAS" = "#f35b04"),
                     labels = c("ECSAS","PIROP"))+
  #ylim(0, ymax)+
  xlab("Year")+
  ylab(expression(Density ~ (individual~birds~per~km^2)))+
  theme_classic()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 18))
a

