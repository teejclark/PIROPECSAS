#### Analyzing ECSAS data... ####

library(tidyverse)
library(ggplot2)
library(Distance)
library(lubridate)
library(sf)
library(here)

# load data
the.data <- readRDS(here("Data/RData/all_data_nontest_surveys.RDS"))

# Create observations for detection function fitting
# data <- readRDS("./Data/RData/Seabird_detections_nontest_surveys_transects.RDS") %>%
data <- the.data$distdata %>% 
  filter(DistanceCode == "A" | DistanceCode == "B" |
           DistanceCode == "C" | DistanceCode == "D") %>%
  filter(InTransect == "-1") %>%
  rename(distance = Distance) %>%
  rename(size = Count)

data$year <- year(data$Date)
data$month <- month(data$Date)
data <- arrange(data, Date)

# other set-up
# distance bins
dbins <- c(0, 50, 100, 200, 300)

mcp_sf_gannet <- st_read("mcp_sf_gannet.shp") %>% st_set_crs("+proj=longlat +datum=WGS84 +no_defs")

# Create tr_table for dht.
transects <- st_as_sf(the.data$transects, coords = c("LongStart", "LatStart"), 
                      crs = "+proj=longlat +datum=WGS84 +no_defs") %>% 
  st_intersection(mcp_sf_gannet) %>% 
  st_drop_geometry() %>%
  mutate(Region.Label = as.factor(year(.$Date))) %>% 
  select(Sample.Label, Region.Label, Effort)

# duplicate tables and add ddf identifier
tr_table <- rbind(transects, transects)
tr_table$ddf_id <- c(rep(1, nrow(transects)),
                     rep(2, nrow(transects)))
tr_table$Effort <- tr_table$Effort/2

###############################################################################################################################

#### GANNETS ####

mcp_sf_gannet <- st_read("mcp_sf_gannet.shp") %>% st_set_crs("+proj=longlat +datum=WGS84 +no_defs")

gannet <- data %>%
  filter(English == "Northern Gannet")

# clip spatial data
gannet <- st_as_sf(gannet, coords = c("LongStart", "LatStart"), 
                  crs = "+proj=longlat +datum=WGS84 +no_defs")
gannet <- st_intersection(gannet, mcp_sf_gannet)
gannet <- st_drop_geometry(gannet)

table(gannet$year)

# set-up
# unique ID for each observation
gannet$object <- seq(1, nrow(gannet), 1)

# set strata
gannet$Region.Label <- as.factor(year(gannet$Date))

# just density estimates, convert units
gannet$Area <- 0
uf_gannet <- unflatten(gannet)
cu <- convert_units("metre","kilometre","square kilometer")

# split up detection function
gannet_F <- subset(gannet, FlySwim == "F")
gannet_W <- subset(gannet, FlySwim == "W")

# fit detection functions without size...
gannet.hnF.nosize <- ds(data=gannet_F, key = "hn", adjustment = "cos", nadj = 2, cutpoints = dbins) # AIC: 8691.32
plot(gannet.hnF.nosize)
gannet.hnW.nosize <- ds(data=gannet_W, key = "hn", adjustment = "cos", nadj = 2, cutpoints = dbins) # AIC: 3055.512
plot(gannet.hnW.nosize)

# first get g(x) (i.e., probability of detection conditional on distance)
gx.F <- mrds::detfct(gannet_F$distance, gannet.hnF.nosize$ddf$ds$aux$ddfobj, width = 300)
gx.W <- mrds::detfct(gannet_W$distance, gannet.hnW.nosize$ddf$ds$aux$ddfobj, width = 300)

# do a regression
# 
plot(gx.F, log(gannet_F$size))
plot(gx.W, log(gannet_W$size))
reg.F <- lm(log(gannet_F$size)~gx.F)
reg.W <- lm(log(gannet_W$size)~gx.W)

summary(reg.F)
summary(reg.W)

# now estimate E(s) -- the expected group size for both groups
reg_sig.F <- summary(reg.F)$sigma
n.F <- length(gannet_F$size)
gbar.F <- sum(gx.F)/n.F
varz.F <- reg_sig.F^2 * (1 + 1/n.F + (1-gbar.F)/sum((gx.F - gbar.F)^2))
Es.F <- exp(coef(reg.F)[1] + coef(reg.F)[2] + varz.F/2)

reg_sig.W <- summary(reg.W)$sigma
n.W <- length(gannet_W$size)
gbar.W <- sum(gx.W)/n.W
varz.W <- reg_sig.W^2 * (1 + 1/n.W + (1-gbar.W)/sum((gx.W - gbar.W)^2))
Es.W <- exp(coef(reg.W)[1] + coef(reg.W)[2] + varz.W/2)

# compare this to the mean group size: slightly different with cos adj in ddf
Es.F;Es.W # 1.509065 ; 1.314151

mean(gannet_F$size)
mean(gannet_W$size)

# calculate vasize# calculate variance of group size (12/16/22)
varEs.F <- exp(2*(coef(reg.F)[1] + coef(reg.F)[2]) + varz.F) * (1+varz.F/2) * varz.F/n.F
varEs.W <- exp(2*(coef(reg.W)[1] + coef(reg.W)[2]) + varz.W) * (1+varz.W/2) * varz.W/n.W

# replace size column in data.frame with estimated average size
gannet2 <- gannet %>%
  mutate(size = case_when(FlySwim == "F" ~ Es.F,
                          FlySwim == "W" ~ Es.W))

# split up detection function
gannet2_F <- subset(gannet2, FlySwim == "F")
gannet2_W <- subset(gannet2, FlySwim == "W")

# Choose final model with cos adjustments
gannet.hnF.nosize <- ds(data=gannet2_F, key = "hn", adjustment = "cos", nadj = 2, cutpoints = dbins) # AIC: 8482.488 
plot(gannet.hnF.nosize)
gannet.hnW.nosize <- ds(data=gannet2_W, key = "hn", adjustment = "cos", nadj = 2, cutpoints = dbins) # AIC: 2013.442

uf_gannet2 <- unflatten(gannet2)

# with group regression
gannet_dht2b <- dht2(ddf = list(gannet.hnF.nosize$ddf, gannet.hnW.nosize$ddf),
                     observations = uf_gannet2$obs.table,
                     transects = tr_table, geo_strat = uf_gannet2$region.table,
                     stratification = "geographical",
                     strat_formula = ~Region.Label,
                     convert_units = cu,
                     sample_fraction = 0.5)

