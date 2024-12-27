#### Calculating ECSAS/PIROP conversion factor

library(tidyverse)
library(ggplot2)
library(Distance)
library(lubridate)
library(sf)
library(sp)
library(here)

#### GANNETS ####
the.data <- readRDS(here("Data/RData/all_data_test_surveys.RDS"))

# data <- readRDS("./Data/RData/Seabird_detections_test_surveys_transects.RDS")
data <- the.data$distdata
mcp_sf_gannet <- st_read("mcp_sf_gannet.shp") %>% st_set_crs("+proj=longlat +datum=WGS84 +no_defs")

# in transect if = -1, out of transect if = 0

# add in spatial stuff
library(ggmap)
spatial_data <- data %>%
  filter(English == "Northern Gannet") %>%
  dplyr::select(LongStart, LatStart)

spatial_data <- spatial_data %>%
  {SpatialPointsDataFrame(cbind(.$LongStart, .$LatStart), .,
                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))}

spatial_data1 <- subset(spatial_data, LongStart <= -40)
spatial_data1 <- subset(spatial_data1, LatStart >= 40)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

spatial_data2 <- spatial_data1 %>%
  st_as_sf(coords = c("LongStart", "LatStart")) 

ggplot(data = world) +
  geom_sf() +
  geom_point(data = spatial_data2, aes(x = LongStart, y = LatStart))+
  coord_sf(xlim = c(-70, -48), ylim = c(40, 57), expand = FALSE)

# sort by date, and create year and month variables
data$year <- year(data$Date)
data$month <- month(data$Date)
data <- arrange(data, Date)

# clip data frame for gannets
data2 <- st_as_sf(data, coords = c("LongStart", "LatStart"), 
                  crs = "+proj=longlat +datum=WGS84 +no_defs")
data3 <- st_intersection(data2, mcp_sf_gannet)
data3 <- st_drop_geometry(data3)

gannet_test_PIROP <- data3 %>%
  filter(English == "Northern Gannet")

table(gannet_test_PIROP$DistanceCode, useNA = "a")

# Just need to remove observations on other side of ship ("T"), (and keep
# all others at all distances,) although there aren't any "T"s in the  current data
# set, adding this anyway for future. 
gannet_test_PIROP <- gannet_test_PIROP %>% 
  filter(DistanceCode != "T" ) %>% 
  # filter(DistanceCode == "A" | DistanceCode == "B" |
  #          DistanceCode == "C" | DistanceCode == "D" |
  #          DistanceCode == "E") %>% # include all estimated distances
  rename(distance = Distance)%>%
  rename(size = Count)

gannet_test_ECSAS <- gannet_test_PIROP %>%
  filter(DistanceCode == "A" | DistanceCode == "B" |
           DistanceCode == "C" | DistanceCode == "D" | DistanceCode == "3") %>%
  filter(InTransect == "-1")
table(gannet_test_ECSAS$DistanceCode, useNA = "a")
table(gannet_test_ECSAS$distance, useNA = "a")

# set up distance bins
dbins_ECSAS <- c(0, 50, 100, 200, 300)

# Create tr_table to get the effort for dht. The effort is the same
# regardless of species of PIROP vs ECSAS - it's the same transects so done once
# here.
transects <- st_as_sf(the.data$transects, coords = c("LongStart", "LatStart"), 
                      crs = "+proj=longlat +datum=WGS84 +no_defs") %>% 
  st_intersection(mcp_sf_gannet) %>%
  st_drop_geometry() %>%
  mutate(Region.Label = as.factor(year(.$Date))) %>% 
  dplyr::select(Sample.Label, Region.Label, Effort)

# duplicate tables and add ddf identifier for doing two detection functions
tr_table <- rbind(transects, transects) %>%
  mutate(ddf_id = c(rep(1, nrow(transects)),
                    rep(2, nrow(transects))),
         Effort = Effort / 2)

# need unique ID for each observation
gannet_test_ECSAS$object <- seq(1, nrow(gannet_test_ECSAS), 1)

# for flatfile, need "Region.Label" and "Area" for the strata
gannet_test_ECSAS$Region.Label <- as.factor(year(gannet_test_ECSAS$Date))

# 0 => just density estimates
gannet_test_ECSAS$Area <- 0
uf_test_ECSAS <- unflatten(gannet_test_ECSAS)

# unit conversion, since distances are in m, but Effort is in km
cu <- convert_units("metre","kilometre","square kilometer")

# split data up for detection functions
# NOTE: don't need to split for PIROP data
gannet_splitF_ECSAS <- subset(gannet_test_ECSAS, FlySwim == "F")
gannet_splitW_ECSAS <- subset(gannet_test_ECSAS, FlySwim == "W")


################# Do size bias regression #################
gannet.hnF.nosize <- ds(data=gannet_splitF_ECSAS, key = "hn", adjustment = NULL, cutpoints = dbins_ECSAS) # AIC: 3291.146
plot(gannet.hnF.nosize, main = "Test Flying")
gannet.hnW.nosize <- ds(data=gannet_splitW_ECSAS, key = "hn", adjustment = NULL, cutpoints = dbins_ECSAS) # AIC: 1268.063
plot(gannet.hnW.nosize, main = "Test Water")

# look at plot of size vs detection distance
plot(gannet_splitF_ECSAS$distance, gannet_splitF_ECSAS$size)
plot(gannet_splitW_ECSAS$distance, gannet_splitW_ECSAS$size)

# first get g(x) (i.e., probability of detection conditional on distance)
gx.F <- mrds::detfct(gannet_splitF_ECSAS$distance, gannet.hnF.nosize$ddf$ds$aux$ddfobj, width = 300)
gx.W <- mrds::detfct(gannet_splitW_ECSAS$distance, gannet.hnW.nosize$ddf$ds$aux$ddfobj, width = 300)

# do a regression
# 
plot(gx.F, log(gannet_splitF_ECSAS$size))
plot(gx.W, log(gannet_splitW_ECSAS$size))
reg.F <- lm(log(gannet_splitF_ECSAS$size)~gx.F)
reg.W <- lm(log(gannet_splitW_ECSAS$size)~gx.W)

summary(reg.F)
summary(reg.W)

# now estimate E(s) -- the expected group size for both groups
reg_sig.F <- summary(reg.F)$sigma
n.F <- length(gannet_splitF_ECSAS$size)
gbar.F <- sum(gx.F)/n.F
varz.F <- reg_sig.F^2 * (1 + 1/n.F + (1-gbar.F)/sum((gx.F - gbar.F)^2))
Es.F <- exp(coef(reg.F)[1] + coef(reg.F)[2] + varz.F/2)

reg_sig.W <- summary(reg.W)$sigma
n.W <- length(gannet_splitW_ECSAS$size)
gbar.W <- sum(gx.W)/n.W
varz.W <- reg_sig.W^2 * (1 + 1/n.W + (1-gbar.W)/sum((gx.W - gbar.W)^2))
Es.W <- exp(coef(reg.W)[1] + coef(reg.W)[2] + varz.W/2)

# compare this to the mean group size:
Es.F;Es.W # 1.531668 ; 1.470021

mean(gannet_splitF_ECSAS$size)
mean(gannet_splitW_ECSAS$size)

# calculate vasize# calculate variance of group size (12/16/22)
varEs.F <- exp(2*(coef(reg.F)[1] + coef(reg.F)[2]) + varz.F) * (1+varz.F/2) * varz.F/n.F
varEs.W <- exp(2*(coef(reg.W)[1] + coef(reg.W)[2]) + varz.W) * (1+varz.W/2) * varz.W/n.W

# replace size column in data.frame with estimated average size
gannet_test_ECSAS <- gannet_test_ECSAS %>%
  mutate(size = case_when(FlySwim == "F" ~ Es.F,
                          FlySwim == "W" ~ Es.W))

# ... and resplit them
gannet_splitF_ECSAS <- subset(gannet_test_ECSAS, FlySwim == "F")
gannet_splitW_ECSAS <- subset(gannet_test_ECSAS, FlySwim == "W")

###############################################################

# Final detection function models
gannet.hnF_ECSAS <- ds(data = gannet_splitF_ECSAS, key = "hn", adjustment = "cos", nadj = 2, cutpoints = dbins_ECSAS, 
                       formula = ~ 1) #AIC 3194.035
plot(gannet.hnF_ECSAS, main = "Test Flying")

gannet.hnW_ECSAS <- ds(data=gannet_splitW_ECSAS, key = "hn", adjustment = "cos", nadj = 1, cutpoints = dbins_ECSAS,
                       formula = ~1) # AIC: 1252.757
plot(gannet.hnW_ECSAS, main = "Test Water")

# Get ECSAS densities
gannet_dht2_ECSAS <- dht2(ddf = list(gannet.hnF_ECSAS$ddf, gannet.hnW_ECSAS$ddf),
                              observations = uf_test_ECSAS$obs.table,
                              transects = tr_table, geo_strat = uf_test_ECSAS$region.table,
                              stratification = "geographical",
                              strat_formula = ~Region.Label, convert_units=cu,
                              sample_fraction = 0.5)

# Now that we have ECSAS square densities we need to get PIROP test data linear
# densities and create a conversion factor between them. Note that we want a
# conversion factor that can convert between the original PIROP and current
# ECSAS, but original PIROP surveys were double-sided, whereas the test surveys
# were only single sided.
PIROP_abund <- gannet_test_PIROP %>% 
  group_by(year) %>% 
  summarize(tot_birds = sum(size)) %>% 
  mutate(year = as.factor(year))

yearly_effort <- transects %>% 
  group_by(Region.Label) %>% 
  summarize(Effort = sum(Effort))

# The multiplication by two 2 is b/c PIROP test surveys were one-sided 
# whereas original PIROP were double-sided. So, on average, original PIROP
# surveys would have seen twice as many birds that the PIROP test surveys.
densities <- left_join(PIROP_abund, yearly_effort, by = join_by(year == Region.Label)) %>% 
  mutate(PIROP_lindens = 2 * (tot_birds/Effort) )

densities <- dplyr::select(gannet_dht2_ECSAS, Region.Label, Abundance) %>% 
  as.data.frame() %>% 
  filter(Region.Label != "Total") %>% 
  rename(year = Region.Label,
         ECSAS_sqdens = Abundance) %>% 
  full_join(densities %>% dplyr::select(year, PIROP_lindens), by = "year") %>% 
  # Create conversion factor between ECSAS surveys and ORIGINAL PIROP surveys.
  mutate(conv_fact = PIROP_lindens/ECSAS_sqdens) %>% 
  droplevels()

# compare examples graphically
densities_new <- data.frame(year = c(densities$year, densities$year),
                            density = c(densities$ECSAS_sqdens, densities$PIROP_lindens),
                            program = c(rep("ECSAS",8),rep("PIROP",8)))

ggplot(densities_new, aes(x = year, y = density, color = program, group = program)) +
  geom_line(lwd = 1.5)+
  #geom_line(aes(year, conv_fact, group = 1), color = "green", linewidth = 1.5) +
  #ggtitle("Comparison of ECSAS SQUARE densities to PIROP LINEAR densities")+
  scale_color_manual("Program", values = c("#f35b04","#3d348b"))+
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

# Compute conversion factor. This version converts PIROP linear densities directly
# to ECSAS square densities 
hist(densities$conv_fact)
mean(densities$conv_fact)
median(densities$conv_fact)

(conversion.factor <- median(densities$conv_fact))
(conversion.factor.ci <- quantile(densities$conv_fact, c(0.05, 0.95)))

#========================================================================================================
# Plot conversion factor by year
ggplot(data = densities, aes(year, conv_fact, group = 1)) +
  geom_line()

# since the equivalent ECSAS square density depends upon the magnitude of 
# the a given PIROP linear density, can use model to choose the best conversion factor
m <- lm(ECSAS_sqdens ~ PIROP_lindens, data = densities)
summary(m)

m1 <- mgcv::gam(ECSAS_sqdens ~ s(PIROP_lindens, bs = "ts", k = 8), data = densities)
summary(m1)
plot(m1)

# Show fitted relationship between PIROP linear densities and equivalent
# ECSAS square densities.
pred.data = data.frame(PIROP_lindens = seq(0, 0.5, length.out = 50))
pred.data$equivECSAS = predict(m1, newdata = pred.data)
ggplot(pred.data, aes(PIROP_lindens, equivECSAS)) + 
  geom_line()

ggplot(densities, aes(PIROP_lindens, conv_fact)) +
  geom_line() +
  ggtitle("Yearly conversion factors")

ggplot(densities, aes(ECSAS_sqdens, conv_fact)) +
  geom_line()
