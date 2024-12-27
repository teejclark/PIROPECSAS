#### Testing to see if there's spatial bias, and where to clip extent

library(tidyverse)
library(ggplot2)
library(lubridate)
library(rnaturalearth)
library(sf)
library(adehabitatHR)
library(ggmap)
library(patchwork)

###############################################################################################################################

#### Test ECSAS/PIROP data
ECSAS <- readRDS("Data/RData/all_data_nontest_surveys.RDS")$distdata %>% 
  mutate(label = "ECSAS")

# ECSAS <- readRDS("./Data/RData/Seabird_detections_nontest_surveys_transects.RDS") %>%
#   mutate(label = "ECSAS")
PIROP <- readRDS("./Data/RData/all_PIROP_data.RDS")
PIROP <- PIROP$flatfiledata
PIROP$label <- "PIROP"

# plot obs.lat and obs.long of watches, and see what it looks like?
# create spatial data frame
spatial_data <- ECSAS %>%
  filter(English == "Northern Gannet") %>%
  dplyr::select(LongStart, LatStart, label)

PIROP_sp <- PIROP %>%
  filter(English == "Northern Gannet") %>%
  dplyr::select(LongStart, LatStart, label)

spatial_data <- spatial_data %>%
  add_row(PIROP_sp) %>%
  {SpatialPointsDataFrame(cbind(.$LongStart, .$LatStart), .,
                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))}
  
# graph
ne_countries(scale = "medium", returnclass = "sf") %>%
  st_geometry() %>%
  plot(xlim = c(-106, 0), ylim = c(30,80),
       axes = T, main = "Northern Gannet")

# create date format for data
ECSAS$month <- as.factor(lubridate::month(ECSAS$Date))

ECSAS %>%
  filter(English == "Northern Gannet") %>%
  ggplot(aes(x = month))+
  geom_bar()

###############################################################################################################################

#### NORTHERN GANNETS

# combine PIROP/ECSAS
spatial_data <- ECSAS %>%
  filter(English == "Northern Gannet") %>%
  dplyr::select(LongStart, LatStart, label)

PIROP_sp <- PIROP %>%
  filter(English == "Northern Gannet") %>%
  dplyr::select(LongStart, LatStart, label)

spatial_data <- spatial_data %>%
  add_row(PIROP_sp) %>%
  {SpatialPointsDataFrame(cbind(.$LongStart, .$LatStart), .,
                          proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))}

# plot
plot(spatial_data, axes = T)
  
# STEP 1: cut off data using 47W and 41N
# NOTE: visually chose to do this at 40W and 40N
spatial_data1 <- subset(spatial_data, LongStart <= -40)
spatial_data1 <- subset(spatial_data1, LatStart >= 40)

# graph again
spatial_data1 %>%
  st_as_sf(coords = c("LongStart", "LatStart")) %>%
  ggplot(aes(color = label)) + geom_sf()

# STEP 2: constrain by 95% minimum polygon
# add animal id
spatial_data2 <- spatial_data1
spatial_data2$id <- as.factor(rep(1, nrow(spatial_data2)))

# calculate mcp
mcp <- mcp(spatial_data2[,4], percent = 95)

# graph with google
map <-
  get_map(
    location = c(-75, 40, -40, 65),
    source = "stadia",
    maptype =  "stamen_terrain_background"
  )

a <- ggmap(map)+
  geom_polygon(data = fortify(mcp), aes(long, lat), alpha = 0.3, color = "black", lwd = 1.25)+
  geom_point(data = data.frame(spatial_data2), aes(x = LongStart, y = LatStart, color = label), alpha = 0.7)+
  ggtitle("Northern Gannets") #+
  #theme(legend.position="none")
a  
# clip data by object
spatial_data3 <- st_as_sf(spatial_data2)
mcp_sf_gannet <- st_as_sf(mcp)
spatial_data3 <- st_intersection(spatial_data3, mcp_sf_gannet)

ggmap(map)+
  geom_polygon(data = fortify(mcp), aes(long, lat), alpha = 0.3, color = "black", lwd = 1.25)+
  geom_point(data = data.frame(spatial_data3), aes(x = LongStart, y = LatStart, color = rev(label)), alpha = 0.7)+
  ggtitle("Northern Gannets")+
theme(legend.position="none")

# save polygon to clip by
st_write(mcp_sf_gannet, "mcp_sf_gannet.shp", delete_layer = TRUE)
  
# look at time
time <- ECSAS %>%
  filter(English == "Northern Gannet") %>%
  dplyr::select(Date) %>%
  mutate(label = "ECSAS")

PIROP_time <- PIROP %>%
  filter(English == "Northern Gannet") %>%
  dplyr::select(Date) %>%
  mutate(label = "PIROP")
  
time %>%
  add_row(PIROP_time) %>%
  mutate(month = lubridate::month(Date)) %>%
  ggplot(aes(x=month, fill = label), ) +
  geom_bar()


#### FIGURE 2 - reproduce figure for paper!
# NOTE: purple is PIROP, orange is ECSAS

# space

map <- get_map(location = c(-70, 40, -48, 57),
               source = "stadia",
               maptype =  "stamen_terrain_background")


noga_col.sf <- st_read(here::here("GIS/Shapefiles/NOGA_colonies.shp"))
noga_col <- data.frame(long = st_coordinates(noga_col.sf)[, 1],
                       lat =  st_coordinates(noga_col.sf)[, 2])
a <- ggmap(map)+
  geom_polygon(data = fortify(mcp), aes(long, lat), alpha = 0, color = "black", lwd = 1.15)+
  geom_point(data = data.frame(spatial_data3), aes(x = LongStart, y = LatStart, 
                                                   # not sure why I had to change this to get colors right
                                                   # color = rev(label)), alpha = 1)+
                                                   color = label), alpha = 1)+
  geom_point(data = noga_col, aes(long, lat), color = "Yellow", size = 4) +
  # scale_color_manual(values = c("#3d348b","#f35b04"))+
  # DAF Not sure why bu colors were now coming out wrong
  scale_color_manual(values = c("#f35b04", "#3d348b"))+
  ggsn::scalebar(location = "bottomleft",
                 x.min = -54, x.max = -51,
                 y.min = 41, y.max = 42.5,
                 transform = TRUE, dist_unit = "km",
                 dist = 200, height = .2, 
                 st.dist = .3, st.size = 3)+
  ggspatial::annotation_north_arrow(location = "tr", which_north = "grid",
                                    height = unit(.75, "cm"), width = unit(.75, "cm"))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_classic()+
  theme(
    legend.position = "none",
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"))
a

# time
b <- time %>%
  add_row(PIROP_time) %>%
  mutate(month = lubridate::month(Date, label = T, abbr = T)) %>%
  ggplot(aes(x=month, fill = label)) +
  geom_bar(position = "dodge", stat = "count")+
  scale_fill_manual("Program", values = c("#f35b04","#3d348b"))+
  scale_x_discrete(breaks = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  xlab("Month")+
  ylab("Count")+
  theme_classic()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 18))

a+b + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 32))

ggsave("figure2_alt.png", width = 14, height = 8, dpi = 600)

ggsave("figure2_alt_alt.png", width = 8, height = 7, dpi = 600)

# add graph of effort - check with Dave
effort <- ECSAS %>%
  filter(English == "Northern Gannet") %>%
  dplyr::select(Date, Effort) %>%
  mutate(label = "ECSAS")

PIROP_effort <- PIROP %>%
  filter(English == "Northern Gannet") %>%
  dplyr::select(Date, Effort) %>%
  mutate(label = "PIROP")

a <- effort %>%
  add_row(PIROP_effort) %>%
  #mutate(month = lubridate::month(Date)) %>%
  mutate(month = lubridate::month(Date, label = T, abbr = T)) %>%
  group_by(month, label) %>%
  summarize(total_effort = sum(Effort, na.rm = T)) %>%
  ggplot(aes(x = month, y = total_effort, fill = label)) +
  geom_bar(stat = "identity")+
  scale_fill_manual("Program", values = c("#f35b04","#3d348b"))+
  scale_x_discrete(breaks = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  xlab("Month")+
  ylab("Effort (km)")+
  theme_classic()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    legend.position = "none")

b <- effort %>%
  add_row(PIROP_effort) %>%
  #mutate(month = lubridate::month(Date)) %>%
  #mutate(month = lubridate::month(Date, label = T, abbr = T)) %>%
  mutate(year = lubridate::year(Date)) %>%
  group_by(year, label) %>%
  summarize(total_effort = sum(Effort, na.rm = T)) %>%
  ggplot(aes(x = year, y = total_effort, fill = label)) +
  geom_bar(stat = "identity")+
  scale_fill_manual("Program", values = c("#f35b04","#3d348b"))+
  #scale_x_discrete(breaks = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  xlab("Year")+
  ylab("Effort (km)")+
  theme_classic()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=24),
    axis.text = element_text(size=18, color="black"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 18))

