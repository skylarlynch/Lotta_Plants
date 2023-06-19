#' ---
#' title: Get BioClim variables for list of lat/longs
#' author: Sierra Jech
#' date: 29 Nov 2021
#' output:
#'     github_document:
#'         pandoc_args: --webtex
#' ---

#' Purpose: Get the BioClim variables for the study sites in the database using online tutorial: https://www.worldclim.org/data/bioclim.html and C. Havrilla's meta-analysis code 
#' 
#' Bio-Climatic Variables
#' BIO1 = Annual Mean Temperature**
#' BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#' BIO3 = Isothermality (BIO2/BIO7) (×100)
#' BIO4 = Temperature Seasonality (standard deviation ×100)
#' BIO5 = Max Temperature of Warmest Month
#' BIO6 = Min Temperature of Coldest Month
#' BIO7 = Temperature Annual Range (BIO5-BIO6)
#' BIO8 = Mean Temperature of Wettest Quarter
#' BIO9 = Mean Temperature of Driest Quarter
#' BIO10 = Mean Temperature of Warmest Quarter
#' BIO11 = Mean Temperature of Coldest Quarter
#' BIO12 = Annual Precipitation **
#' BIO13 = Precipitation of Wettest Month
#' BIO14 = Precipitation of Driest Month
#' BIO15 = Precipitation Seasonality (Coefficient of Variation)
#' BIO16 = Precipitation of Wettest Quarter
#' BIO17 = Precipitation of Driest Quarter
#' BIO18 = Precipitation of Warmest Quarter
#' BIO19 = Precipitation of Coldest Quarter
#' 
library(raster) #do not load tidyr after loading this package. Both have extract() functions, so you will need to use raster::extract() if that happens
library(dismo)
library(dplyr)
library(maptools)
library(maps)
library(mapdata)
library(ggplot2)
library(rasterVis)

# Load data
# Load lat/longs without other info
#lat_long <- read.csv("data/lat_long.csv")

# # Load more data with lat_longs
# locations <- read.csv("data/biocrust_restoration_getBIOCLIM.csv")
# str(locations)
# locations$Latitude_x <- as.numeric(locations$Latitude_x)
# # Get rid of NA's
# locations <- locations %>% 
#   filter(locations$Longitude_y != "NA")
# # Get rid of NA's
# locations <- locations %>% 
#   filter(locations$Latitude_x != "NA")

# LOAD DATA DATABASE
locations <- read.csv("data/Restoration_database_long_data.csv")
str(locations)
str(locations$x)
locations$x <- as.numeric(locations$x)
str(locations$x)
# Get rid of NA's
locations <- locations %>% 
  filter(locations$y != "NA")
# Get rid of NA's
locations <- locations %>% 
  filter(locations$x != "NA")

# get distinct locations
locations_minimal <- locations %>% 
  distinct(y, x, .keep_all = TRUE)

# how many distinct locations are there?
nrow(locations_minimal)

# Map it on a map of the world
world <- map_data("world")
#head(world)
#tail(world)

# Add my points
ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  theme_bw() + 
  geom_point(data = locations_minimal, 
             mapping = aes(x = x, y = y, 
                           #color = year_experiment_initiated
                           color = country # when color is by country, you can tell there are some errors in the lat/long info
             )) +
  theme(legend.position = "none")


#Use R to extract data from WorldClim - MAT, MAP
r <- getData("worldclim",var="bio",res=10)

#Bio 1 and Bio12 are mean annual temperature and annual precipitation:
r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")

#load in locations of papers as points - extract from metalocs
lats <- c(locations$x)
longs <- c(locations$y)

#make dataframe 
coords<- data.frame(x=lats, y=longs)

# make a spatial polygons object 
points <- SpatialPoints(coords, proj4string = r@crs)

#extract climate info from BioClim dataframe
values <- extract(r,points)

#make data frame 
df <- cbind.data.frame(coordinates(points),values)
head(df)
#df

##############################
#Aridity analysis in "ClimClass"
# Carrie did not finish this code and instead used this website for aridity: https://www.uea.ac.uk/groups-and-centres/climatic-research-unit
# This analysis uses this website: https://cgiarcsi.community/2019/01/24/global-aridity-index-and-potential-evapotranspiration-climate-database-v2/
# Go to the website, download the aridity zip file and find the tiff file of values. Put it in the data folder

# Instead of using stars package, load tif as RasterStack as Carrie did in the BioClim code above
rastlist <- list.files(path = "data", pattern='.tif$', all.files=TRUE, full.names=TRUE)
aridity <- lapply(rastlist, raster)
# Access and use only the first RasterStack
aridity <- aridity[[1]]
#if not done above, load in locations of papers as points 
#lats <- c(locations_minimal$x)
#longs <- c(locations_minimal$y)
#make dataframe 
#coords<- data.frame(x=lats, y=longs)
# make a spatial polygons object 
#points <- SpatialPoints(coords, proj4string = aridity@crs)
#extract climate info from Aridity data
ai_values <- extract(aridity,points) # this works if the crs match for both points and aridity
# Aridity data needs to be rescaled
ai_values <- ai_values * 0.0001
# add to the dataframe above with Temp and Precip
df_temp <- cbind.data.frame(df, ai_values)
# check it before export 
head(df_temp)
df_temp

# Make sure data looks right 
distinct(df_temp) # across the whole database, there are 65 locations with distinct Temp and 
hist(df_temp$Temp)
df_temp$Temp <- df_temp$Temp / 10
hist(df_temp$Temp)
hist(df_temp$Prec)
hist(df_temp$ai_values)

#write df into a .csv file, if you haven't already
#write.csv(df_temp, "data/Meta_Analysis_WorldClim.csv")


##################################### Plotting Code ########################
# Get environmental variables
bioclim.data <- getData(name = "worldclim",
                        var = "bio",
                        res = 2.5,
                        path = "data/")

#res(bioclim.data$bio1)

#Global Annual Temperature Plot
# might want to re-do the scaling of the legend because it seems wrong
plot(bioclim.data$bio1/10, main ="Annual Mean Temperature (\u00B0C) ")
points(x=locations$x, y=locations$y, pch=19, cex = .5)

#Global Annual Precip Plot
# the axis is incorrect for precipitation - change the scale 
plot(bioclim.data$bio12, main ="Annual Mean Precipitation (mm)")
points(x=locations$x, y=locations$y, pch=19, cex = .5)
#levelplot(bioclim.data$bio12, par.settings = viridisTheme)
#histogram(bioclim.data$bio12)

### Plotting from Seb's Class Project
# data in this raster is stored with incorrect units, so multiply by 0.0001 so that you do not have to do that operation every time you want to use the data
#library(stars)
#img_vis  <- file.path("data/ai_et0.tif")
#sat_vis <- read_stars(img_vis, RasterIO = list(nBufXSize = 600, nBufYSize = 600))
#plot(sat_vis, axes = FALSE)

# Plot Aridity! - This one works but is not perfect
breaks <- c(0, 0.03, 0.2, 0.5, 0.65, 10)
ggplot() + 
  geom_stars(data = sat_vis) +
  scale_fill_gradientn(
    breaks = breaks,
    labels=c("", "0.03", "0.2", "0.5", "0.65", ""),
    colours = hcl.colors(6, palette = "RdYlBu"),
    values = c(0, 0.00625, 0.025, 0.05, 0.075, 1.0), # telling it what percent each color represents to help highlight the really dry areas instead of the wetter areas
    #na.value = "white",
    na.value = "transparent", 
    limits = range(breaks)
  ) + 
  #coord_cartesian() + # I should be using this but I don't remember how
  coord_equal()+
  theme_bw() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggtitle("Study Locations and Global Aridity") +
  geom_point(data = coords, 
             mapping = aes(x = x, y = y),
             pch = 1, alpha = 0.75, cex = 0.75) + 
  labs(fill = "Aridity", y = "", x = "") 
#theme(legend.position = "none")

ggsave(filename = "global_aridity.jpeg", path = "output", device='jpeg', width = 5, height = 5, dpi=700)