ramen_ratings <- readr::read_csv("https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2019/2019-06-04/ramen_ratings.csv")

library(tidyverse)
library(scales)
library(viridis)
library(maps)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(data.table)


#quick peek at the data
summary(ramen_ratings)
n_distinct(ramen_ratings$brand); n_distinct(ramen_ratings$variety); n_distinct(ramen_ratings$style)

ramen_ratings = ramen_ratings %>%
  filter(stars != "NA")

unique(ramen_ratings$country) #so we can see there are some countries here that a) have typos and b) aren't the proper names

## let's make the map to see what it looks like then join some data to it to make our chloropleth
world = ne_countries(scale = 'medium', returnclass = 'sf')

#let's check that my countries are in the world countries
my_countries = unique(ramen_ratings$country) %in% world$sovereignt

#some of my countries aren't names properly, so let's fix that - there are misspellings, duplicates and false names
ramen_ratings = ramen_ratings %>% 
  mutate(country = replace(country, country == 'Dubai', 'United Arab Emirates')) %>% 
  mutate(country = replace(country, country == 'Holland', 'Netherlands')) %>% 
  mutate(country = replace(country, country == 'Hong Kong', 'China')) %>% 
  mutate(country = replace(country, country == 'Phlippines', 'Philippines')) %>% 
  mutate(country = replace(country, country == 'United States', 'United States of America')) %>% 
  mutate(country = replace(country, country == 'Sarawak', 'Malaysia')) %>% 
  mutate(country = replace(country, country == 'USA', 'United States of America')) %>% 
  mutate(country = replace(country, country == 'UK', 'United Kingdom'))

#now that thats done, lets see what this world map looks like
map = ggplot(data = world) +
  geom_sf()

#now, let's join our ratings data to the world data so we can colour by that data
ratings_by_c = ramen_ratings %>% 
  group_by(country) %>% 
  dplyr::summarise(Prop_Stars = mean(stars)) %>% 
  rename(sovereignt = country)

#use data.table to set keys so we can left-join
setDT(world); setkey(world, sovereignt)
setDT(ratings_by_c); setkey(ratings_by_c, sovereignt)

#perform merge 
world = merge(world, ratings_by_c, all = TRUE)

#now let's plot our map again but let's make it coloured by our mean stars category
map = ggplot(data = world) +
        geom_sf(colour = 'grey20', aes(fill = Prop_Stars)) +
        scale_fill_viridis_c(option = 'plasma', na.value = 'grey80', name = 'Mean Star Rating \nOut of Five') +
        coord_sf(ylim = c(-55,80), xlim = c(-170, 170)) + 
        theme_bw() +
        labs(title = 'Average Ramen Rating by Country', x = 'Latitude', y = 'Longitude') +
        ggplot2::theme(legend.position = c(0.071,0.267),
                       legend.background = element_rect(colour = 'black'),
                       legend.title = element_text(size = 8),
                       legend.text = element_text(size = 6.5))

