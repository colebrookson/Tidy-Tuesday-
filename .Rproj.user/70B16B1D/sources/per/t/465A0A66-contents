library(tidyverse)
library(cowplot)
library(bbmle)
library(lubridate)
library(rsample)
library(tibble)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(MuMIn)
library(here)
library(mapdata)
library(maps)
library(maptools)
library(ggmap)
library(ggrepel)
library(raster)
library(ggthemes)
library(ggsn)
library(rgeos)
library(rgdal)

#test <- read_csv("C:/Users/brookson/Salmon_Work/SalmonWork-master/Hakai_lice_data_CB_edits.csv")

#make vars into factors
mainlice$year <- as.factor(mainlice$year);mainlice$collection <- as.factor(mainlice$collection)

## Figure 1: Map
salmonsites <- read_csv('Hakai_sampling_site_coordinates.csv')
names(salmonsites)


BC_shp = readOGR('C:/Users/brookson/Documents/GitHub/SalmonWork/SpatialData/COAST_TEST2.shp')
coords <- data.frame(cbind(salmonsites$Xlatitude,salmonsites$Ylongitutde))
colnames(coords) <- c('lat', 'long')

sort(unique(ggplot2::map_data("world")$region))
westcoast <- map_data('canada')

#make two themes, one for the inset, one for the big map
fte_theme_map_small <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_rect(colour = 'black')) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 12, color = color.axis.text, angle = 90)) + 
    theme(axis.text.y = element_text(size = 12, color = color.axis.text)) + 
    theme(axis.title.x = element_text(size = 14, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 14, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(color="black", size = 0.15)) 
}
fte_theme_map_big <- function(){
  color.background = 'grey75'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill = 'white', color = 'white')) +
    theme(plot.background = element_rect(fill=color.background,color = color.background)) +
    theme(panel.border = element_rect(colour = 'black')) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_blank()) + 
    theme(axis.text.y = element_blank()) + 
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(plot.title = element_blank()) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(color="black", size = 0.15)) 
}
#create the two maps
states    <- c('Washington')
provinces <- c("British Columbia")

us <- getData("GADM",country="USA",level=1)
canada <- getData("GADM",country="CAN",level=1)

us.states <- us[us$NAME_1 %in% states,]
ca.provinces <- canada[canada$NAME_1 %in% provinces,]

biggermap = ggplot()+
  geom_polygon(data = us.states,aes(x=long,y=lat,group=group), colour = 'black', size = 0.01, fill = 'grey75')+
  geom_polygon(data=ca.provinces, aes(x=long,y=lat,group=group), colour = 'black', size = 0.01, fill = 'grey75')+
  coord_cartesian(xlim = c(-128.5,-119.5), ylim = c(48.1,51.0)) +
  fte_theme_map_big()+
  annotate("rect", xmin = -127.3, xmax = -125, ymin = 49.8, ymax = 50.9, alpha = .7)+
  annotate('text', x = -121, y = 50.7, label = 'British Columbia', size = 4)+
  annotate('text', x = -120.7, y = 48.5, label = 'Washington', size = 4)+
  annotate('text', x = -126.7, y = 48.5, label = 'Vancouver Island', size = 4)+
  annotate('segment',x=-126.7, y=48.6, xend=-125.5, yend=49.5, arrow=arrow(length = unit(0.04, "npc")), 
           alpha = 0.8, size=1.1, color="black")
biggermap

smallermap = ggplot()+
  geom_polygon(data=ca.provinces,aes(x=long,y=lat,group=group), colour = 'black', size = 0.01, fill = 'grey75')+
  coord_cartesian(xlim = c(-127.0,-125), ylim = c(50,50.75))+
  geom_point(data = coords, aes(long,lat), color = 'black', size = 4, shape = 21, fill = 'red3')+
  fte_theme_map_small()+
  labs(x = 'Longitude', y = 'Latitude')+
  annotate("rect", xmin = -125.51, xmax = -125.05, ymin = 50.1, ymax = 50.5, alpha = .65)+
  annotate("rect", xmin = -126.9, xmax = -126.55, ymin = 50.45, ymax = 50.7, alpha = .65)+
  annotate('text', x= -125.8, y = 50.33, label = 'Discovery Islands', size = 4)+
  annotate('text', x = -126.6, y = 50.38, label = 'Johnstone Strait', size = 4)+
  annotate('segment',x=-125.8, y=50.302, xend=-125.565, yend=50.27, arrow=arrow(length = unit(0.04, "npc")), 
           alpha = 0.8, size=1.1, color="black")+
  annotate('segment',x=-126.6, y=50.4, xend=-126.65, yend=50.44, arrow=arrow(length = unit(0.04, "npc")), 
           alpha = 0.8, size=1.1, color="black")
smallermap

#make the one plot inset with the other
insetmap = ggdraw()+
  draw_plot(smallermap) + 
  draw_plot(biggermap, x=0.105, y=0.161, width=0.5, height=0.3) 
insetmap
ggsave('study_map.png', plot = insetmap,
       width = 8, height = 7.5,
       dpi = 300)

## Figure 2: 3x2 of sal and lice species
fig_2_data = mainlice %>% 
  group_by(year, spp) %>% 
  summarize(mean_cal = mean(all.cal), mean_lep = mean(all.leps))
lep_y_title = expression(paste("Mean Motile ", italic("L. salmonis"), " Per Fish"))
cal_y_title = expression(paste("Mean Motile ", italic("C. clemensi"), " Per Fish"))

#make each plot then stitch them together
fig_2_1 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'PI')), aes(x = year, y = mean_lep), 
           colour = 'black', fill = '#ff9999', width = 0.7) +
  labs(x = 'Year', y = lep_y_title, title = 'Pink') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 0.45))
fig_2_2 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'CU')), aes(x = year, y = mean_lep), 
           colour = 'black', fill = '#59AE7F', width = 0.7) +
  labs(x = 'Year', y = lep_y_title, title = 'Chum') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0.0, 0.45))
fig_2_3 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'SO')), aes(x = year, y = mean_lep), 
           colour = 'black', fill = '#23359d', width = 0.7) +
  labs(x = 'Year', y = lep_y_title, title = 'Sockeye') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0.0, 0.45))
fig_2_4 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'PI')), aes(x = year, y = mean_cal), 
           colour = 'black', fill = '#ff9999', width = 0.7) +
  labs(x = 'Year', y = cal_y_title, title = 'Pink') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 0.85))
fig_2_5 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'CU')), aes(x = year, y = mean_cal), 
           colour = 'black', fill = '#59AE7F', width = 0.7) +
  labs(x = 'Year', y = cal_y_title, title = 'Chum') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 0.85)) 
fig_2_6 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'SO')), aes(x = year, y = mean_cal), 
           colour = 'black', fill = '#23359d', width = 0.7) +
  labs(x = 'Year', y = cal_y_title, title = 'Sockeye') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 0.85)) 

fig_2 = plot_grid(fig_2_1, fig_2_2, fig_2_3,
                  fig_2_4, fig_2_5, fig_2_6, nrow = 2, rel_heights = c(0.85, 1), rel_widths = c(1, 0.8, 0.8))
ggsave('lice_per_fish_sp.png', plot = fig_2,
       width = 8, height = 7.5,
       dpi = 300)

## Figure 3: average number of lice per fish (fish sp and by year)
fig_3_data = mainlice %>% 
  group_by(year, spp) %>% 
  summarize(mean_lice = mean(all.lice))

fig_3_1 = ggplot() +
  geom_col(data = (fig_3_data %>% filter(spp == 'PI')), aes(x = year, y = mean_lice), 
           colour = 'black', fill = '#ff9999', width = 0.7) +
  labs(x = 'Year', y = 'Mean Number of Motile Lice Per Fish', title = 'Pink') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 1.15))
fig_3_2 = ggplot() +
  geom_col(data = (fig_3_data %>% filter(spp == 'CU')), aes(x = year, y = mean_lice), 
           colour = 'black', fill = '#59AE7F', width = 0.7) +
  labs(x = 'Year', y = cal_y_title, title = 'Chum') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 1.15)) 
fig_3_3 = ggplot() +
  geom_col(data = (fig_3_data %>% filter(spp == 'SO')), aes(x = year, y = mean_lice), 
           colour = 'black', fill = '#23359d', width = 0.7) +
  labs(x = 'Year', y = cal_y_title, title = 'Sockeye') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 1.15))

fig_3 = plot_grid(fig_3_1, fig_3_2, fig_3_3, nrow = 1,
                  rel_widths = c(1, 0.8, 0.8))
ggsave('lice_per_fish.png', plot = fig_3,
       scale = 1, width = 8, height = 6.7,
       dpi = 300)

#get the necessary model selection things:
lepmod.crossed <- glmmTMB(all.leps ~ spp * site.region + spp * year + 
                            site.region * year + (1 | collection), 
                          data = mainlice, family=nbinom2)
calmod.crossed <- glmmTMB(all.cal ~ spp * site.region + spp * year + 
                            site.region * year + (1 | collection), 
                          data = mainlice, family=nbinom2)
lepmod.crossed_dredge = MuMIn::dredge(lepmod.crossed, subset = (`cond(site.region)` && `cond(year)`))
calmod.crossed_dredge = MuMIn::dredge(calmod.crossed, subset = (`cond(site.region)` && `cond(year)`))

#so the goal here is to bootstrap the data (parametric) by resampling the hierarchical levels, then 
#run the model averaging process with the new data and use that to get our model-averaged CI's. Essentially, 
#we're wrapping our esitmator of mu in a function and bootstrapping it

bootlice =  mainlice %>% 
  dplyr::select(all.cal, all.leps, spp, site.region, collection, year, ufn)

sock2015D <- bootlice %>% filter(spp == 'SO' & year == '2015' & site.region == 'D')
sock2015J <- bootlice %>% filter(spp == 'SO' & year == '2015' & site.region == 'J')
sock2016D <- bootlice %>% filter(spp == 'SO' & year == '2016' & site.region == 'D')
sock2016J <- bootlice %>% filter(spp == 'SO' & year == '2016' & site.region == 'J')
sock2017D <- bootlice %>% filter(spp == 'SO' & year == '2017' & site.region == 'D')
sock2017J <- bootlice %>% filter(spp == 'SO' & year == '2017' & site.region == 'J')
sock2018D <- bootlice %>% filter(spp == 'SO' & year == '2018' & site.region == 'D')
sock2018J <- bootlice %>% filter(spp == 'SO' & year == '2018' & site.region == 'J')
sock2019D <- bootlice %>% filter(spp == 'SO' & year == '2019' & site.region == 'D')
sock2019J <- bootlice %>% filter(spp == 'SO' & year == '2019' & site.region == 'J')

chum2015D <- bootlice %>% filter(spp == 'CU' & year == '2015' & site.region == 'D')
chum2015J <- bootlice %>% filter(spp == 'CU' & year == '2015' & site.region == 'J')
chum2016D <- bootlice %>% filter(spp == 'CU' & year == '2016' & site.region == 'D')
chum2016J <- bootlice %>% filter(spp == 'CU' & year == '2016' & site.region == 'J')
chum2017D <- bootlice %>% filter(spp == 'CU' & year == '2017' & site.region == 'D')
chum2017J <- bootlice %>% filter(spp == 'CU' & year == '2017' & site.region == 'J')
chum2018D <- bootlice %>% filter(spp == 'CU' & year == '2018' & site.region == 'D')
chum2018J <- bootlice %>% filter(spp == 'CU' & year == '2018' & site.region == 'J')
chum2019D <- bootlice %>% filter(spp == 'CU' & year == '2019' & site.region == 'D')
chum2019J <- bootlice %>% filter(spp == 'CU' & year == '2019' & site.region == 'J')

pink2015D <- bootlice %>% filter(spp == 'PI' & year == '2015' & site.region == 'D')
pink2015J <- bootlice %>% filter(spp == 'PI' & year == '2015' & site.region == 'J')
pink2016D <- bootlice %>% filter(spp == 'PI' & year == '2016' & site.region == 'D')
pink2016J <- bootlice %>% filter(spp == 'PI' & year == '2016' & site.region == 'J')
pink2017D <- bootlice %>% filter(spp == 'PI' & year == '2017' & site.region == 'D')
pink2017J <- bootlice %>% filter(spp == 'PI' & year == '2017' & site.region == 'J')
pink2018D <- bootlice %>% filter(spp == 'PI' & year == '2018' & site.region == 'D')
pink2018J <- bootlice %>% filter(spp == 'PI' & year == '2018' & site.region == 'J')
pink2019D <- bootlice %>% filter(spp == 'PI' & year == '2019' & site.region == 'D')
pink2019J <- bootlice %>% filter(spp == 'PI' & year == '2019' & site.region == 'J')

bootintervalcal = matrix(nrow = 30, ncol = 1000)
bootintervallep = matrix(nrow = 30, ncol = 1000)

#Run the Loop to populate the matrix - CAUTION:: THIS TAKES OVER 2 FULL DAYS TO RUN ON A POWERFUL DESKTOP COMPUTER
pb = txtProgressBar(min = 0, max = 1000, initial = 0) 
start_time <- Sys.time()
for(i in 1:1000) {
  sock2015Dboot = matrix(nrow = nrow(sock2015D), ncol = 7)
  n = 1
  for(k in unique(sock2015D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2015D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2015Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2015Jboot = matrix(nrow = nrow(sock2015J), ncol = 7)
  n = 1
  for(k in unique(sock2015J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2015J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2015Jboot[rows,] = res
    n = n + nrow(res)
  }
  sock2016Dboot = matrix(nrow = nrow(sock2016D), ncol = 7)
  n = 1
  for(k in unique(sock2016D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2016D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2016Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2016Jboot = matrix(nrow = nrow(sock2016J), ncol = 7)
  n = 1
  for(k in unique(sock2016J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2016J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2016Jboot[rows,] = res
    n = n + nrow(res)
  }
  sock2017Dboot = matrix(nrow = nrow(sock2017D), ncol = 7)
  n = 1
  for(k in unique(sock2017D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2017D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2017Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2017Jboot = matrix(nrow = nrow(sock2017J), ncol = 7)
  n = 1
  for(k in unique(sock2017J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2017J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2017Jboot[rows,] = res
    n = n + nrow(res)
  }
  sock2018Dboot = matrix(nrow = nrow(sock2018D), ncol = 7)
  n = 1
  for(k in unique(sock2018D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2018D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2018Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2018Jboot = matrix(nrow = nrow(sock2018J), ncol = 7)
  n = 1
  for(k in unique(sock2018J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2018J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2018Jboot[rows,] = res
    n = n + nrow(res)
  }  
  sock2019Dboot = matrix(nrow = nrow(sock2019D), ncol = 7)
  n = 1
  for(k in unique(sock2019D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2019D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2019Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2019Jboot = matrix(nrow = nrow(sock2019J), ncol = 7)
  n = 1
  for(k in unique(sock2019J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2019J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2019Jboot[rows,] = res
    n = n + nrow(res)
  } 
  
  #pink
  
  chum2015Dboot = matrix(nrow = nrow(chum2015D), ncol = 7)
  n = 1
  for(k in unique(chum2015D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2015D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2015Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2015Jboot = matrix(nrow = nrow(chum2015J), ncol = 7)
  n = 1
  for(k in unique(chum2015J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2015J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2015Jboot[rows,] = res
    n = n + nrow(res)
  }
  chum2016Dboot = matrix(nrow = nrow(chum2016D), ncol = 7)
  n = 1
  for(k in unique(chum2016D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2016D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2016Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2016Jboot = matrix(nrow = nrow(chum2016J), ncol = 7)
  n = 1
  for(k in unique(chum2016J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2016J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2016Jboot[rows,] = res
    n = n + nrow(res)
  }
  chum2017Dboot = matrix(nrow = nrow(chum2017D), ncol = 7)
  n = 1
  for(k in unique(chum2017D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2017D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2017Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2017Jboot = matrix(nrow = nrow(chum2017J), ncol = 7)
  n = 1
  for(k in unique(chum2017J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2017J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2017Jboot[rows,] = res
    n = n + nrow(res)
  }
  chum2018Dboot = matrix(nrow = nrow(chum2018D), ncol = 7)
  n = 1
  for(k in unique(chum2018D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2018D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2018Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2018Jboot = matrix(nrow = nrow(chum2018J), ncol = 7)
  n = 1
  for(k in unique(chum2018J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2018J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2018Jboot[rows,] = res
    n = n + nrow(res)
  }  
  chum2019Dboot = matrix(nrow = nrow(chum2019D), ncol = 7)
  n = 1
  for(k in unique(chum2019D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2019D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2019Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2019Jboot = matrix(nrow = nrow(chum2019J), ncol = 7)
  n = 1
  for(k in unique(chum2019J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2019J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2019Jboot[rows,] = res
    n = n + nrow(res)
  } 
  
  #pink
  
  pink2015Dboot = matrix(nrow = nrow(pink2015D), ncol = 7)
  n = 1
  for(k in unique(pink2015D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2015D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2015Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2015Jboot = matrix(nrow = nrow(pink2015J), ncol = 7)
  n = 1
  for(k in unique(pink2015J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2015J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2015Jboot[rows,] = res
    n = n + nrow(res)
  }
  pink2016Dboot = matrix(nrow = nrow(pink2016D), ncol = 7)
  n = 1
  for(k in unique(pink2016D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2016D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2016Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2016Jboot = matrix(nrow = nrow(pink2016J), ncol = 7)
  n = 1
  for(k in unique(pink2016J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2016J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2016Jboot[rows,] = res
    n = n + nrow(res)
  }
  pink2017Dboot = matrix(nrow = nrow(pink2017D), ncol = 7)
  n = 1
  for(k in unique(pink2017D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2017D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2017Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2017Jboot = matrix(nrow = nrow(pink2017J), ncol = 7)
  n = 1
  for(k in unique(pink2017J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2017J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2017Jboot[rows,] = res
    n = n + nrow(res)
  }
  pink2018Dboot = matrix(nrow = nrow(pink2018D), ncol = 7)
  n = 1
  for(k in unique(pink2018D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2018D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2018Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2018Jboot = matrix(nrow = nrow(pink2018J), ncol = 7)
  n = 1
  for(k in unique(pink2018J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2018J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2018Jboot[rows,] = res
    n = n + nrow(res)
  }  
  pink2019Dboot = matrix(nrow = nrow(pink2019D), ncol = 7)
  n = 1
  for(k in unique(pink2019D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2019D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2019Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2019Jboot = matrix(nrow = nrow(pink2019J), ncol = 7)
  n = 1
  for(k in unique(pink2019J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2019J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2019Jboot[rows,] = res
    n = n + nrow(res)
  }  
  
  #bind the matrices so we can have our resampled dataframe
  bootdata = data.frame(rbind(sock2015Dboot,sock2015Jboot,sock2016Dboot,sock2016Jboot,sock2017Dboot,sock2017Jboot,sock2018Dboot,sock2018Jboot,sock2019Dboot,sock2019Jboot,
                              chum2015Dboot,chum2015Jboot,chum2016Dboot,chum2016Jboot,chum2017Dboot,chum2017Jboot,chum2018Dboot,chum2018Jboot,chum2019Dboot,chum2019Jboot,
                              pink2015Dboot,pink2015Jboot,pink2016Dboot,pink2016Jboot,pink2017Dboot,pink2017Jboot,pink2018Dboot,pink2018Jboot,pink2019Dboot,pink2019Jboot)) %>% 
    rename(all.cal = X1, all.leps = X2, spp = X3, site.region = X4, collection = X5, year = X6, ufn = X7)
  bootdata$all.cal = as.integer(as.character(bootdata$all.cal))
  bootdata$all.leps = as.integer(as.character(bootdata$all.leps))
  
  #now run our set of models       
  #omited the last two models in each set since their weights were 0. 
  lep1 = glmmTMB(all.leps ~ site.region + year + spp + 
                   site.region * year + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  lep2 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * year + site.region * year + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  lep3 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * site.region + site.region * year + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  lep4 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * site.region + spp * year + site.region * year + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  lep5 = glmmTMB(all.leps ~ site.region + year + spp + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  lep6 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * year +  
                   (1 | collection), data = bootdata, family = nbinom2) 
  lep7 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * site.region + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  lep8 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * site.region + spp * year +  
                   (1 | collection), data = bootdata, family = nbinom2) 
  
  cal1 = glmmTMB(all.cal ~ site.region + year + spp + 
                   site.region * year + site.region * spp +
                   (1 | collection), data = bootdata, family = nbinom2) 
  cal2 = glmmTMB(all.cal ~ site.region + year + spp + 
                   site.region * year + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  cal3 = glmmTMB(all.cal ~ site.region + year + spp + 
                   spp * site.region + site.region * year + spp * year + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  cal4 = glmmTMB(all.cal ~ site.region + year + spp + 
                   spp * year + site.region * year + 
                   (1 | collection), data = bootdata, family = nbinom2)
  cal5 = glmmTMB(all.cal ~ site.region + year + spp + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  cal6 = glmmTMB(all.cal ~ site.region + year + spp + 
                   site.region * year +  
                   (1 | collection), data = bootdata, family = nbinom2)
  cal7 = glmmTMB(all.cal ~ site.region + year + spp + 
                   spp * year + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  cal8 = glmmTMB(all.cal ~ site.region + year + spp + 
                   spp * site.region + spp * year +  
                   (1 | collection), data = bootdata, family = nbinom2) 
  
  #get the predictions of the estiamtes
  lep1pred <- ggpredict(lep1, terms = c('spp', 'year', 'site.region'))
  lep2pred <- ggpredict(lep2, terms = c('spp', 'year', 'site.region')) 
  lep3pred <- ggpredict(lep3, terms = c('spp', 'year', 'site.region')) 
  lep4pred <- ggpredict(lep4, terms = c('spp', 'year', 'site.region')) 
  lep5pred <- ggpredict(lep5, terms = c('spp', 'year', 'site.region')) 
  lep6pred <- ggpredict(lep6, terms = c('spp', 'year', 'site.region')) 
  lep7pred <- ggpredict(lep7, terms = c('spp', 'year', 'site.region')) 
  lep8pred <- ggpredict(lep8, terms = c('spp', 'year', 'site.region')) 
  
  cal1pred <- ggpredict(cal1, terms = c('spp', 'year', 'site.region'))
  cal2pred <- ggpredict(cal2, terms = c('spp', 'year', 'site.region')) 
  cal3pred <- ggpredict(cal3, terms = c('spp', 'year', 'site.region')) 
  cal4pred <- ggpredict(cal4, terms = c('spp', 'year', 'site.region')) 
  cal5pred <- ggpredict(cal5, terms = c('spp', 'year', 'site.region')) 
  cal6pred <- ggpredict(cal6, terms = c('spp', 'year', 'site.region')) 
  cal7pred <- ggpredict(cal7, terms = c('spp', 'year', 'site.region')) 
  cal8pred <- ggpredict(cal8, terms = c('spp', 'year', 'site.region')) 
  
  ###start by getting them all in one dataframe with the weights
  
  #pull the predicted values from each one
  lepallpred = data.frame(cbind(lep1pred$predicted, lep2pred$predicted, lep3pred$predicted, lep4pred$predicted,
                                lep5pred$predicted, lep6pred$predicted, lep7pred$predicted, lep8pred$predicted)) %>% 
    rename(lep1 = X1, lep2 = X2, lep3 = X3, lep4 = X4, lep5 = X5, lep6 = X6, lep7 = X7, lep8 = X8)
  calallpred = data.frame(cbind(cal1pred$predicted, cal2pred$predicted, cal3pred$predicted, cal4pred$predicted,
                                cal5pred$predicted, cal6pred$predicted, cal7pred$predicted, cal8pred$predicted)) %>% 
    rename(cal1 = X1, cal2 = X2, cal3 = X3, cal4 = X4, cal5 = X5, cal6 = X6, cal7 = X7, cal8 = X8)
  
  #add the weights from the model selection object
  lepallpred = lepallpred %>% 
    mutate(w1 = rep(lepmod.crossed_dredge$weight[1], nrow(lepallpred)), 
           w2 = rep(lepmod.crossed_dredge$weight[2], nrow(lepallpred)),
           w3 = rep(lepmod.crossed_dredge$weight[3], nrow(lepallpred)),
           w4 = rep(lepmod.crossed_dredge$weight[4], nrow(lepallpred)),
           w5 = rep(lepmod.crossed_dredge$weight[5], nrow(lepallpred)),
           w6 = rep(lepmod.crossed_dredge$weight[6], nrow(lepallpred)),
           w7 = rep(lepmod.crossed_dredge$weight[7], nrow(lepallpred)),
           w8 = rep(lepmod.crossed_dredge$weight[8], nrow(lepallpred)))
  
  calallpred = calallpred %>% 
    mutate(w1 = rep(calmod.crossed_dredge$weight[1], nrow(calallpred)), 
           w2 = rep(calmod.crossed_dredge$weight[2], nrow(calallpred)),
           w3 = rep(calmod.crossed_dredge$weight[3], nrow(calallpred)),
           w4 = rep(calmod.crossed_dredge$weight[4], nrow(calallpred)),
           w5 = rep(calmod.crossed_dredge$weight[5], nrow(calallpred)),
           w6 = rep(calmod.crossed_dredge$weight[6], nrow(calallpred)),
           w7 = rep(calmod.crossed_dredge$weight[7], nrow(calallpred)),
           w8 = rep(calmod.crossed_dredge$weight[8], nrow(calallpred)))
  
  #now make averaged predictions!
  lepallpred = lepallpred %>% 
    mutate(lep1w = lep1*w1, lep2w = lep2*w2, lep3w = lep3*w3, lep4w = lep4*w4, lep5w = lep5*w5, lep6w = lep6*w6, lep7w = lep7*w7, lep8w = lep8*w8) %>% 
    mutate(avg = lep1w + lep2w + lep3w + lep4w + lep5w + lep6w + lep7w + lep8w)
  
  calallpred = calallpred %>% 
    mutate(cal1w = cal1*w1, cal2w = cal2*w2, cal3w = cal3*w3, cal4w = cal4*w4, cal5w = cal5*w5, cal6w = cal6*w6, cal7w = cal7*w7, cal8w = cal8*w8) %>% 
    mutate(avg = cal1w + cal2w + cal3w + cal4w + cal5w + cal6w + cal7w + cal8w)
  
  #keep just the averaged predictions and the relevant grouping info 
  lepavgpred = lepallpred %>% 
    dplyr::select(avg) %>% 
    mutate(sal = lep1pred$x, reg = lep1pred$facet, yr = lep1pred$group)
  lepavgpred$sal = factor(lepavgpred$sal, levels = c(1, 2, 3), labels = c('CU', 'PI', 'SO'))
  
  calavgpred = calallpred %>% 
    dplyr::select(avg) %>% 
    mutate(sal = cal1pred$x, reg = cal1pred$facet, yr = cal1pred$group)
  calavgpred$sal = factor(calavgpred$sal, levels = c(1, 2, 3), labels = c('CU', 'PI', 'SO'))
  
  bootintervalcal[,i] = calavgpred$avg
  bootintervallep[,i] = lepavgpred$avg
  
  setTxtProgressBar(pb,i)
}
end_time <- Sys.time()
end_time - start_time

boot_int_cal = data.frame(bootintervalcal)
boot_int_lep = data.frame(bootintervallep)
write_csv(boot_int_cal, 'boot_int_cal.csv')
write_csv(boot_int_lep, 'boot_int_lep.csv')

#####NOTE: to make life easy, I just manually go and change these to long version in the csv itself - have to do that 
#before readng in the next lines

#pull the data
interval_cal_long = read_csv('boot_int_cal_long.csv')
interval_lep_long = read_csv('boot_int_lep_long.csv')

#name the columns, sort them, and transpose them
names_lep = lepavgpred %>% 
  unite(., col = 'names', sal:yr, sep = '_')
names_cal = calavgpred %>% 
  unite(., col = 'names', sal:yr, sep = '_')

names_lep = as.vector(names_lep$names)
names_cal = as.vector(names_cal$names)

colnames(interval_cal_long) = names_cal
colnames(interval_lep_long) = names_lep

interval_cal_long_sorted <- apply(interval_cal_long,2,sort,decreasing=F)
interval_lep_long_sorted <- apply(interval_lep_long,2,sort,decreasing=F)

upci_cal = interval_cal_long_sorted[995, ]
loci_cal = interval_cal_long_sorted[5, ]
upci_lep = interval_lep_long_sorted[995, ]
loci_lep = interval_lep_long_sorted[5, ]

#fit the models and do the model averaging process again to get the actual estimates themselves
lep1 = glmmTMB(all.leps ~ site.region + year + spp + 
                 site.region * year + 
                 (1 | collection), data = mainlice, family = nbinom2) 
lep2 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * year + site.region * year + 
                 (1 | collection), data = mainlice, family = nbinom2) 
lep3 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * site.region + site.region * year + 
                 (1 | collection), data = mainlice, family = nbinom2) 
lep4 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * site.region + spp * year + site.region * year + 
                 (1 | collection), data = mainlice, family = nbinom2) 
lep5 = glmmTMB(all.leps ~ site.region + year + spp + 
                 (1 | collection), data = mainlice, family = nbinom2) 
lep6 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * year +  
                 (1 | collection), data = mainlice, family = nbinom2) 
lep7 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * site.region + 
                 (1 | collection), data = mainlice, family = nbinom2) 
lep8 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * site.region + spp * year +  
                 (1 | collection), data = mainlice, family = nbinom2) 

cal1 = glmmTMB(all.cal ~ site.region + year + spp + 
                 site.region * year + site.region * spp +
                 (1 | collection), data = mainlice, family = nbinom2) 
cal2 = glmmTMB(all.cal ~ site.region + year + spp + 
                 site.region * year + 
                 (1 | collection), data = mainlice, family = nbinom2) 
cal3 = glmmTMB(all.cal ~ site.region + year + spp + 
                 spp * site.region + site.region * year + spp * year + 
                 (1 | collection), data = mainlice, family = nbinom2) 
cal4 = glmmTMB(all.cal ~ site.region + year + spp + 
                 spp * year + site.region * year + 
                 (1 | collection), data = mainlice, family = nbinom2)
cal5 = glmmTMB(all.cal ~ site.region + year + spp + 
                 (1 | collection), data = mainlice, family = nbinom2) 
cal6 = glmmTMB(all.cal ~ site.region + year + spp + 
                 site.region * year +  
                 (1 | collection), data = mainlice, family = nbinom2)
cal7 = glmmTMB(all.cal ~ site.region + year + spp + 
                 spp * year + 
                 (1 | collection), data = mainlice, family = nbinom2) 
cal8 = glmmTMB(all.cal ~ site.region + year + spp + 
                 spp * site.region + spp * year +  
                 (1 | collection), data = mainlice, family = nbinom2) 

#get the predictions of the estiamtes
lep1pred <- ggpredict(lep1, terms = c('spp', 'year', 'site.region'))
lep2pred <- ggpredict(lep2, terms = c('spp', 'year', 'site.region')) 
lep3pred <- ggpredict(lep3, terms = c('spp', 'year', 'site.region')) 
lep4pred <- ggpredict(lep4, terms = c('spp', 'year', 'site.region')) 
lep5pred <- ggpredict(lep5, terms = c('spp', 'year', 'site.region')) 
lep6pred <- ggpredict(lep6, terms = c('spp', 'year', 'site.region')) 
lep7pred <- ggpredict(lep7, terms = c('spp', 'year', 'site.region')) 
lep8pred <- ggpredict(lep8, terms = c('spp', 'year', 'site.region')) 

cal1pred <- ggpredict(cal1, terms = c('spp', 'year', 'site.region'))
cal2pred <- ggpredict(cal2, terms = c('spp', 'year', 'site.region')) 
cal3pred <- ggpredict(cal3, terms = c('spp', 'year', 'site.region')) 
cal4pred <- ggpredict(cal4, terms = c('spp', 'year', 'site.region')) 
cal5pred <- ggpredict(cal5, terms = c('spp', 'year', 'site.region')) 
cal6pred <- ggpredict(cal6, terms = c('spp', 'year', 'site.region')) 
cal7pred <- ggpredict(cal7, terms = c('spp', 'year', 'site.region')) 
cal8pred <- ggpredict(cal8, terms = c('spp', 'year', 'site.region')) 

###start by getting them all in one dataframe with the weights

#pull the predicted values from each one
lepallpred = data.frame(cbind(lep1pred$predicted, lep2pred$predicted, lep3pred$predicted, lep4pred$predicted,
                              lep5pred$predicted, lep6pred$predicted, lep7pred$predicted, lep8pred$predicted)) %>% 
  rename(lep1 = X1, lep2 = X2, lep3 = X3, lep4 = X4, lep5 = X5, lep6 = X6, lep7 = X7, lep8 = X8)
calallpred = data.frame(cbind(cal1pred$predicted, cal2pred$predicted, cal3pred$predicted, cal4pred$predicted,
                              cal5pred$predicted, cal6pred$predicted, cal7pred$predicted, cal8pred$predicted)) %>% 
  rename(cal1 = X1, cal2 = X2, cal3 = X3, cal4 = X4, cal5 = X5, cal6 = X6, cal7 = X7, cal8 = X8)

#add the weights from the model selection object
lepallpred = lepallpred %>% 
  mutate(w1 = rep(lepmod.crossed_dredge$weight[1], nrow(lepallpred)), 
         w2 = rep(lepmod.crossed_dredge$weight[2], nrow(lepallpred)),
         w3 = rep(lepmod.crossed_dredge$weight[3], nrow(lepallpred)),
         w4 = rep(lepmod.crossed_dredge$weight[4], nrow(lepallpred)),
         w5 = rep(lepmod.crossed_dredge$weight[5], nrow(lepallpred)),
         w6 = rep(lepmod.crossed_dredge$weight[6], nrow(lepallpred)),
         w7 = rep(lepmod.crossed_dredge$weight[7], nrow(lepallpred)),
         w8 = rep(lepmod.crossed_dredge$weight[8], nrow(lepallpred)))

calallpred = calallpred %>% 
  mutate(w1 = rep(calmod.crossed_dredge$weight[1], nrow(calallpred)), 
         w2 = rep(calmod.crossed_dredge$weight[2], nrow(calallpred)),
         w3 = rep(calmod.crossed_dredge$weight[3], nrow(calallpred)),
         w4 = rep(calmod.crossed_dredge$weight[4], nrow(calallpred)),
         w5 = rep(calmod.crossed_dredge$weight[5], nrow(calallpred)),
         w6 = rep(calmod.crossed_dredge$weight[6], nrow(calallpred)),
         w7 = rep(calmod.crossed_dredge$weight[7], nrow(calallpred)),
         w8 = rep(calmod.crossed_dredge$weight[8], nrow(calallpred)))

#now make averaged predictions!
lepallpred = lepallpred %>% 
  mutate(lep1w = lep1*w1, lep2w = lep2*w2, lep3w = lep3*w3, lep4w = lep4*w4, lep5w = lep5*w5, lep6w = lep6*w6, lep7w = lep7*w7, lep8w = lep8*w8) %>% 
  mutate(avg = lep1w + lep2w + lep3w + lep4w + lep5w + lep6w + lep7w + lep8w)

calallpred = calallpred %>% 
  mutate(cal1w = cal1*w1, cal2w = cal2*w2, cal3w = cal3*w3, cal4w = cal4*w4, cal5w = cal5*w5, cal6w = cal6*w6, cal7w = cal7*w7, cal8w = cal8*w8) %>% 
  mutate(avg = cal1w + cal2w + cal3w + cal4w + cal5w + cal6w + cal7w + cal8w)

#keep just the averaged predictions and the relevant grouping info 
lepavgpred = lepallpred %>% 
  dplyr::select(avg) %>% 
  mutate(sal = lep1pred$x, reg = lep1pred$facet, yr = lep1pred$group)

calavgpred = calallpred %>% 
  dplyr::select(avg) %>% 
  mutate(sal = cal1pred$x, reg = cal1pred$facet, yr = cal1pred$group)

#put the up and lo CI's into the df's 
calavgpred$conf.high = upci_cal
calavgpred$conf.low = loci_cal
lepavgpred$conf.high = upci_lep
lepavgpred$conf.low = loci_lep

#Make the plots
fte_theme1 <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 12, color = color.axis.text)) + 
    theme(axis.text.y = element_text(size = 12, color = color.axis.text)) + 
    theme(axis.title.x = element_text(size = 14, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 14, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(color="black", size = 0.15)) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 12))+
    theme(legend.position = c(0.8,0.82),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
} 

leg_title <- 'Salmon Species'
lep_y_effects = expression(paste("Mean Number of Motile  ", italic("L. salmonis"), " Per Fish"))
cal_y_effects = expression(paste("Mean Number of Motile  ", italic("C. clemensi"), " Per Fish"))
lepsmodplot_avg <- lepavgpred %>% 
  group_by(., yr,sal,reg) %>% 
  ggplot(aes(x = sal, y = avg, colour = sal, shape = reg)) +
  scale_shape_manual(values = c(15,17), labels = c('Discovery Islands', 'Johnstone Strait')) +
  geom_errorbar(aes(ymin=conf.low, ymax = conf.high,width = 0), position = position_dodge(width = 0.8),colour = 'Black')+
  geom_point(size = 4,position = position_dodge(width = 0.8)) +
  facet_wrap(~yr,nrow=1,strip.position = "bottom")+
  theme(strip.background = element_blank(), strip.placement = "outside") + 
  scale_color_manual(leg_title,values=c('#59AE7F', '#ff9999', '#23359d'), labels = c('Chum', 'Pink', 'Sockeye'))+
  labs(x = 'Salmon Species/Year', y = lep_y_effects) +
  guides(shape = guide_legend(title = 'Region', override.aes = list(shape = c(0,2)), type = 'b')) +
  fte_theme1()
lepsmodplot_avg
ggsave('model_ests_lep.png', plot = lepsmodplot_avg,
       width = 8, height = 7.5,
       dpi = 300)
calmodplot_avg <- calavgpred %>% 
  group_by(., yr,sal,reg) %>% 
  ggplot(aes(x = sal, y = avg, colour = sal, shape = reg)) +
  scale_shape_manual(values = c(15,17), labels = c('Discovery Islands', 'Johnstone Strait')) +
  geom_errorbar(aes(ymin=conf.low, ymax = conf.high,width = 0), position = position_dodge(width = 0.8),colour = 'Black')+
  geom_point(size = 4,position = position_dodge(width = 0.8)) +
  facet_wrap(~yr,nrow=1,strip.position = "bottom")+
  theme(strip.background = element_blank(), strip.placement = "outside") + 
  scale_colour_manual(leg_title,values=c('#59AE7F', '#ff9999', '#23359d'), labels = c('Chum', 'Pink', 'Sockeye'))+
  labs(x = 'Salmon Species/Year', y = cal_y_effects) +
  guides(shape = guide_legend(title = 'Region', override.aes = list(shape = c(0,2)), type = 'b')) +
  fte_theme1()
calmodplot_avg

ggsave('model_ests_cal.png', plot = calmodplot_avg,
       width = 8, height = 7.5,
       dpi = 300)


mainlice %>% 
  filter(spp == PI) %>% 
  filter(year == 2019) %>% 
  summarize(all = mean(all.cal))

































