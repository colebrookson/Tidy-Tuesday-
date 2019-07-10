#Libraries
library(tidyverse)
library(data.table)
library(imager)
library(gridGraphics)
library(ggimage)
library(RColorBrewer)
library(cowplot)
library(grid)


#Get the data
wwc_outcomes <- readr::read_csv("https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2019/2019-07-09/wwc_outcomes.csv")
squads <- readr::read_csv("https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2019/2019-07-09/squads.csv")
codes <- readr::read_csv("https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2019/2019-07-09/codes.csv")

#Let's look at the concept of goals scored vs. goals scored against and compare the teams from each world 
#cup who had the most goals/game and the fewest goals against/game and see how often either of those teams win

#make and populate a matrix to get the information for a single game in one line
goals_both = matrix(nrow = 568, ncol = 9)

x = 1
for(yr in unique(wwc_outcomes$year)) {
  
  temp = wwc_outcomes %>% 
    filter(year == yr) %>% 
    rename(team1 = team, score1 = score, win_status1 = win_status)
  
  for(i in unique(temp$yearly_game_id)) {
    
    temp1 = temp %>% 
      filter(yearly_game_id == i) 
    
    focal_row = temp1 %>% 
      filter(team_num == 1)
    
    op_row = temp1 %>% 
      filter(team_num == 2)
    
    focal_row$team2 = op_row$team1
    focal_row$score2 = op_row$score1
    focal_row$win_status2 = op_row$win_status1
    
    op_row$team2 = focal_row$team1
    op_row$score2 = focal_row$score1
    op_row$win_status2 = focal_row$win_status1
    
    op_row = op_row %>% 
      select(-team_num)
    
    focal_row = focal_row %>% 
      select(-team_num)
    
    focal_row[2,] = op_row
    
    focal_row = as.matrix(focal_row)
    
    goals_both[x:(x+1),] = focal_row[,]
    
    x = x+2
    
  }
}

#clean up the resulting data
col_names = colnames(focal_row)
goals_both = data.frame(goals_both)
colnames(goals_both) = col_names
goals_both$score1 = as.numeric(as.character(goals_both$score1))
goals_both$score2 = as.numeric(as.character(goals_both$score2))

#extract goals for and against for the focal team
goals_for = goals_both %>% 
  group_by(year, team1) %>% 
  summarize(gf = mean(score1))
goals_against = goals_both %>% 
  group_by(year, team1) %>% 
  summarize(ga = mean(score2))

goals_for_against = goals_for
goals_for_against$ga = goals_against$ga

#figure out the goals for-against rate
goals_for_against = goals_for_against %>% 
  mutate(gfa = gf/ga)

n_distinct(goals_for_against$year)

#get the final game scores
finals = goals_both %>% 
  filter(round == 'Final') %>% 
  filter(win_status1 == 'Won')

#put in winning team or not into the for_against dataframe
goals_for_against = goals_for_against %>% 
  mutate(Won = ifelse(team1 == 'USA' & year %in% c(1991, 1999, 2015, 2019), 1,
                      ifelse(team1 == 'GER' & year %in% c(2003, 2007), 1, 
                             ifelse(team1 == 'JPN' & year == 2011, 1, 
                                    ifelse(team1 == 'NOR' & year == 1995, 1, 0)))))

#now we can do some plotting
all_gf = goals_for_against %>% 
  group_by(year) %>% 
  filter(gf == max(gf))
all_ga = goals_for_against %>% 
  group_by(year) %>% 
  filter(ga == min(ga))
all_gf$year = as.numeric(as.character(all_gf$year))
all_gf$gf = as.numeric(as.character(all_gf$gf))
all_ga$year = as.numeric(as.character(all_ga$year))
all_ga$ga = as.numeric(as.character(all_ga$ga))

#get a trophy to use as the winning points
img = load.image("C:/Users/coleb/Documents/GitHub/Tidy-Tuesday-/trophy.png")
img = rasterGrob(img, interpolate = FALSE)

#get megan rapinoe because why not
megan = load.image("C:/Users/coleb/Documents/GitHub/Tidy-Tuesday-/megan1.png")
megan = rasterGrob(megan, interpolate = FALSE)

#get a theme
theme1 <- function(){
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
    theme(axis.text.x = element_text(size = 12, color = color.axis.text)) + 
    theme(axis.text.y = element_text(size = 12, color = color.axis.text)) + 
    theme(axis.title.x = element_text(size = 18, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 18, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(color="black", size = 0.15)) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 10))+
    theme(legend.position = 'none')
}

#split the dataframes for easy plotting
all_ga_w = all_ga %>% 
  filter(Won == 1)
all_gf_w = all_gf %>% 
  filter(Won == 1)

#initial plot
ga = ggplot() +
  geom_point(data = all_ga, aes(x = year, y = ga), colour = 'red3', size = 4)

#add trophies
ga = ga + geom_blank(data = all_ga_w, aes(year, ga)) +
  mapply(function(xx, yy, ii) {
    img$name <- ii
    annotation_custom(img, xmin=xx-0.6, xmax=xx+0.6, ymin=yy-0.15, ymax=yy+0.15)},
  all_ga_w$year, all_ga_w$ga, seq_len(nrow(all_ga_w))) +
  labs(x = 'Year', y = "Lowest 'Goals Against' Per Game", title = "Fewest 'Goals Against' Per Game Each Year in WWC History") +
  theme1() +
  scale_x_continuous(breaks = c(1991, 1995, 1999, 2003, 2007, 2011, 2015, 2019))

#initial plot
gf = ggplot() +
  geom_point(data = all_gf, aes(x = year, y = gf), colour = 'red3', size = 4)

#add trophies
gf = gf + geom_blank(data = all_gf_w, aes(year, gf)) +
  mapply(function(xx, yy, ii) {
    img$name <- ii
    annotation_custom(img, xmin=xx-0.6, xmax=xx+0.6, ymin=yy-0.15, ymax=yy+0.15)},
    all_gf_w$year, all_gf_w$gf, seq_len(nrow(all_gf_w))) +
  labs(x = 'Year', y = "Highest 'Goals For' Per Game", title = "Most 'Goals For' Per Game Each Year in WWC History") +
  theme1() +
  scale_x_continuous(breaks = c(1991, 1995, 1999, 2003, 2007, 2011, 2015, 2019))

gf_ga = plot_grid(gf, ga, rel_widths = c(0.95, 1))



