########################## 1. Load the needed packages #########################
rm(list=ls())

library(data.table) # For data wrangling and aggregation 
library(ggplot2) # For plotting data
library(leaflet) # For mapping
library(mice)
library(readxl)
library(janitor) #clean_names
library(viridis)
library(sf)
library(rgdal)
library(sp)
library(stringr)
library(tmap)    # for static and interactive maps
library(gifski)  # for animated maps
library(PostcodesioR) # for postcode to lon+lat conversion
library(fuzzyjoin) # for fuzzy string match
library(stringdist)
library(ggspatial)
library(gridExtra) #for arranging ggplot objects
library(grid)      #for grid.rect: Draw rectangles
library(dplyr)
library(geepack)
library(gee)
library(tsibble)
library(splitstackshape)
library(stargazer) #for formatting outcome tables
library(biscale)   #for bivariate choropleth maps
library(cowplot)   #for the ggdraw function: Set up a drawing layer on top of a ggplot
library(viridisLite) #for colour scales
library(scales)   #for pretty_breaks function
library(colourvalues)  #maps viridis colours (by default) to values 
library(mapview)   #for the mapshot() function
library(Hmisc)    #for rcorr() function
library(MuMIn)    #for the model.sel command with geepack models
library(DescTools)  #for the function UncertCoef
library(ggcorrplot) #for the model.matrix function
library(targeted)   #for constructing the constant risk differences
library(brm)       #for fitting Binary Regression Model
library(knitr)     #for the kable() function
library(kableExtra) #for formatting output tables
library(spdep)      #for Moran I calculation
library(data.table) # For data wrangling and aggregation 
library(ggplot2) # For plotting data
library(janitor) #clean_names
library(stringr)
library(tsibble)
library(cowplot)   #for the ggdraw function: Set up a drawing layer on top of a ggplot
library(plm)
library(AER)
library(microsynth)
library(survey)
library(imputeTS) # for na_locf: Missing Value Imputation by Last Observation Carried Forward
library(survminer)
library(survival)
library(MatchIt)
library(Matching)
library(lubridate)
library(biostat3)
library(tidyr)    # for the unnest function
library(margins)


# File path to project folders
a<-"/Sharing Folder a"
b<-"/Sharing Folder b"

########################## 2. Load the vaccine bus site data ########################
# load the bus locations
bus_location <- read_xlsx(normalizePath(file.path(a,"Diary.xlsx")),sheet = "Sheet1") %>%
  as.data.table() %>%
  clean_names(.,"lower_camel")


bus_location[grepl("\\b(?:[A-Z]{1,2}[0-9]{1,2}\\s[0-9][A-Z]{2})\\b",venue),
             postcd:=regmatches(venue, regexpr("\\b(?:[A-Z]{1,2}[0-9]{1,2}\\s[0-9][A-Z]{2})\\b",venue))]


bus_location[grepl("\\b(?:[L][1-2 4-9]{1,2}\\s[0-9][A-Z]{2})\\b",postcd),
             day:= sequence(NROW(date))]
bus_location[day>11,day:= NA]


bus_location$date <- trimws(bus_location$date, "left", "\\w")
bus_location$date<-trimws(bus_location$date)
bus_location$date<-gsub("th|st|rd|nd", "", bus_location$date)
bus_location$date<-gsub("Aug|Augu","August",bus_location$date)
bus_location$date<-paste(bus_location$date,"2021",sep=" ")
bus_location$date<-format(as.Date(bus_location$date,format = "%d %B %Y"), "%B %d %Y")
class(bus_location$date)
table(bus_location$date)

# St. Helens
st.helens <- read_xlsx(normalizePath(file.path(a,"st.helens.xlsx")),sheet = "Sheet1") %>%
  as.data.table() %>%
  clean_names(.,"lower_camel")
st.helens$date <- format(as.Date(st.helens$date,format = "%d/%m/%Y"), "%B %d %Y")

bus_location <- rbind(bus_location, st.helens[,.(date, venue, postcd)], fill=TRUE)
bus_location$date<-as.Date(bus_location$date,format ="%B %d %Y")
class(bus_location$date)

# CWC - Cheshire West and Chester
cwc <- read_xlsx(normalizePath(file.path(a,"CheshireWestAndChester.xlsx")),
                 sheet = "Sheet1", range = "A4:F24") %>%
  as.data.table() %>%
  clean_names(.,"lower_camel")
cwc$date <- as.Date(cwc$date,format = "%Y-%m-%d")

bus_location <- rbind(bus_location, cwc[,.(date, venue, postcd, popUpClinic,vaccine)], 
                      fill=TRUE)

NROW(unique(bus_location[date<=as.Date("2021-06-30")]))
NROW(unique(bus_location[date<=as.Date("2021-06-30")]$postcd))

for (i in 1:nrow(bus_location)) {
  bus_location$eastings[i] <- postcode_lookup(bus_location$postcd[i])[3]
  bus_location$northings[i] <- postcode_lookup(bus_location$postcd[i])[4]
  bus_location$longitude[i] <- postcode_lookup(bus_location$postcd[i])[7]
  bus_location$latitude[i] <- postcode_lookup(bus_location$postcd[i])[8]
  bus_location$lsoa_code[i] <- postcode_lookup(bus_location$postcd[i])[33]
  bus_location$msoa_code[i] <- postcode_lookup(bus_location$postcd[i])[34]
  bus_location$lau2_code[i] <- postcode_lookup(bus_location$postcd[i])[35]
}


bus_location <- st_as_sf(bus_location, coords = c("longitude", "latitude"), crs = 4326)
bus_location %<>%
  group_by(postcd) %>%
  mutate(freq=n())

##only switch on the next line to exclude the pop-in sites for the sensitivity test
#bus_location <- bus_location[is.na(bus_location$popUpClinic),]

## check distribution of the date of the visits
ggplot(bus_location[bus_location$date<=as.Date("2021-06-28"),], aes(x = date)) +
  geom_histogram(binwidth=.5) +
  labs(title = "Distribution Box Plot by Date",
     x = "Date",
     y = "Frequency")


############################ 3. Load the lsoa centroids ########################

lsoa_cen <- read_sf(normalizePath(file.path(a,"lsoa_Population_Weighted_Centroids.shp")))

lsoa_cen <- lsoa_cen[grepl(paste(c("Liverpool","Cheshire East",
            "Cheshire West and Chester","Halton","Knowsley","Sefton",
            "St. Helens","Warrington","Wirral"), collapse="|"),lsoa_cen$lsoa11nm),]

lsoa_cen <- st_transform(lsoa_cen, 4326)

###################### 4. Calculate the distance buffer ########################
#set the time limit for the first dose
bus_dist <- as.data.table(st_distance(bus_location[bus_location$date>="2021-04-12"&bus_location$date<="2021-06-28",]$geometry,
                                      lsoa_cen$geometry))
bus_dist<- bus_dist[, lapply(.SD, function(x) as.numeric(x))]
colnames(bus_dist) = make.names(lsoa_cen$lsoa11cd, unique=TRUE)


bus_site0<-as.data.table(bus_location)

bus_site0[, id:=1:.N]
bus_site0<-cbind(bus_site0[date>="2021-04-12"& date<="2021-06-28", 
                         .(date,venue,postcd, eastings, northings,
                                               site_lsoa_code=lsoa_code)], bus_dist)

bus_site0 <- melt(bus_site0, id=c("date","venue","postcd", "eastings", "northings","site_lsoa_code"),
                 variable.name = "exp_lsoa_code")

## check distribution of the distance between bus locations and lsoa centriods
ggplot(bus_site0, aes(x = value)) +
  geom_histogram(binwidth=.5) +
  labs(title = "Distribution Box Plot by Distance",
       x = "Distance",
       y = "Frequency")

# select those sites withi 1 KM or site is within lsoa 
bus_site<-bus_site0[value<=1000 | exp_lsoa_code==site_lsoa_code]



# # switch the following code on to select those sites with in 500M or site is within lsoa 
# bus_site<-bus_site0[value<=500 | exp_lsoa_code==site_lsoa_code]

# # switch the following code on to select those sites with in 1.5 KM or site is within lsoa 
# bus_site<-bus_site0[value<=1500 | exp_lsoa_code==site_lsoa_code] 

# # switch the following code on to select those sites with in 2 KM or site is within lsoa 
# bus_site<-bus_site0[value<=2000 | exp_lsoa_code==site_lsoa_code] 


# # switch the following code on to select those sites withi 3 KM or site is within lsoa 
# bus_site<-bus_site0[value<=3000 | exp_lsoa_code==site_lsoa_code]



#  number of times LSOA exposed to vaccine bus
bus_site[, num_exp:=.N , by=.(exp_lsoa_code)]
bus_site[, week:=yearweek(date)]
table(bus_site$exp_lsoa_code)
table(bus_site$week)
table(bus_site$num_exp)
table(bus_site$week,bus_site$num_exp)

#only switch the following line of code on for the 1KM threshold
# bus_site[num_exp==6,table(exp_lsoa_code,week)]

#only switch the following line of code on for the 1KM threshold
#unique(bus_site[num_exp==6,.(exp_lsoa_code,week)])

unique(bus_site[exp_lsoa_code=="E01006563" | exp_lsoa_code=="E01006548",
                .(num_exp,week)])
NROW(bus_site[exp_lsoa_code=="E01006563"])
NROW(unique(bus_site[,.(exp_lsoa_code,week,num_exp)]))
NROW(unique(bus_site[,.(exp_lsoa_code,week)]))

bus_site[order(week), visit_num:=1:.N , by=.(exp_lsoa_code)]
table(bus_site$visit_num)
bus_site[, visit_num:=NULL]


bus_site_week<-bus_site[, list(num_exp=min(num_exp)), by=.(exp_lsoa_code,week)]
summary(bus_site_week[,-2])
bus_site_week[exp_lsoa_code=="E01006563" | exp_lsoa_code=="E01006548"]
NROW(unique(bus_site_week))


bus_site_week[, num_exp:=min(num_exp), by=.(exp_lsoa_code)]
summary(bus_site_week[,-2])
bus_site_week[exp_lsoa_code=="E01006563" | exp_lsoa_code=="E01006548"]
NROW(unique(bus_site_week))


bus_site_week[order(week), visit_num:=1:.N , by=.(exp_lsoa_code)]
summary(bus_site_week[,-2])
bus_site_week[,table(visit_num,num_exp)]
bus_site_week[,table(visit_num,week)]

#only switch the following line of code on for the 1KM threshold
#bus_site_week[num_exp==6,table(visit_num,week)]


########################### 5.  Derive look up table  ###########################
#giving the LSOAs exposed, number of times exposed and week of each exposure. 
look_up<-dcast(bus_site_week, exp_lsoa_code+num_exp~visit_num, value.var = "week")
look_up<-clean_names(look_up)
look_up[, vacc_bus:=1]
class(look_up$exp_lsoa_code)
look_up[,exp_lsoa_code:=as.character(exp_lsoa_code)]

#excluding Chester East from the analysis
look_up<-merge(as.data.table(lsoa_cen[grepl(paste(c("Liverpool",
                                                    #"Cheshire East",
                                                    "Cheshire West and Chester","Halton","Knowsley","Sefton",
                                                    "St. Helens","Warrington","Wirral"), collapse="|"),
                                            lsoa_cen$lsoa11nm),])[,c("objectid","geometry"):=NULL], look_up, 
               by.x="lsoa11cd", by.y="exp_lsoa_code", all.x=T)
look_up[is.na(vacc_bus)==T, vacc_bus:=0]
look_up[is.na(num_exp)==T, num_exp:=0]


# impute the pseudo vacc_bus week number for the control group

vacc_prop <- prop.table(table(look_up[vacc_bus==1,.(x1)])) %>% 
  as.data.table()

set.seed(617618)
#detach("package:Tidyverse", unload=TRUE)
look_up[vacc_bus==0,x1:=yearweek(sample(sort(vacc_prop$x1, decreasing = FALSE), 
                                        size=.N, replace=T, 
                                        prob=vacc_prop$N))] 

prop.table(table(look_up$x1, look_up$vacc_bus),2)

check<-look_up[vacc_bus==1 & x1<yearweek("2021 W24")]

table(look_up$vacc_bus)


# use the following line for generating the 1KM threshold lookup table.
write.csv(look_up, normalizePath(file.path(a,"look_up1k.csv")))

# #switch the following line on for excluding popin sites
# write.csv(look_up, normalizePath(file.path(a,"look_up_no_popin.csv")))


# #switch the following line on for 500m threshold
# write.csv(look_up, normalizePath(file.path(a,"look_up5h.csv")))
# #switch the following line on for 500m threshold
# write.csv(look_up, normalizePath(file.path(a,"look_up15h.csv")))
# #switch the following line on for 2KM threshold
# write.csv(look_up, normalizePath(file.path(a,"look_up2k.csv")))
# #switch the following line on for 3KM threshold
# write.csv(look_up, normalizePath(file.path(a,"look_up3k.csv")))

