########################## 0. Load the needed packages #########################
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
library(survminer) # for survival analysis and the '.buildggsurvplot' object
library(survival)
library(MatchIt)
library(Matching)
library(lubridate)
library(biostat3)
library(tidyr)    # for the unnest function
library(margins)
library(jtools)  # for the summ function
library(gtsummary) # for the tbl_regression() function, which takes a regression model object in R and returns a formatted table of regression model results that is publication-ready 
library(aweek)
library(ISOweek)

`%notin%` <- Negate(`%in%`)

# File path to project folders
a<-"/Sharing Folder a"
b<-"/Sharing Folder b"


########################## 1. Load the CIPHA vaccine data  ########################
load(normalizePath(file.path(b,"datasets.RData")))
data7 <- as.data.table(data7)
data5n <- as.data.table(data5n)
rm(data8a)


############## 2. Set up the microdataset   ##################
data5n[, COPD:=as.numeric(grepl("COPD_REG", QOFRegisters )) ]
#CHD, chronic heart disease
data5n[, CHD:=as.numeric(grepl("CHD_REG", QOFRegisters )) ]
#DM, Diabetes
data5n[, DIAB:=as.numeric(grepl("DM_REG", QOFRegisters )) ]
#CKD, Chronic Kidney Disease
data5n[, CKD:=as.numeric(grepl("CKD_REG", QOFRegisters )) ]
#AST, Asthma
data5n[, AST:=as.numeric(grepl("AST_REG", QOFRegisters )) ]
#CAN, Cancer
data5n[, CAN:=as.numeric(grepl("CAN_REG", QOFRegisters )) ]
#OBES, Obesity
data5n[, OBES:=as.numeric(grepl("OBES_REG", QOFRegisters )) ]
#DEP1, Depression
data5n[, DEP:=as.numeric(grepl("DEP1_REG", QOFRegisters )) ]
#STIA, Stroke/TIA
data5n[, STRK:=as.numeric(grepl("STIA_REG", QOFRegisters )) ]


dist <- read.delim(normalizePath(file.path(b,"Distance.txt")), 
                   header = TRUE, sep = "\t", dec = ".")
dist <- as.data.table(dist)
dist <- unique(dist[FK_Patient_Link_ID %in% data5n$FK_Patient_Link_ID,])
colnames(dist)[1]<- "distance_min"
dist[, freq:= .N, by = .(FK_Patient_Link_ID)]
dist[, num:= 1:.N, by = .(FK_Patient_Link_ID)]
dist<-dist[num==1, .(distance_min,distance_km,FK_Patient_Link_ID)]
data5n <- merge(data5n, dist, by="FK_Patient_Link_ID", all.x=T)
rm(dist)


max(data5n$dose1dt, na.rm = T)
test<-data5n[Deceased=="N", .(pid=FK_Patient_Link_ID,age=Patient.Age,sex=Sex,lsoa11=LSOA_Code,ethn=EthnicMainGroup,
                              dose1dt,dose1, dose2,dose2dt, ch_flag=NursingCareHomeFlag,ch_match=MatchType,
                              ch_conf=MatchConfidence,COPD,CHD,DIAB,CKD,AST,CAN,OBES,DEP,STRK,
                              staff=StaffFlag,social_care=SocialCareFlag,
                              carer=cohort_carersDWP,distance_min,distance_km, prevposdt, posdt)]
table(as.numeric(test$distance_km)==0)
table(as.numeric(test$distance_min)==0)
test[,log_distkm := ifelse(as.numeric(test$distance_km)==0,log(as.numeric(distance_km)+0.001),
                           log(as.numeric(distance_km)))]
test[,log_distmin := ifelse(as.numeric(test$distance_min)==0,log(as.numeric(distance_min)+0.001),
                            log(as.numeric(distance_min)))]
test[, posdt:=as.Date(posdt, origin = "1970-01-01")]
test[, dose1:=as.numeric(dose1)]
test[is.na(dose1), dose1:=0]
test[, dose2:=as.numeric(dose1)]
test[is.na(dose2), dose2:=0]
with(test, table(dose1, useNA = "ifany"))
summary(test$age)


test<-test[age>=18]
test[, age_group:=cut(age, breaks=c(18,20,30,40,50,55,60,70,75,80,110), tr.lowest=T, right=F)]
invisible(table(test$age_group))
test[, id:=as.numeric(as.factor(lsoa11))]
table(test$sex)
test[sex=="U", sex:=NA]
test[, female:=ifelse(sex=="F",1,0)]
table(test$female)
with(test, table(sex))


test[, ch_flag:=as.numeric(ch_flag=="Y")]
with(test, prop.table(table(ch_flag)))
with(test, prop.table(table(ethn)))


test[ethn=="NULL", ethn:=NA]
# convert to weekly format
test[, d1_week:=yearweek(dose1dt)]
test[, d2_week:=yearweek(dose2dt)]

table(test$carer)
test[is.na(carer), carer:=0]
table(test$carer)

table(test$social_care)
test[, social_care:=ifelse(social_care=="Y",1,0)]
table(test$social_care)

table(test$staff)
class(test$staff)

max(data5n$dose1dt, na.rm = T)
max(test$dose1dt, na.rm = T)

summary(test$age)

############## 3.   derive cross-sectional LSOA data   ##################

load(normalizePath(file.path(a,"lsoa_risk.RData")))
class(ag_risk_lsoa)
ag_risk_lsoa[,geometry:=NULL]
rural_urban <- fread(normalizePath(file.path(a,"Rural_Urban_Classification_2011_lsoa.csv")), 
                     sep = ",", header= TRUE) %>%
  clean_names(.,"lower_camel")
ag_risk_lsoa <- merge(ag_risk_lsoa,rural_urban[,.(lsoa11Cd,ruc11Cd,ruc11)], 
                      by.x="lsoa11",by.y="lsoa11Cd", all.x=T)


cipha_lsoa_data<-test[, list(pos_past_28=sum(posdt>(max(test$dose1dt, na.rm = T)-28), na.rm=T),
                             mdist_min=mean(as.numeric(distance_min), na.rm = T),
                             mdist_km=mean(as.numeric(distance_km), na.rm = T),
                             prop_female=mean(female, na.rm = T),
                             pop5n_18_29=sum(as.numeric(age>=18 & age<=29), na.rm = T),
                             pop5n_30_39=sum(as.numeric(age>=30 & age<=39), na.rm = T),
                             pop5n_40_49=sum(as.numeric(age>=40 & age<=49), na.rm = T),
                             pop5n_50_59=sum(as.numeric(age>=50 & age<60), na.rm = T),
                             pop5n_60_69=sum(as.numeric(age>=60 & age<70), na.rm = T),
                             pop5n_70_79=sum(as.numeric(age>=70 & age<80), na.rm = T),
                             pop5n_80plus=sum(as.numeric(age>=80 & age<110), na.rm = T),
                             prop5n_18_29=mean(as.numeric(age>=18 & age<=29), na.rm = T),
                             prop5n_30_39=mean(as.numeric(age>=30 & age<=39), na.rm = T),
                             prop5n_40_49=mean(as.numeric(age>=40 & age<=49), na.rm = T),
                             prop5n_50_59=mean(as.numeric(age>=50 & age<60), na.rm = T),
                             prop5n_60_69=mean(as.numeric(age>=60 & age<70), na.rm = T),
                             prop5n_70_79=mean(as.numeric(age>=70 & age<80), na.rm = T),
                             prop5n_80plus=mean(as.numeric(age>=80 & age<110), na.rm = T),
                             pop_asian5n=sum(as.numeric(ethn=="Asian or Asian British"), na.rm = T),
                             pop_black5n=sum(as.numeric(ethn=="Black or Black British"), na.rm = T),
                             pop_mixed5n=sum(as.numeric(ethn=="Mixed"), na.rm = T),
                             pop_other5n=sum(as.numeric(ethn=="Other Ethnic Groups"), na.rm = T),
                             pop_white5n=sum(as.numeric(ethn=="White"), na.rm = T),
                             prop_asian5n=mean(as.numeric(ethn=="Asian or Asian British"), na.rm = T),
                             prop_black5n=mean(as.numeric(ethn=="Black or Black British"), na.rm = T),
                             prop_mixed5n=mean(as.numeric(ethn=="Mixed"), na.rm = T),
                             prop_other5n=mean(as.numeric(ethn=="Other Ethnic Groups"), na.rm = T),
                             prop_white5n=mean(as.numeric(ethn=="White"), na.rm = T),
                             prop_carer=mean(as.numeric(carer), na.rm = T),
                             prop_social_care=mean(as.numeric(social_care), na.rm = T),
                             prop_staff=mean(as.numeric(staff), na.rm = T),
                             cipha_pop=.N
), by=.(lsoa11)]


ag_risk_lsoa<-merge(ag_risk_lsoa,cipha_lsoa_data, by=c("lsoa11") , all.x.y=T)
ag_risk_lsoa[, log_mdist_min:=log(mdist_min)]
ag_risk_lsoa[, log_mdist_km:=log(mdist_km)]
ag_risk_lsoa[,pop_densx:=scale(pop_dens, center = TRUE, scale = TRUE)]
ag_risk_lsoa <- merge(ag_risk_lsoa, unique(data7[,.(mean_age,LSOA_Code)]),
                      by.x="lsoa11",by.y="LSOA_Code",all.x=T)


# merge in with microdata

test<-merge(test,ag_risk_lsoa, by="lsoa11", all.x=T)
test[, prop_1less:=prop_1less*100]
test[, pc_imd:=ntile(imd_score, n=5)]

########################## 4  derive an lsoa balanced panel ##################
# derive balanced panel from first dose vaccine to last dose vaccine in data

min(test$d1_week, na.rm = T)
max(test$d1_week, na.rm = T)
max(test$dose1dt, na.rm = T)
num_weeks<-as.numeric(max(test$d1_week, na.rm = T)-min(test$d1_week, na.rm = T))+1
lsoa_panel <- as.data.table(ag_risk_lsoa[rep(1:.N,each=num_weeks), ])
lsoa_panel[, week:=min(test$d1_week, na.rm = T)+(1:.N)-1, by=.(lsoa11)]

table(lsoa_panel$week)

# should still be a balanced panel 1562 per week at the end of this 
min(lsoa_panel$week, na.rm = T)
max(lsoa_panel$week, na.rm = T)
max(test$d1_week, na.rm = T)

dose1_panel<-test[is.na(d1_week)==F , list(dose1=sum(dose1, na.rm = T)), 
                  by=.(lsoa11,d1_week )]

dose2_panel<-test[is.na(d2_week)==F , list(dose2=sum(dose2, na.rm = T)), 
                  by=.(lsoa11,d2_week )]


lsoa_panel<-merge(lsoa_panel, dose1_panel, by.x=c("week", "lsoa11"), 
                  by.y=c("d1_week", "lsoa11"), all.x=T)
lsoa_panel<-merge(lsoa_panel, dose2_panel, by.x=c("week", "lsoa11"), 
                  by.y=c("d2_week", "lsoa11"), all.x=T)

lsoa_panel[is.na(dose1),dose1:=0 ]
lsoa_panel[is.na(dose2),dose2:=0 ]
table(lsoa_panel$week)
options(max.print=999999)
table(lsoa_panel$week, lsoa_panel$dose1)

########################## 5. Load the vaccine bus site data ########################

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

##exclude the pop-in sites for the sensitivity test
#bus_location <- bus_location[is.na(bus_location$popUpClinic),]

## check distribution of the date of the visits
ggplot(bus_location[bus_location$date<=as.Date("2021-06-28"),], aes(x = date)) +
  geom_histogram(binwidth=.5) +
  labs(title = "Distribution Box Plot by Date",
       x = "Date",
       y = "Frequency")


######################## 6. Load the look up table #############################
look_up<- fread(normalizePath(file.path(a,"look_up1k.csv")), 
                sep = ",", header= TRUE)

# only switch on the following code for sensitivity test of no-popins
# look_up<- fread(normalizePath(file.path(a,"look_up_no_popin.csv")), 
#                 sep = ",", header= TRUE)

# only switch on the following code for sensitivity test of 500M
# look_up<- fread(normalizePath(file.path(a,"look_up5h.csv")), 
#                 sep = ",", header= TRUE)

# only switch on the following code for sensitivity test of 1500M
# look_up<- fread(normalizePath(file.path(a,"look_up15h.csv")), 
#                 sep = ",", header= TRUE)

# only switch on the following code for sensitivity test of 2000M
# look_up<- fread(normalizePath(file.path(a,"look_up2k.csv")), 
#                 sep = ",", header= TRUE)
# only switch on the following code for sensitivity test of 3KM
# look_up<- fread(normalizePath(file.path(a,"look_up3k.csv")), 
#                 sep = ",", header= TRUE)


# only switch on the following code for sensitivity test of spatial spillover effect with distance <= 1KM & !(1KM < distance <= 1.5KM) 
# look_up1.5h<- fread(normalizePath(file.path(a,"look_up15h.csv")), 
#                 sep = ",", header= TRUE)
# unique(look_up[(vacc_bus==1)]$lsoa11cd)[unique(look_up[(vacc_bus==1)]$lsoa11cd) %notin% unique(look_up1.5h[(vacc_bus==1)]$lsoa11cd)]
# table(look_up$vacc_bus)
# look_up[lsoa11cd %in% unique(look_up[(vacc_bus==1)]$lsoa11cd)[unique(look_up[(vacc_bus==1)]$lsoa11cd) %notin% 
#                                                                 unique(look_up1.5h[(vacc_bus==1)]$lsoa11cd)],
#         vacc_bus:=0]
# table(look_up$vacc_bus)



look_up[,V1:=NULL]

look_up[,x1:=yearweek(x1)]
look_up[!is.na(x2),x2:=yearweek(x2)]
look_up[!is.na(x3),x3:=yearweek(x3)]
look_up[!is.na(x4),x4:=yearweek(x4)]

# need to switch the next line off for distance threshold 500M
look_up[!is.na(x5),x5:=yearweek(x5)]

################### 7. merge look up with lsoa panel ###################

min(as.Date(look_up[vacc_bus==1]$week), na.rm=T)
yearweek("2021 W15")
as.Date(yearweek("2021 W15"))

yearweek("2021 W26")
as.Date(yearweek("2021 W26"))

# the earliest vaccine site is 12th April 2021 - "2021 W15"
# latest date is 28th June 2021 - ("2021 W26")

# for all intervention LSOAs we would need data from w 8 to W 29
lsoa_panel<-lsoa_panel[is.na(dose1)==F]
# we have it from 2019 to to week 35

#  so all intervention lsoas should be included 
length(unique(look_up[vacc_bus==1]$lsoa11cd))
length(unique(look_up[vacc_bus==0]$lsoa11cd))
# that's what we have - 216 and 1112

lsoa_panel<-merge(lsoa_panel, look_up, by.x="lsoa11", by.y="lsoa11cd", all.x=T)
summary(lsoa_panel$vacc_bus)
lsoa_panel[, vacc_rate:=dose1*100/cipha_pop]
summary(lsoa_panel$vacc_rate)


lsoa_panel[order(week), cum_dose1:=cumsum(dose1), by=.(lsoa11)]
lsoa_panel[, cum_vacc_rate:=cum_dose1*100/cipha_pop]
summary(lsoa_panel$cum_vacc_rate)

lsoa_panel[order(week),cum_dose1_growth:=(cum_vacc_rate - shift(cum_vacc_rate,1))/shift(cum_vacc_rate,1),
           by=.(lsoa11)]
summary(lsoa_panel$cum_dose1_growth)


hist(lsoa_panel$cum_vacc_rate)


lsoa_panel[, week_from_vaccbus:=week-x1]
table(lsoa_panel$week_from_vaccbus)
# lsoa_panel<-lsoa_panel[week_from_vaccbus>-20]
table(lsoa_panel$week_from_vaccbus)


#excluding Chester East
fig1<-lsoa_panel[!is.na(vacc_bus) &week_from_vaccbus>-20 , list(dose1=sum(dose1), cipha_pop=sum(cipha_pop)), 
                 by=.(week_from_vaccbus, vacc_bus)]
fig1[, vacc_rate:=dose1*100/cipha_pop]
dev.off()

## visualise the vacc_bus effect
ggplot(data=fig1,aes(x=week_from_vaccbus,y=vacc_rate,group=as.factor(vacc_bus), 
                     color=as.factor(vacc_bus)))+
  geom_line(size=2,alpha=0.5) +
  geom_vline(xintercept=0, colour="red", linetype = "dashed") +
  scale_color_manual(values = c("#74D055FF", "#50127CFF"),name="Mobile vaccine units") + 
  labs(x ="vacc_bus week number",          
       y="Weekly vaccination rate") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face="bold"))



ggplot(data=fig1,aes(x=week_from_vaccbus,y=dose1,group=as.factor(vacc_bus), 
                     color=as.factor(vacc_bus)))+
  geom_line(size=2,alpha=0.5) +
  geom_vline(xintercept=0, colour="red", linetype = "dashed") +
  scale_color_manual(values = c("#74D055FF", "#50127CFF"),name="Mobile vaccine units") + 
  labs(x ="vacc_bus week number",          
       y="Weekly number of dosages (dose 1)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face="bold"))


ggplot(data=fig1[vacc_bus==1],aes(x=week_from_vaccbus,y=dose1,group=as.factor(vacc_bus), 
                                  color=as.factor(vacc_bus)))+
  geom_line(size=2,alpha=0.5) +
  geom_vline(xintercept=0, colour="red", linetype = "dashed") +
  scale_color_manual(values = "#50127CFF",name="Mobile vaccine units") + 
  labs(x ="vacc_bus week number",          
       y="Weekly number of dosages (dose 1)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face="bold"))



########### 8. Microsynth analysis  of the impact of the vaccine bus ###########
table(lsoa_panel[vacc_bus==1]$week_from_vaccbus)

test2<-lsoa_panel[week_from_vaccbus>-8 & week_from_vaccbus<4]


table(test2$week_from_vaccbus, test2$vacc_bus)
#  those with at least 10 weeks of data 
test2[, weeks:=.N, by=.(lsoa11)]
table(test2$weeks)
test2<-test2[weeks==max(test2$weeks)]
table(test2[vacc_bus==1]$week_from_vaccbus)

table(test2$week_from_vaccbus, test2$vacc_bus)

test2[, time:=as.numeric(as.factor(week_from_vaccbus))]
table(test2$time)
table(test2$week_from_vaccbus, test2$time)
test2[, id:=as.factor(lsoa11)]
test2[week_from_vaccbus==0, table(week)]

NROW(unique(data5n[MSOA=="Cheshire East",LSOA_Code]))
summary(test2$vacc_bus)
NROW(unique(test2[week_from_vaccbus==0, lsoa11]))

##number of LSOAs excluded from the analysis
NROW(unique(data5n[,LSOA_Code])) - 
  NROW(unique(data5n[MSOA=="Cheshire East",LSOA_Code])) - 
  NROW(unique(test2[week_from_vaccbus==0, lsoa11]))


fig2<-test2[, list(dose1=sum(dose1), cipha_pop=sum(cipha_pop)), 
            by=.(week_from_vaccbus, vacc_bus)]

fig2[, vacc_rate:=dose1*100/cipha_pop]


## visualise the vacc_bus effect
ggplot(data=fig2,aes(x=week_from_vaccbus,y=vacc_rate,group=as.factor(vacc_bus), 
                     color=as.factor(vacc_bus)))+
  geom_line(size=2,alpha=0.5) +
  geom_vline(xintercept=0, colour="red", linetype = "dashed") +
  scale_color_manual(values = c("#74D055FF", "#50127CFF"),name="Mobile vaccine units") + 
  scale_x_continuous(breaks=seq(-7,1,2)) +
  labs(x ="vacc_bus week number",          
       y="Weekly vaccination rate") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face="bold"))

table(test2[vacc_bus==1]$ladnm)
table(test2$time)
table(test2$week_from_vaccbus, test2$time)
class(test2$mean_age)
test2$mean_age <- as.numeric(test2$mean_age)

#  I think we have cumulative doses at the start of study period the have weekly dose as the outcome - cumulative doses will be a product of this
test2[, cum_dose_b:=cum_dose1[time=1], by=.(lsoa11)]
summary(test2[, c("cipha_pop","prop_female","pop_dens","mean_age","prop_asian5n", 
                  "prop_black5n", "prop_mixed5n","ch_prop","imd_score","gpp_dist",
                  "prop_1less","mdist_min","med_gb_use", "num_exp","prop_carer",
                  "prop_social_care","prop_staff","vacc_bus","id","time","cum_dose1", "cum_dose_b")])


cov.var <- c("cipha_pop","prop_female","pop_dens","mean_age","prop_asian5n", 
             "prop_black5n", "prop_mixed5n","imd_score","cum_dose_b", 
             "prop_1less","mdist_min")

summary(test2[,c("cipha_pop","prop_female","pop_dens","mean_age","prop_asian5n", 
                 "prop_black5n", "prop_mixed5n","imd_score","cum_dose_b", 
                 "prop_1less","mdist_min")])


end.post<-max(test2$time)

table(test2[time==8]$vacc_bus)
table(as.Date(test2[time==7]$week))
get_date(week = 25, year = 2021)
set.seed(617618)
sea2 <- microsynth(as.data.frame(test2), 
                   idvar="id", timevar="time", intvar="vacc_bus", 
                   start.pre=1, end.pre=8, end.post=end.post, 
                   match.out=c("dose1"), match.covar=cov.var, 
                   result.var=c("dose1"), 
                   test="twosided",
                   use.backup = TRUE,
                   confidence = 0.95,
                   perm=250,
                   omnibus.var = NULL,
                   jack=F, 
                   n.cores = 1)
sea2
par(mar=c(2,2,2,2))

table(test2$ladnm)

plot_microsynth(sea2,height = 11,sep = F,main.tc="",main.diff = "",
                ylab.tc="Number of first doses",
                width = 8.5,
                legend.spot = "topleft", plot.var = "dose1")

######### 9. Results table #######
table1 <- as.data.table(sea2$Results$`11`)
write.csv(table1, normalizePath(file.path(a,"main_result_lsoa_1k.csv")))

# #only switch the next line on for distance threshold 500M
# write.csv(table1, normalizePath(file.path(a,"main_result_lsoa_5h.csv")))

# #only switch the next line on for distance threshold 1.5KM
# write.csv(table1, normalizePath(file.path(a,"main_result_lsoa_15h.csv")))

# #only switch the next line on for distance threshold 2KM
# write.csv(table1, normalizePath(file.path(a,"main_result_lsoa_2K.csv")))

# #only switch the next line on for distance threshold 3KM
# write.csv(table1, normalizePath(file.path(a,"main_result_lsoa_3K.csv")))

# #only switch the next line on for testing spatial spillover effect with distance threshold [1,1.5]KM
# write.csv(table1, normalizePath(file.path(a,"main_result_spill_1_1.5K.csv")))



#  add weights into main dataset. 
weights2<-cbind(ID=row.names(sea2$w$Weights),as.data.frame(sea2$w$Weights),
                row.names=NULL)
weights2<-cbind(weights2,interv2=as.data.frame(sea2$w$Intervention), 
                row.names=NULL)
colnames(weights2)[1]<-"id"

summary(weights2$Main)

test2$Main<-NULL
test2<-as.data.table(merge(test2,weights2[,1:2], by="id", 
                           all.x=T ))


num1<-sum(test2[time>8 & vacc_bus==1]$dose1*test2[time>8 & vacc_bus==1]$Main)

num2<-sum(test2[time>8 & vacc_bus==0]$dose1*test2[time>8 & vacc_bus==0]$Main)

num1/num2

num1-num2

num1-(num1/1.25)

num1_ex_war<-sum(test2[time>8 & vacc_bus==1 &ladnm!="Warrington"]$dose1*test2[time>8 & vacc_bus==1&ladnm!="Warrington"]$Main)

(num1_ex_war-(num1_ex_war/1.25))/3800

(num2*1.25)-num2

num1<-sum(test2[time>8 & vacc_bus==1]$dose1*test2[time>8 & vacc_bus==1]$Main)

table(test2$ladnm, test2$ladnm)

svyfull<-svydesign(ids=~0, weights=~Main, data=test2)



prop2<-cbind(svyby(~dose1, denominator=~cipha_pop, by=~vacc_bus+week_from_vaccbus,
                   design=svyfull,FUN=svyratio),
             confint(svyby(~dose1, denominator=~cipha_pop, 
                           by=~vacc_bus+week_from_vaccbus ,design=svyfull,FUN=svyratio)))


names(prop2)<-c("vacc_bus", "weeks","uptake", "se_rate","lcl", "ucl")


ggplot(data=prop2, aes(x = weeks, y = uptake*100, color=as.factor(vacc_bus)))+ 
  geom_line(aes(lty=as.factor(vacc_bus)), size=1)+
  geom_vline(xintercept=0, colour="red", linetype = "dashed") +
  geom_ribbon(aes(fill=as.factor(vacc_bus), 
                  ymin=lcl*100, ymax=ucl*100),size=0, alpha=0.3)+
  scale_colour_manual(values = c("#440154FF","#FDE725FF"), 
                      labels=c("Synthetic Control", "Intervention")) +
  scale_linetype_manual(values = c("solid","dashed"), 
                        labels=c("Synthetic Control", "Intervention")) +
  scale_fill_manual(values = c("#440154FF","#FDE725FF"), 
                    labels=c("Synthetic Control", "Intervention")) +
  scale_x_continuous(breaks=seq(-7,3,1)) +
  xlab("Weeks since vaccine bus") + 
  ylab("Weekly vaccination rate %")+
  theme_classic()+
  theme(text = element_text(size=16),legend.position="bottom", 
        legend.box = "horizontal", legend.text = element_text(size=24),
        legend.title = element_blank(),legend.key.size = unit(1, "cm"), 
        axis.text.x = element_text(angle = 45,vjust=0.5))+ 
  theme(strip.background = element_blank())

ggsave(normalizePath(file.path(a,"lsoa_ci.png")), width = 16, height = 14,device="png",
       units = "in", dpi=700)

ggsave(normalizePath(file.path(a,"lsoa_ci.jpeg")), width = 16, height = 14,device="jpeg",
       units = "in", dpi=400)

ggsave(normalizePath(file.path(a,"lsoa_ci.svg")), width = 16, height = 14,device="svg",
       units = "in", dpi=700)


######################## 10. Map out the study areas of sea2  ######################
# 
# load LAD boundaries that are downloaded from https://geoportal.statistics.gov.uk/datasets/local-authority-districts-december-2017-boundaries-gb-bfc?geometry=-35.772%2C51.103%2C30.981%2C59.783
lad <- read_sf(dsn = normalizePath(file.path(b,"lad.shp")))
#lad <- lad[lad$LAD17NM %in% data5n$MSOA,]
lad <- lad[substr(lad$LAD17CD, 1,1)=="E",]
lad <- st_transform(lad, 4326)


# map out the bus locations
tmap_mode("view")
tm_shape(lad) +
  tm_borders() +
  tm_text("LAD17NM", size="Shape__Are") +
  tm_shape(bus_location) +
  tm_bubbles(col = "postcd",
             size = 0.1,
             palette = "-viridis")

# load LSOA boundaries that are downloaded from https://geoportal.statistics.gov.uk/datasets/lower-layer-super-output-areas-december-2011-boundaries-full-clipped-bfc-ew-v3?geometry=-35.549%2C48.013%2C31.203%2C57.298
lsoas <- read_sf(dsn =normalizePath(file.path(b, "lsoa.shp")))
lsoas <- lsoas[substr(lsoas$LSOA11CD, 1,1)=="E",]
lsoas <- st_transform(lsoas, 4326)

table(test2[week_from_vaccbus==0]$vacc_bus)
sp_lsoa <- merge(lsoas,test2[week_from_vaccbus==0,.(week_from_vaccbus,vacc_bus,Main,lsoa11)],
                 by.x="LSOA11CD", by.y="lsoa11", all.x=T)
table(test2$week_from_vaccbus,test2$vacc_bus)
table(sp_lsoa$vacc_bus)
class(sp_lsoa$vacc_bus)

str(bus_location)
class(bus_location$lsoa_code)
bus_location <- unnest(bus_location, lsoa_code)
str(bus_location)
class(bus_location)
bus_location <- merge(bus_location,test2[week_from_vaccbus==0,.(week_from_vaccbus,vacc_bus,Main,lsoa11)],
                      by.x="lsoa_code", by.y="lsoa11", all.x=T)
table(test2$week_from_vaccbus,test2$vacc_bus)
table(bus_location$vacc_bus)
summary(bus_location[bus_location$date<="2021-06-30",]$vacc_bus)
unique(bus_location[bus_location$date<="2021-06-30" & is.na(bus_location$vacc_bus),
                    "postcd"])
NROW(unique(bus_location[bus_location$date<="2021-06-30" & is.na(bus_location$vacc_bus),
                         "postcd"]))
NROW(bus_location[bus_location$date<="2021-06-30" & is.na(bus_location$vacc_bus),
                  "postcd"])


###plotting the C&M areas### 
c_m<-c("Sefton","Knowsley","St. Helens", "Warrington","Halton","Cheshire East",
       "Cheshire West and Chester","Wirral","Liverpool")

ggplot() +
  geom_sf(data = sp_lsoa,aes(fill = as.factor(vacc_bus)),
          colour = NA,
          alpha=1,
          lwd = 0) +
  scale_fill_manual(values = c("#440154FF","#FDE725FF","white"), na.value = "white",
                    guide = "none") +
  geom_sf(data = lad[lad$LAD17NM %notin% c_m,],aes(fill = "grey"), colour = "grey",
          na.rm = F,alpha=0,size = 0.5) +
  geom_sf(data = st_union(lad[lad$LAD17NM %in% c_m,]),
          aes(fill = NA), na.rm = F,alpha=0,size = 1) +
  guides(color = FALSE, size = FALSE,scale = "none") +
  theme(legend.position = "none") +
  theme_void()

ggsave(normalizePath(file.path(a,"CM.png")), width = 16, height = 14,device="png",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"CM.svg")), width = 16, height = 14,device="svg",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"CM.jpeg")), width = 16, height = 14,device="jpeg",
       units = "in", dpi=400)


###plotting the study areas### 
table(sp_lsoa$vacc_bus)
table(bus_location$vacc_bus)

table(bus_location$popUpClinic)
bus_location[is.na(bus_location$popUpClinic),]$popUpClinic <- "No"
table(bus_location$popUpClinic)
table(bus_location[bus_location$date<=as.Date("2021-06-30"),]$popUpClinic)#
table(as.factor(bus_location[bus_location$date<=as.Date("2021-06-30"),]$popUpClinic))
levels(as.factor(bus_location[bus_location$date<=as.Date("2021-06-30"),]$popUpClinic))

ggplot() +
  geom_sf(data = sp_lsoa[sp_lsoa$LSOA11CD%in% unique(data5n$LSOA_Code),],
          aes(fill = as.factor(vacc_bus)),
          colour = NA,
          lwd = 0) +
    scale_fill_manual(values = c("#440154FF","#FDE725FF","white"), 
                      labels=c("Non-intervention", "Intervention"),
                      na.value = "white",
                      na.translate = F,
                      guide_legend(title = "")) +
   geom_sf(data = bus_location[bus_location$date<=as.Date("2021-06-30"),],
           aes(fill = NA,colour = as.factor(popUpClinic)),
           fill = NA,
           size = 1) +
   scale_colour_manual(values = c("No" = "red", "Yes" = "cyan"),
                       labels= c("Mobile sites", "Static sites"),
                       na.value = "white",
                       guide_legend(),
                       guide = guide_legend(title = "",
                                            override.aes = list(linetype = c("blank", "blank"), 
                                                                size = 2))) +
  geom_sf(data = lad[lad$LAD17NM%in%c_m,],aes(fill = NA), na.rm = F,alpha=0,size = 0.1) +
  guides(fill = guide_legend(title = "")) +
  theme_void() +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=20)) +
        # legend.position = c(0.8, 0.05))+
  annotation_scale(location = "br", width_hint = 0.25, text_cex=1) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)

ggsave(normalizePath(file.path(a,"sea2_areas_1k_clean_updated.png")), width = 16, height = 14,device="png",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"sea2_areas_1k_clean_updated.svg")), width = 16, height = 14,device="svg",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"sea2_areas_1k_clean_updated.jpeg")), width = 16, height = 14,device="jpeg",
       units = "in", dpi=400)


###plotting out the weights### 
summary(sp_lsoa$Main)

ggplot() +
  geom_sf(data = sp_lsoa[sp_lsoa$Main!=0&sp_lsoa$vacc_bus==0,],
          aes(fill = as.numeric(Main)),
          colour = NA,
          lwd = 0) +
  scale_fill_viridis(
    trans = "sqrt", 
    direction = -1,
    option = "C",
    breaks = c(0,0.1,2.0,4.0,6.0,8.0), 
    labels = c("0.0","0.1","2.0","4.0", "6.0", "8.0"),
    na.value = "white",
    name="Synthetic control weights") +
  geom_sf(data = sp_lsoa[sp_lsoa$vacc_bus==1,],
          aes(fill = vacc_bus),
          colour = "black",
          fill = "black",
          lwd = 0) +
  geom_sf(data = bus_location[bus_location$date<=as.Date("2021-06-30"),],
          aes(fill = NA,colour = as.factor(popUpClinic)),
          fill = NA,
          size = 1) +
  scale_colour_manual(values = c("No" = "red", "Yes" = "cyan"),
                      labels= c("Mobile sites", "Static sites"),
                      na.value = "white",
                      guide_legend(),
                      guide = guide_legend(title = "Vaccination sites",
                                           override.aes = list(linetype = c("blank", "blank"), 
                                                               size = 2))) +
  geom_sf(data = lad[lad$LAD17NM%in%c_m,],
          aes(fill = NA), na.rm = F,alpha=0,size = 0.1) +
  theme_void() +
  theme(legend.text=element_text(size=16),
        legend.title=element_text(size=20),
        legend.position = c(0.9, 0.1))+
  annotation_scale(location = "bl", width_hint = 0.25, text_cex=1) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         pad_x = unit(0.1, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)

ggsave(normalizePath(file.path(a,"sea2_areas_1k_clean_updated_wts.png")), width = 16, height = 14,device="png",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"sea2_areas_1k_clean_updated_wts.svg")), width = 16, height = 14,device="svg",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"sea2_areas_1k_clean_updated_wts.jpeg")), width = 16, height = 14,device="jpeg",
       units = "in", dpi=400)

#
######################## 11. Indidvual level analysis  #########################
# construct survey data
# people unvaccinated a week before the vacc_bus date, including:
# week_from_vaccbus>=-1, people hadn't vaccinated before the vacc_bus date but received the vaccine by the end of the study period;
# is.na(week_from_vaccbus), people hadn't vaccinated either before the vacc_bus date or by the end of the study period;

# merge in look up

look_up2<-test2[week_from_vaccbus<0, 
                .(lsoa11, cum_vacc_rate, vacc_bus, int_week=x1,cum_dose1_growth, 
                  week_from_vaccbus, Main)]
test2[, panel:=as.numeric(as.factor(lsoa11))]


look_up2<-dcast(look_up2, lsoa11+vacc_bus+int_week+Main~week_from_vaccbus, 
                value.var = "cum_vacc_rate")


look_up2<-clean_names(look_up2)

prop.table(table(test$ethn, useNA = "ifany"))
surv_data<-merge(test, 
                 look_up2, by="lsoa11", all.x.y=T)

full_sample <- surv_data
# drop if vaccinated before bus first arrived 

surv_data<-surv_data[d1_week>int_week | is.na(d1_week)==T]

NROW(full_sample)-NROW(surv_data)


# recode dose1 for people who got vaccinated 3 weeks after the vacc_bus as 0.

surv_data[, week_from_vaccbus:=d1_week-int_week]
table(surv_data$week_from_vaccbus)

surv_data[, outcome2:=dose1]
surv_data[week_from_vaccbus>3 & dose1==1,outcome2:=0]


table(surv_data$ethn)
surv_data[, ethn_f:=as.factor(ethn)]
table(surv_data$ethn_f)
surv_data[,ethn_f:=relevel(ethn_f, ref=5)]
surv_data[,pop_densx:=scale(pop_dens, center = TRUE, scale = TRUE)]


surv_data[, pc_imd:=ntile(imd_score, n=5)]

table(surv_data$pc_imd)

class(surv_data$COPD)
table(surv_data$COPD)
surv_data[,chronic:=sum(COPD+CHD+DIAB+CKD+AST+CAN+OBES+DEP+STRK),by=pid]
table(surv_data$chronic)

table(surv_data$ch_flag,surv_data$ch_conf)
surv_data[,ch_weighted:= ifelse(as.numeric(ch_flag)==1 & as.numeric(ch_conf)<=4,1,0)]
table(surv_data$ch_weighted)
class(surv_data$ch_weighted)

#summary(surv_data$cum_vacc_rate)
#summary(surv_data$cum_dose1_growth)


cols <- c('sex','ethn_f','pc_imd','ch_weighted','ruc11','carer','social_care','staff')
surv_data[,(cols):=lapply(.SD, as.factor),.SDcols=cols]

table(surv_data$vacc_bus)
check<-surv_data[, list(outcome2=sum(outcome2), pop=.N), 
                 by=.(lsoa11,cipha_pop,main, vacc_bus)]

check[, panel:=as.numeric(as.factor(lsoa11))]



# just remove those with zero weights to speed up - these are not used in the calculation anyway
summary(surv_data$main)
NROW(surv_data[main==0])

prop.table(table(complete.cases(surv_data[,.(vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1)])))
summary(surv_data[,.(vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1)])
prop.table(table(complete.cases(surv_data$sex)))
prop.table(table(complete.cases(surv_data$ethn_f)))

surv0 <- surv_data[,.(vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1)]


surv_data<-surv_data[main>0]


table(complete.cases(surv_data[,.(vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1)]))
prop.table(table(complete.cases(surv_data[,.(vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1)])))
summary(surv_data[,.(vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1)])
prop.table(table(complete.cases(surv_data$sex)))
prop.table(table(complete.cases(surv_data$ethn_f)))


table(surv_data$ch_weighted, useNA = "ifany")
table(surv_data$carer, useNA = "ifany")
table(surv_data$social_care, useNA = "ifany")
prop.table(table(surv_data$ethn, useNA = "ifany"))
prop.table(table(surv_data$ethn, useNA = "ifany"))
prop.table(table(surv_data$ethn, useNA = "ifany"))
prop.table(table(is.na(surv_data$ethn)==T))
prop.table(table(is.na(surv_data$age)==T))
prop.table(table(is.na(surv_data$sex)==T))
prop.table(table(is.na(surv_data$carer)==T))
prop.table(table(is.na(surv_data$imd_score)==T))

options(scipen=999)

model1<-svyglm(outcome2 ~ vacc_bus+age+sex+ethn_f+imd_score+chronic+carer+
                 social_care+mdist_min, 
               family=poisson(link=log), 
               design=svydesign(ids=~pid, weights =~main, data=surv_data))
summary(model1)
lincom(model1, "vacc_bus", eform=T)

rr1<-as.data.table(tidy(model1,conf.int = T))
rr1[, estimate:=exp(estimate)]
rr1[, conf.low:=exp(conf.low)]
rr1[, conf.high:=exp(conf.high)]
rr1<-rr1[, .(term, estimate,conf.low, conf.high, p.value)]
rr1[,result:="with confounder"]


#only switch it on for the sensitivity test of excluding the popin sites
#write.csv(rr1, normalizePath(file.path(a,"main_result_no_popin.csv")))


# by imd - cut it to three groups 
surv_data[, pc_imd:=as.factor(ntile(imd_score, n=3))]


###################### 12. All interactions simultaneously #####################
########## using main weights
surv_data[,age18_64:=ifelse(age>=18 & age<65,1,0)]
summary(surv_data$age18_64)
table(surv_data$age18_64)
table(surv_data$pc_imd)
surv_data[, pc_imd:=as.factor(ntile(imd_score, n=3))]
table(surv_data$pc_imd)
surv_data[, agegrp:=cut(age, breaks=c(18,30,65,200))]
table(surv_data$agegrp)


surv_data[,agesq2:=age*age]
summary(surv_data[,.(age,agesq2)])
surv_data[,imdsc_sq2:=imd_score*imd_score]
summary(surv_data[,.(imd_score,imdsc_sq2)])
summary(surv_data[pc_imd==3,.(imd_score)])
all_inter<-svyglm(outcome2 ~ pc_imd*vacc_bus+ethn_f*vacc_bus+agegrp*vacc_bus+
                    sex+chronic+carer+
                    social_care+mdist_min, 
                  family=poisson(link=log), 
                  design=svydesign(ids=~pid, weights =~main, data=surv_data))
summary(all_inter)

# The tbl_regression() function takes a regression model object in R and returns a formatted table of regression model results 
# that is publication-ready. It is a simple way to summarize and present your analysis results using R! Like tbl_summary(), 
# tbl_regression() creates highly customizable analytic tables with sensible defaults.
tbl_regression(all_inter, exponentiate = TRUE)

write.csv(summ(all_inter, model.info = FALSE, model.fit = FALSE,digits = 3)$coeftable, 
          normalizePath(file.path(a,"all_inter_reg.csv")))


lincom(all_inter, c("vacc_bus+pc_imd2:vacc_bus"),eform=T)
lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British"),eform=T)
lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British"),eform=T)

res_inter<-rbindlist(list(as.data.table(lincom(all_inter, c("vacc_bus"), eform=T)),
                          as.data.table(lincom(all_inter, c("vacc_bus+pc_imd2:vacc_bus"),eform=T)),
                          as.data.table(lincom(all_inter, c("vacc_bus+pc_imd3:vacc_bus"),eform=T)),
                          as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British"),eform=T)),
                          as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British"),eform=T)),
                          as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed"),eform=T)),
                          as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups"),eform=T)),
                          as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:agegrp(30,65]"),eform=T)),
                          as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:agegrp(65,200]"),eform=T))))


names(res_inter)<-c('estimate','conf.low', 'conf.high','chi', 'p.value')
res_inter[2, p2:=lincom(all_inter, c("pc_imd2:vacc_bus"))[5]]
res_inter[3, p2:=lincom(all_inter, c("pc_imd3:vacc_bus"))[5]]
res_inter[4, p2:=lincom(all_inter, c("vacc_bus:ethn_fAsian or Asian British"))[5]]
res_inter[5, p2:=lincom(all_inter, c("vacc_bus:ethn_fBlack or Black British"))[5]]
res_inter[6, p2:=lincom(all_inter, c("vacc_bus:ethn_fMixed"))[5]]
res_inter[7, p2:=lincom(all_inter, c("vacc_bus:ethn_fOther Ethnic Groups"))[5]]
res_inter[8, p2:=lincom(all_inter, c("vacc_bus:agegrp(30,65]"))[5]]
res_inter[9, p2:=lincom(all_inter, c("vacc_bus:agegrp(65,200]"))[5]]

#res_inter[1, result:="individual - IMD1"]
res_inter[2, result:="individual - IMD2"]
res_inter[3, result:="individual - IMD3"]
#res_inter[4, result:="individual - White"]
res_inter[4, result:="individual - Asian"]
res_inter[5, result:="individual - Black"]
res_inter[6, result:="individual - Mixed"]
res_inter[7, result:="individual - Other Ethnic Groups"]
#res_inter[8, result:="individual - aged 65 and above"]
res_inter[8, result:="individual - aged between 30 and 65"]
res_inter[9, result:="individual - aged between 66 and 100"]


write.csv(res_inter, normalizePath(file.path(a,"appendix_inter.csv")))

#only switch the following line on for the sensity test of excluding the two popin sites
#write.csv(res_inter, normalizePath(file.path(a,"papers/vaccine_bus/table4_updated_no_popin.csv")))



####################### 13. Heat map for the result table ######################
res_heat<-rbindlist(list(
  # white
  as.data.table(lincom(all_inter, c("vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:agegrp(30,65]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:agegrp(30,65]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:agegrp(30,65]+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:agegrp(65,200]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:agegrp(65,200]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:agegrp(65,200]+pc_imd3:vacc_bus"), eform=T)),
  # Asian or Asian British
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British+vacc_bus:agegrp(30,65]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British+vacc_bus:agegrp(30,65]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British+vacc_bus:agegrp(30,65]+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British+vacc_bus:agegrp(65,200]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British+vacc_bus:agegrp(65,200]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fAsian or Asian British+vacc_bus:agegrp(65,200]+pc_imd3:vacc_bus"), eform=T)),
  # Black or Black British
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British+vacc_bus:agegrp(30,65]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British+vacc_bus:agegrp(30,65]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British+vacc_bus:agegrp(30,65]+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British+vacc_bus:agegrp(65,200]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British+vacc_bus:agegrp(65,200]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fBlack or Black British+vacc_bus:agegrp(65,200]+pc_imd3:vacc_bus"), eform=T)),
  # Mixed
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed+vacc_bus:agegrp(30,65]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed+vacc_bus:agegrp(30,65]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed+vacc_bus:agegrp(30,65]+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed+vacc_bus:agegrp(65,200]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed+vacc_bus:agegrp(65,200]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fMixed+vacc_bus:agegrp(65,200]+pc_imd3:vacc_bus"), eform=T)),
  # Other Ethnic Groups
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups+vacc_bus:agegrp(30,65]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups+vacc_bus:agegrp(30,65]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups+vacc_bus:agegrp(30,65]+pc_imd3:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups+vacc_bus:agegrp(65,200]"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups+vacc_bus:agegrp(65,200]+pc_imd2:vacc_bus"), eform=T)),
  as.data.table(lincom(all_inter, c("vacc_bus+vacc_bus:ethn_fOther Ethnic Groups+vacc_bus:agegrp(65,200]+pc_imd3:vacc_bus"), eform=T))
))


names(res_heat)<-c('estimate','conf.low', 'conf.high','chi', 'p.value')
res_heat[1:9,ethn_f:="White or White British"]
res_heat[10:18,ethn_f:="Asian or Asian British"]
res_heat[19:27,ethn_f:="Black or Black British"]
res_heat[28:36,ethn_f:="Mixed"]
res_heat[37:45,ethn_f:="Other Ethnic Groups"]
res_heat[,agegrp:=rep(c("(18,30]","(30,65]","(65,110]"),each=3,times=5)]
res_heat[,pc_imd:=rep(1:3,times=15)]
res_heat[,vacc_bus:=1]


res_heat[,estp:=as.character(format(round(estimate, 2), nsmall = 2))]


ggplot(data = res_heat, mapping = aes(x = pc_imd,
                                      y = agegrp,
                                      fill = estimate)) +
  scale_fill_viridis(name="Estimate",
                     begin = 0.1, end = 1,
                     na.value="white") +
  geom_tile(aes(fill = estimate)) +
  geom_text(aes(label = estp), fontface = "bold", size = 12,hjust="left") +
  facet_grid(ethn_f~vacc_bus)+
  theme_classic() +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"))+
  xlab(label = "Tercile of deprivation (from least to most deprived)") +
  ylab(label = "Age group")
#
#
ggsave(normalizePath(file.path(a,"res_heatmap.png")), width = 16, height = 14,device="png",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"res_heatmap.svg")), width = 16, height = 14,device="svg",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"res_heatmap.jpeg")), width = 16, height = 14,device="jpeg",
       units = "in", dpi=400)


##################################  14. Summary statistics for the main analysis ##########################
summary(test2$prop_1less)
summary(test2$prop_female)


table5a <- test2[week_from_vaccbus<0,list(cipha_pop=sum(cipha_pop),
                                          prop_female=weighted.mean(prop_female,w=cipha_pop)*100,
                                          pop_dens=weighted.mean(pop_dens,w=cipha_pop),
                                          mean_age=weighted.mean(mean_age,w=cipha_pop),
                                          prop_asian5n=weighted.mean(prop_asian5n,w=cipha_pop)*100,
                                          prop_black5n=weighted.mean(prop_black5n,w=cipha_pop)*100,
                                          prop_mixed5n=weighted.mean(prop_mixed5n,w=cipha_pop)*100,
                                          imd_score=weighted.mean(imd_score,w=cipha_pop),
                                          prop_1less=weighted.mean(prop_1less,w=cipha_pop)*100,
                                          mdist_min=weighted.mean(mdist_min),
                                          cum_vacc_rate=weighted.mean(cum_vacc_rate[week_from_vaccbus==-1], w=cipha_pop[week_from_vaccbus==-1]),
                                          #   dose1=weighted.mean(dose1),
                                          vacc_rate=weighted.mean(vacc_rate,w=cipha_pop),
                                          x_numlsoa=.N),by=.(vacc_bus)]

table5a[, cipha_pop:= cipha_pop/7]
table5a[, x_numlsoa:=x_numlsoa/7]
table5a<-as.data.table(t(table5a), keep.rownames=T)


cols <- c('V1','V2')
# table5a[,(cols):=lapply(.SD, round),.SDcols=cols]
table5a[,(cols) := round(.SD,2), .SDcols=cols]
## this is to report more decimal places as Dr Philip Britteon requested



write.csv(table5a,
          file= normalizePath(file.path(a,"var_summary.csv")))

##################################  15. Check missing ##########################
table(test$ethn)
test[, ethn_f:=as.factor(ethn)]
table(test$ethn_f)
test[,ethn_f:=relevel(ethn_f, ref=5)]

test[,chronic:=sum(COPD+CHD+DIAB+CKD+AST+CAN+OBES+DEP+STRK),by=pid]
table(test$chronic)


test0<-merge(test, look_up2, by="lsoa11", all.x.y=T)


test0 <- test0[, .(pid,vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1,dose1dt)]
summary(test0$age)
summary(test0)

table(test0$sex)
table(data5n$Sex)
table(complete.cases(data5n$Sex))

NROW(data5n[Deceased=="Y",])
NROW(data5n[Patient.Age<18,])

prop.table(table(complete.cases(test0[,.(vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1)])))


full_sample[, ethn_f:=as.factor(ethn)]
full_sample[,ethn_f:=relevel(ethn_f, ref=5)]
full_sample[,chronic:=sum(COPD+CHD+DIAB+CKD+AST+CAN+OBES+DEP+STRK),by=pid]

full_sample<-full_sample[,.(pid,vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1,dose1dt)]
summary(full_sample)
full_sample$sex <- as.factor(full_sample$sex)
table(full_sample$sex)
prop.table(table(complete.cases(full_sample[,.(vacc_bus,age,sex,ethn_f,imd_score,chronic,carer,social_care, mdist_min,dose1)])))

#############################  16. Draw figures for accumulated uptake over time ##########################
test0<-data5n[Deceased=="N", .(age=Patient.Age, sex=Sex,lsoa11=LSOA_Code, ethn=EthnicMainGroup,
                              dose1dt,dose1, dose2,dose2dt, ch_flag=NursingCareHomeFlag,ch_match=MatchType,
                              ch_conf=MatchConfidence,COPD,CHD,DIAB,CKD,AST,CAN,OBES,DEP,STRK,
                              distance_min,distance_km, prevposdt, posdt)]
table(as.numeric(test0$distance_km)==0)
table(as.numeric(test0$distance_min)==0)
test0[,log_distkm := ifelse(as.numeric(test0$distance_km)==0,log(as.numeric(distance_km)+0.001),
                           log(as.numeric(distance_km)))]
test0[,log_distmin := ifelse(as.numeric(test0$distance_min)==0,log(as.numeric(distance_min)+0.001),
                            log(as.numeric(distance_min)))]

test0[, posdt:=as.Date(posdt, origin = "1970-01-01")]

test0<-merge(test0,ag_risk_lsoa, by="lsoa11", all.x=T)
test0[, dose1:=as.numeric(dose1)]
test0[is.na(dose1), dose1:=0]

with(test0, table(dose1, useNA = "ifany"))

summary(test0$age)

## Set up age threshold for the sample
test0<-test0[age>=18]

test0[, age_group:=cut(age, breaks=c(18,20,30,40,50,55,60,70,75,80,110), include.lowest=T, right=F)]
table(test0$age_group)
test0[, pc_imd:=ntile(imd_score, n=5)]
test0[, id:=as.numeric(as.factor(lsoa11))]
test0[sex=="U", sex:=NA]
with(test0, table(sex))
test0[, prop_1less:=prop_1less*100]

test0[, ch_flag:=as.numeric(ch_flag=="Y")]
with(test0, prop.table(table(ch_flag)))
with(test0, prop.table(table(ethn)))


test0[ethn=="NULL", ethn:=NA]

# convert to weekly format
test0[, d1_week:=yearweek(dose1dt)]
# first week anyone vaccinated
st_week=min(test0$d1_week, na.rm = T)
# number of weeks in the data se
num_week=max(test0$d1_week, na.rm = T)-min(test0$d1_week, na.rm = T)

#  get age specific population for the whole population 
with(test0, prop.table(table(age_group)))


#  derive age specific uptake for each lsoa
lsoa_age_adj_uptake<-test0[, list(dose1=sum(dose1), dose2=sum(as.numeric(!is.na(dose2))),
                                 uptake=sum(dose1)/.N, uptake2=sum(as.numeric(!is.na(dose2)))/.N, 
                                 pop=.N), 
                          by=.(age_group, lsoa11)]

# get age specific population for the whole population 
lsoa_age_adj_uptake[, tot_pop:=sum(pop), by=.(age_group)]

# apply age specific rates to total population - gives the expected number vaccinated if each lsoa had the pop distribution of all of C&M
lsoa_age_adj_uptake[,c("exp", "exp2") := list(uptake*tot_pop, uptake2*tot_pop)]

lsoa_age_adj_uptake<-lsoa_age_adj_uptake[, list(dose1=sum(dose1), exp=sum(exp), 
                                                dose2=sum(dose2), exp2=sum(exp2), 
                                                tot_pop=sum(tot_pop), pop=sum(pop)), 
                                         by=.(lsoa11)]

lsoa_age_adj_uptake[, c("adj_uptake", "adj_uptake2"):=list(exp*100/tot_pop,exp2*100/tot_pop)]
lsoa_age_adj_uptake[, c("crude_uptake", "crude_uptake2"):=list(dose1*100/pop,dose2*100/pop)]

names(lsoa_age_adj_uptake)

lsoa_age_adj_uptake <- merge(lsoa_age_adj_uptake,ag_risk_lsoa,by="lsoa11", all.x=T)
names(lsoa_age_adj_uptake)

lsoa_age_adj_uptake$adj_nonuptake <- 100- lsoa_age_adj_uptake$adj_uptake
lsoa_age_adj_uptake$nonuptake <- 100- lsoa_age_adj_uptake$crude_uptake


max(test0$dose1dt, na.rm = T)-28
max(test0$posdt, na.rm = T)

table(lsoa_age_adj_uptake$mdist_km==0)
table(lsoa_age_adj_uptake$mdist_min==0)
lsoa_age_adj_uptake$log_distmin<-log(lsoa_age_adj_uptake$mdist_min)
lsoa_age_adj_uptake$log_distkm<-log(lsoa_age_adj_uptake$mdist_km)


summary(lsoa_age_adj_uptake[,c("adj_uptake","adj_nonuptake")])
summary(lsoa_age_adj_uptake[,c("nonuptake","crude_uptake")])

## count the number of LSOAa with sub-group pop <= 5 
lapply(lsoa_age_adj_uptake[,c("pop5n_18_29","pop5n_30_39","pop5n_40_49","pop5n_50_59","pop5n_60_69",
                                        "pop5n_70_79","pop5n_80plus","pop_asian5n","pop_black5n","pop_mixed5n",
                                        "pop_other5n","pop_white5n")], 
                 function(x) data.frame(table(x<=5))) 

table(lsoa_age_adj_uptake$crude_uptake*lsoa_age_adj_uptake$tot_pop<5)
table(lsoa_age_adj_uptake$adj_uptake*lsoa_age_adj_uptake$tot_pop<5)

lsoa_age_adj_uptake[,pop_densx:=scale(pop_dens, center = TRUE, scale = TRUE)]

## Figure 1. The crude uptake rate of each age group by level of deprivation
#Deprivation and age - i.e chart y = crude uptake , x= age and color= deprivation quintile
plot01<-test0[,list(uptake=sum(dose1)*100/.N), by=.(age_group,pc_imd)]

##manually creating the palette
magma_m <- c("#74D055FF", "#FB8761FF", "#B63779FF", "#50127CFF", "#000004FF")



## Figure 3. The age-adjusted accumulated uptake rate for each ethnic group and deprivation group over the weeks

table(test0$age_group)
table(test0$ethn)

### total crude & adj
plot1<-test0[, list(pop=.N), by=.(age_group)]

table(test0$num_week)
num_week
# replicate across week
plot1 <- as.data.table(plot1[rep(seq_len(nrow(plot1)),num_week+1),])
plot1[, week:=st_week-1+ (1:.N), by=.(age_group)]


#  get number vaccinated in each segment
week_test<-test0[is.na(d1_week)==F,list(dose1=sum(dose1)), 
                 by=.(week=d1_week, age_group)]


plot1<-merge(plot1,week_test, by=c("week", "age_group"))

plot1[is.na(dose1)==T, dose1:=0]
plot1[order(week), cum_dose1:=cumsum(dose1), by=.(age_group)]

plot1[, cum_uptake:=cum_dose1/pop]

plot1[, tot_pop:=sum(pop), by=.(age_group, week)]

plot1[, exp:=cum_uptake*tot_pop]



plot1<-plot1[, list(cum_dose1=sum(cum_dose1), exp=sum(exp), 
                    tot_pop=sum(tot_pop), pop=sum(pop)), by=.(week)]


plot1[, adj_uptake:=exp*100/tot_pop]
plot1[, crude_uptake:=cum_dose1*100/pop]

# should be the same 
with(plot1, summary(adj_uptake, crude_uptake))
with(plot1, plot(adj_uptake, crude_uptake))

unique(plot1$week)

# crude uptake
ggplot(data=plot1[week>=yearweek("2020 W50") & week<= yearweek("2021 W29"),],aes(x=week, y=crude_uptake))+
  geom_line(size=1) +
  geom_vline(xintercept=as.numeric(as.Date("2021-02-22")), colour="red", linetype = "dashed") +
  scale_color_manual(values = magma_m[2],name="Crude uptaken (%)") + 
  scale_x_yearweek(date_breaks="1 weeks", date_labels = "%Y %W") +
  guides(x =  guide_axis(angle = 90)) +
  labs(x ="Week number", 
       y="Crude accumulated uptake rate (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face="bold"))

ggsave(normalizePath(file.path(a,"crude.png")), width = 10, height = 8,device="png",
       units = "in", dpi=700)

ggsave(normalizePath(file.path(a,"crude.jpeg")), width = 10, height = 8,device="jpeg",
       units = "in", dpi=400)

ggsave(normalizePath(file.path(a,"crude.svg")), width = 10, height = 8,device="svg",
       units = "in", dpi=700)


### crude uptake with CIs
plot2<-test0[, list(pop=.N), by=.(age_group,lsoa11)]

table(test0$num_week)
num_week
# replicate across week
plot2 <- as.data.table(plot2[rep(seq_len(nrow(plot2)),num_week+1),])
plot2[, week:=st_week-1+ (1:.N), by=.(age_group,lsoa11)]


#  get number vaccinated in each segment
week_test<-test0[is.na(d1_week)==F,list(dose1=sum(dose1)), 
                 by=.(week=d1_week, age_group,lsoa11)]


plot2<-merge(plot2,week_test, by=c("week", "age_group","lsoa11"))

plot2[is.na(dose1)==T, dose1:=0]
plot2[order(week), cum_dose1:=cumsum(dose1), by=.(age_group,lsoa11)]

plot2[, cum_uptake:=cum_dose1/pop, by=.(lsoa11)]

plot2[, tot_pop:=sum(pop), by=.(age_group, week, lsoa11)]

plot2[, exp:=cum_uptake*tot_pop]



plot2<-plot2[, list(cum_dose1=sum(cum_dose1), exp=sum(exp), 
                    tot_pop=sum(tot_pop), pop=sum(pop)), by=.(week,lsoa11)]


plot2[, adj_uptake:=exp*100/tot_pop]
plot2[, crude_uptake:=cum_dose1*100/pop]

# should be the same 
summary(plot2[,.(adj_uptake,crude_uptake)])
with(plot2, plot(adj_uptake, crude_uptake))
plot2[, lapply(.SD, function(x) sum(is.na(x))), .SDcols = 7:8]

unique(plot2$week)
#get the CIs for LSOAs
sd(plot2$crude_uptake)
plot2[,c("mcd_uptake", "sdcd_uptake"):=
        list(mean(crude_uptake),sd(crude_uptake)),by=.(week)]
summary(plot2[,.(adj_uptake,crude_uptake,mcd_uptake,sdcd_uptake)])
plot2[is.na(sdcd_uptake),]
plot2[is.na(sdcd_uptake),sdcd_uptake:=0]
summary(plot2[,.(adj_uptake,crude_uptake,mcd_uptake,sdcd_uptake)])

plot2[,c("mcd_ucl","mcd_lcl"):=list(mcd_uptake + (1.96 * sdcd_uptake),
             mcd_uptake - (1.96 * sdcd_uptake)),by=.(week)]
summary(plot2[,.(adj_uptake,crude_uptake,mcd_uptake,sdcd_uptake,mcd_ucl,mcd_lcl)])
names(plot2)
unique(plot2[mcd_lcl<0,.(week,crude_uptake,mcd_uptake,mcd_lcl,mcd_ucl)])

plot2[,c("mcd_max","mcd_min"):=list(max(crude_uptake),min(crude_uptake)),by=.(week)]
summary(plot2[,.(adj_uptake,crude_uptake,mcd_uptake,sdcd_uptake,mcd_ucl,mcd_lcl,mcd_max,mcd_min)])
unique(plot2[,.(week,mcd_uptake,sdcd_uptake,mcd_ucl,mcd_lcl,mcd_max,mcd_min)])

plot2<-merge(unique(plot2[,.(week,mcd_uptake,sdcd_uptake,mcd_ucl,mcd_lcl,mcd_max,mcd_min)]),
             plot1,by="week")
class(plot2$week)

ggplot(data=plot2[week>=yearweek("2020 W50") & week<= yearweek("2021 W29"),],
       aes(x=week, y=crude_uptake))+
  geom_line(size=1, color="#440154FF") +
  geom_vline(xintercept=as.numeric(as.Date("2021-02-22")), colour="red", linetype = "dashed") +
  geom_ribbon(aes(ymin=mcd_min, ymax=mcd_max, fill="#440154FF"),size=0, alpha=0.3)+
  scale_colour_manual(values = c("#440154FF"), 
                      labels=c("Crude uptake")) +
  scale_fill_manual(values = c("#440154FF"), 
                    labels=c("Crude uptake")) +
  xlab("Week number") + ylab("Crude accumulated vaccination rate %")+
  theme_classic()+
  theme(text = element_text(size=16),legend.position="bottom", 
        legend.box = "horizontal", legend.text = element_text(size=24),
        legend.title = element_blank(),legend.key.size = unit(1, "cm"), 
        axis.text.x = element_text(angle = 45,vjust=0.5))+ 
  theme(strip.background = element_blank())

### try 95% CIs again 
plot2[mcd_lcl<0,mcd_lcl:=0]

ggplot(data=plot2[week>=yearweek("2020 W50") & week<= yearweek("2021 W29"),],
       aes(x=week, y=mcd_uptake,
           fill=NA,
           color="#440154FF"))+
  geom_line(size=1) +
  geom_vline(xintercept=as.numeric(as.Date("2021-02-22")), colour="red", linetype = "dashed") +
  geom_ribbon(aes(ymin=mcd_lcl, ymax=mcd_ucl,
                  fill="bands"),size=0, alpha=0.3)+
  scale_fill_manual(name='', values=c("#440154FF"),
                    labels=c("Crude uptake of LSOAs (95% confidence intervals)")) +
  scale_colour_manual(name='', values=c("#440154FF"),
                      labels=c("Average crude uptake across all LSOAs")) +
  scale_x_yearweek(date_breaks="1 weeks", date_labels = "%Y %W") +
  guides(x =  guide_axis(angle = 90)) +
  xlab("Week number") + ylab("Crude accumulated uptake rate (%)")+
  theme_classic()+
  theme(text = element_text(size=16),
        legend.position="bottom", 
        legend.box = "horizontal", legend.text = element_text(size=18),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"))+ 
  theme(strip.background = element_blank())



ggsave(normalizePath(file.path(a,"crude_uptake_lsoa_ci.png")), width = 16, height = 12,device="png",
       units = "in", dpi=700)

ggsave(normalizePath(file.path(a,"crude_uptake_lsoa_ci.jpeg")), width = 16, height = 12,device="jpeg",
       units = "in", dpi=400)

ggsave(normalizePath(file.path(a,"crude_uptake_lsoa_ci.svg")), width = 16, height = 12,device="svg",
       units = "in", dpi=700)




### crude & adj by ethnicity
plot3<-test0[is.na(ethn)==F, list(pop=.N), by=.(age_group,ethn)]
plot3<- plot3[is.na(ethn)==F]

table(test0$num_week)
num_week
# replicate across week
plot3 <- as.data.table(plot3[rep(seq_len(nrow(plot3)),num_week+1),])
plot3[is.na(ethn)==F , week:=st_week-1+ (1:.N), by=.(age_group,ethn)]


#  get number vaccinated in each segment
week_test<-test0[is.na(d1_week)==F & is.na(ethn)==F,list(dose1=sum(dose1)), 
                by=.(week=d1_week, age_group,ethn)]


plot3<-merge(plot3,week_test, by=c("week", "age_group","ethn"))

plot3[is.na(dose1)==T, dose1:=0]
plot3[order(week), cum_dose1:=cumsum(dose1), by=.(ethn, age_group)]

plot3[, cum_uptake:=cum_dose1/pop]

plot3[, tot_pop:=sum(pop), by=.(age_group, week)]

plot3[, exp:=cum_uptake*tot_pop]


plot3<-plot3[, list(cum_dose1=sum(cum_dose1), exp=sum(exp), 
                    tot_pop=sum(tot_pop), pop=sum(pop)), by=.(ethn,week)]


plot3[, adj_uptake:=exp*100/tot_pop]
plot3[, crude_uptake:=cum_dose1*100/pop]

# should be similar 
with(plot3, summary(adj_uptake, crude_uptake))
with(plot3, plot(adj_uptake, crude_uptake))

unique(plot3$week)

# crude uptake
ggplot(data=plot3[week>=yearweek("2020 W50") & week<= yearweek("2021 W29"),],aes(x=week, y=crude_uptake,
                                                   group=as.factor(ethn),
                                                   color=as.factor(ethn)))+
  geom_line(size=1) +
  geom_vline(xintercept=as.numeric(as.Date("2021-02-22")), colour="red", linetype = "dashed") +
  scale_color_manual(values = magma_m,name="Ethnicity") + 
  scale_x_yearweek(date_breaks="1 weeks", date_labels = "%Y %W") +
  guides(x =  guide_axis(angle = 90)) +
  labs(x ="Week number", 
       y="Crude accumulated uptake rate (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face="bold"))

ggsave(normalizePath(file.path(a,"crude_eth.png")), width = 10, height = 8,device="png",
       units = "in", dpi=700)

ggsave(normalizePath(file.path(a,"crude_eth.jpeg")), width = 10, height = 8,device="jpeg",
       units = "in", dpi=400)

ggsave(normalizePath(file.path(a,"crude_eth.svg")), width = 10, height = 8,device="svg",
       units = "in", dpi=700)


#age-adjusted uptake
ggplot(data=plot3[week>=yearweek("2020 W50") & week<= yearweek("2021 W29"),],aes(x=week, y=adj_uptake,
                                                                                 group=as.factor(ethn),
                                                                                 color=as.factor(ethn)))+
  geom_line(size=1) +
  geom_vline(xintercept=as.numeric(as.Date("2021-02-22")), colour="red", linetype = "dashed") +
  scale_color_manual(values = magma_m,name="Ethnicity") + 
  scale_x_yearweek(date_breaks="1 weeks", date_labels = "%Y %W") +
  guides(x =  guide_axis(angle = 90)) +
  labs(x ="Week number", 
       y="Age-adjusted accumulated uptake rate (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face="bold"))
ggsave(normalizePath(file.path(a,"adj_eth.png")), width = 10, height = 8,device="png",
       units = "in", dpi=700)

ggsave(normalizePath(file.path(a,"adj_eth.jpeg")), width = 10, height = 8,device="jpeg",
       units = "in", dpi=400)

ggsave(normalizePath(file.path(a,"adj_eth.svg")), width = 10, height = 8,device="svg",
       units = "in", dpi=700)



#Age adjusted uptake over time (weeks) by IMD  - y= age adjusted uptake , x=week , color = ethnicity
# derive age specific uptake
# get age specific population for the whole population 
table(test0$pc_imd)
test0[, pc_imd:=ntile(imd_score, n=3)]
table(test0$pc_imd)
plot4<-test0[is.na(pc_imd)==F, list(pop=.N), by=.(age_group,pc_imd)]
plot4<- plot4[is.na(pc_imd)==F]


# replicate across week
plot4 <- as.data.table(plot4[rep(seq_len(nrow(plot4)),num_week+1),])
plot4[is.na(pc_imd)==F , week:=st_week-1+ (1:.N), by=.(age_group,pc_imd)]


#  get number vaccinated in each segment
week_test<-test0[is.na(d1_week)==F & is.na(pc_imd)==F,list(dose1=sum(dose1)), 
                by=.(week=d1_week, age_group,pc_imd)]


plot4<-merge(plot4,week_test, by=c("week", "age_group","pc_imd"))

plot4[is.na(dose1)==T, dose1:=0]
plot4[order(week), cum_dose1:=cumsum(dose1), by=.(pc_imd, age_group)]

plot4[, cum_uptake:=cum_dose1/pop]

plot4[, tot_pop:=sum(pop), by=.(age_group, week)]

plot4[, exp:=cum_uptake*tot_pop]


plot4<-plot4[, list(cum_dose1=sum(cum_dose1), exp=sum(exp), 
                    tot_pop=sum(tot_pop), pop=sum(pop)), by=.(pc_imd,week)]


plot4[, adj_uptake:=exp*100/tot_pop]
plot4[, crude_uptake:=cum_dose1*100/pop]

# should be similar 
with(plot4, summary(adj_uptake, crude_uptake))
with(plot4, plot(adj_uptake, crude_uptake))


unique(plot4$week)

ggplot(data=plot4[week>=yearweek("2020 W50"),], aes(x=week, y=crude_uptake,group=as.factor(pc_imd),color=as.factor(pc_imd)))+
  geom_line(size=1) +
  geom_vline(xintercept=as.numeric(as.Date("2021-02-22")), colour="red", linetype = "dashed") +
  scale_color_manual(values = magma_m[1:3],name="Deprivation",
                     labels = c("Least deprived", "Intermediate deprivation", "Most deprived")) + 
  scale_x_yearweek(date_breaks="1 weeks", date_labels = "%Y %W") +
  guides(x =  guide_axis(angle = 90)) +
  labs(x ="Week number", 
       y="Crude accumulated uptake rate (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face = "bold"))
ggsave(normalizePath(file.path(a,"crude_depr.png")), width = 10, height = 8,device="png",
       units = "in", dpi=700)

ggsave(normalizePath(file.path(a,"crude_depr.jpeg")), width = 10, height = 8,device="jpeg",
       units = "in", dpi=400)

ggsave(normalizePath(file.path(a,"crude_depr.svg")), width = 10, height = 8,device="svg",
       units = "in", dpi=700)


ggplot(data=plot4[week>=yearweek("2020 W50"),], aes(x=week, y=adj_uptake,group=as.factor(pc_imd),color=as.factor(pc_imd)))+
  geom_line(size=1) +
  geom_vline(xintercept=as.numeric(as.Date("2021-02-22")), colour="red", linetype = "dashed") +
  scale_color_manual(values = magma_m[1:3],name="Deprivation",
                     labels = c("Least deprived", "Intermediate deprivation", "Most deprived")) + 
  scale_x_yearweek(date_breaks="1 weeks", date_labels = "%Y %W") +
  guides(x =  guide_axis(angle = 90)) +
  labs(x ="Week number", 
       y="Age-adjusted accumulated uptake rate (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face = "bold"))
ggsave(normalizePath(file.path(a,"adj_depr.png")), width = 10, height = 8,device="png",
       units = "in", dpi=700)

ggsave(normalizePath(file.path(a,"adj_depr.jpeg")), width = 10, height = 8,device="jpeg",
       units = "in", dpi=400)

ggsave(normalizePath(file.path(a,"adj_depr.svg")), width = 10, height = 8,device="svg",
       units = "in", dpi=700)

### total crude uptake by LAs
plot5<-test0[, list(pop=.N), by=.(age_group,ladnm)]

table(test0$num_week)
num_week
# replicate across week
plot5 <- as.data.table(plot5[rep(seq_len(nrow(plot5)),num_week+1),])
plot5[, week:=st_week-1+ (1:.N), by=.(age_group,ladnm)]


#  get number vaccinated in each segment
week_test<-test0[is.na(d1_week)==F,list(dose1=sum(dose1)), 
                 by=.(week=d1_week, age_group,ladnm)]


plot5<-merge(plot5,week_test, by=c("week", "age_group", "ladnm"))

plot5[is.na(dose1)==T, dose1:=0]
plot5[order(week), cum_dose1:=cumsum(dose1), by=.(age_group,ladnm)]

plot5[, cum_uptake:=cum_dose1/pop]

plot5[, tot_pop:=sum(pop), by=.(age_group, week,ladnm)]

plot5[, exp:=cum_uptake*tot_pop]



plot5<-plot5[, list(cum_dose1=sum(cum_dose1), exp=sum(exp), 
                    tot_pop=sum(tot_pop), pop=sum(pop)), by=.(week,ladnm)]


plot5[, adj_uptake:=exp*100/tot_pop]
plot5[, crude_uptake:=cum_dose1*100/pop]

# should be the same 
with(plot5, summary(adj_uptake, crude_uptake))
with(plot5, plot(adj_uptake, crude_uptake))

unique(plot5$week)
table(plot5$ladnm)

# crude uptake
magma9 <- c(magma_m, "#fde725","#7e03a8","#0d0887","#21918c")
ggplot(data=plot5[week>=yearweek("2020 W50") & week<= yearweek("2021 W29"),],
       aes(x=week, y=crude_uptake,group=as.factor(ladnm),color=as.factor(ladnm)))+
  geom_line(size=1) +
  geom_vline(xintercept=as.numeric(as.Date("2021-02-22")), colour="red", linetype = "dashed") +
  scale_color_manual(values = magma9,name="Crude uptake (%)") + 
  scale_x_yearweek(date_breaks="1 weeks", date_labels = "%Y %W") +
  guides(x =  guide_axis(angle = 90)) +
  labs(x ="Week number", 
       y="Crude accumulated uptake rate (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, face="bold"))

ggsave(normalizePath(file.path(a,"crude_lad.png")), width = 10, height = 8,device="png",
       units = "in", dpi=700)

ggsave(normalizePath(file.path(a,"crude_lad.jpeg")), width = 10, height = 8,device="jpeg",
       units = "in", dpi=400)

ggsave(normalizePath(file.path(a,"crude_lad.svg")), width = 10, height = 8,device="svg",
       units = "in", dpi=700)

##################################  17. Summary statistics for the full non-intervention population ##########################

### generate summary statistics of the full non-intervention population

summary(test0$cipha_pop)
summary(test0$mean_age)
summary(lsoa_age_adj_uptake)
NROW(unique(lsoa_age_adj_uptake$lsoa11))

NROW(unique(test2$lsoa11))
NROW(unique(test2[vacc_bus==1,.(lsoa11)]))
NROW(unique(test0[lsoa11%in%unique(test2[vacc_bus==1,.(lsoa11)]),.(lsoa11)]))
NROW(unique(test0[lsoa11%notin%unique(test2[vacc_bus==1,.(lsoa11)]),.(lsoa11)]))
NROW(unique(test0[lsoa11%notin%unique(test2[vacc_bus==1]$lsoa11),.(lsoa11)]))

summary(test0$prop_1less)
summary(test0$prop_female)
nonintv <- test0[lsoa11%notin%unique(test2[vacc_bus==1]$lsoa11),
                                          list(cipha_pop=.N,
                                          prop_female=weighted.mean(prop_female,w=cipha_pop)*100,
                                          pop_dens=weighted.mean(pop_dens,w=cipha_pop),
                                          mean_age=weighted.mean(as.numeric(mean_age),w=cipha_pop),
                                          prop_asian5n=weighted.mean(prop_asian5n,w=cipha_pop)*100,
                                          prop_black5n=weighted.mean(prop_black5n,w=cipha_pop)*100,
                                          prop_mixed5n=weighted.mean(prop_mixed5n,w=cipha_pop)*100,
                                          imd_score=weighted.mean(imd_score,w=cipha_pop),
                                          prop_1less=weighted.mean(prop_1less,w=cipha_pop),
                                          mdist_min=weighted.mean(mdist_min),
                                          x_numlsoa=NROW(unique(lsoa11)))]

nonintv<-as.data.table(t(nonintv), keep.rownames=T)


nonintv[,V1 := round(V1,2)]
## this is to report more decimal places as Dr Philip Britteon requested


write.csv(nonintv,
          file= normalizePath(file.path(a,"nonintervention_var_summary.csv")))


NROW(unique(lsoa_panel$lsoa11))
lsoa_panel[, time:=as.numeric(as.factor(week_from_vaccbus))]

lsoa_panel[,list(vacc_rate=weighted.mean(vacc_rate,w=cipha_pop),
                       cum_vacc_rate=weighted.mean(cum_vacc_rate[week_from_vaccbus==-1], w=cipha_pop[week_from_vaccbus==-1]),
                       x_numlsoa=NROW(unique(lsoa11))),by=.(vacc_bus)]

class(lsoa_panel$week)
summary(lsoa_panel[lsoa11%notin%unique(test2[vacc_bus==1]$lsoa11) & week>=yearweek("2021 W08") & week<=yearweek("2021 W29"),
                   .(cum_vacc_rate,vacc_rate,cum_dose1,cipha_pop)])


lsoa_panel[lsoa11%notin%unique(test2[vacc_bus==1]$lsoa11) & week>=yearweek("2021 W08") & week<=yearweek("2021 W29") ,
            list(vacc_rate=weighted.mean(vacc_rate,w=cipha_pop),
                 cum_vacc_rate=weighted.mean(cum_vacc_rate, w=cipha_pop),
                 x_numlsoa=NROW(unique(lsoa11)))]

############################# 18. Check the predictive power on preintervention outcome ##########################

a1_appendix<- as.data.table(test2[time<8,
                               list(tot_dose1=sum(dose1)),by=.(vacc_bus,time)])

a1_appendix <- dcast(a1_appendix, time ~ vacc_bus, value.var = "tot_dose1")

write.csv(a1_appendix,
          file= normalizePath(file.path(a,"a1_appendix.csv")))

####################### 19. Sensitivity - testing distance thresholds and spatial spill over effects #########################
#Simply switch on relevant lines of codes (as commented) to generate the results; 
#please refer to the separate script "mobile_vax_lookup_CDCG.R" for more details on how the lsoa lookup table is generated.


####################### 20. Sensitivity - excluding pop-in sites#########################
#Simply switch on relevant lines of codes (as commented) to generate the results; 
#please refer to the separate script "mobile_vax_lookup_CDCG.R" for more details on how the specific lsoa lookup table is generated.


####################### 21. sensitivity - individual-level survival analysis  #########################
# construct survey data
# people unvaccinated a week before the intervention date, including:
# st_week>=-1, people hadn't vaccinated before the intervention date but received the vaccine by the end of the study period;
# is.na(st_week), people hadn't vaccinated either before the intervention date or by the end of the study period;

table(surv_data$week_from_vaccbus)


fit <- survfit(Surv(week_from_vaccbus, outcome2) ~ vacc_bus, data =surv_data,weights = main)
summary(fit)
class(fit)
str(fit)

# compare with the synthetic control results
summary(sea2)
class(sea2)
str(sea2)

table(surv_data$vacc_bus)
a3_appendix<- as.data.table(test2[time>=8,
                                  list(tot_dose1=sum(dose1*Main)),by=.(vacc_bus,time)])
a3_appendix <- dcast(a3_appendix, time ~ vacc_bus, value.var = "tot_dose1")
table(test2[time>=8,list(tot_dose1=sum(dose1*Main)),by=.(vacc_bus,time)])
names(a3_appendix)[2:3] <- c("no","yes")
a3_appendix[,unvacc_no:=NROW(surv_data[vacc_bus==0])-cumsum(no)]
a3_appendix[,unvacc_yes:=NROW(surv_data[vacc_bus==1])-cumsum(yes)]

write.csv(a3_appendix,
          file= normalizePath(file.path(a,"papers/vaccine_bus/BMJ_open/a3_appendix.csv")))



test2[time>8 & vacc_bus==1]$dose1*test2[time>8 & vacc_bus==1]$Main

p<- ggsurvplot(fit,
           xlim = c(0, 4),
           break.x.by=1,
           xlab = "Weeks since vaccine bus",
           ylim=c(0.5,1),
           surv.scale="percent",
           ylab = "Survival probability\n(the probability of adults to stay unvaccinated)",
           pval = TRUE, 
           pval.method = TRUE,
           conf.int = 0.95,
 #          font.main = 10,
           font.x = 8,
           font.y = 8,
           font.tickslab = 7,
           risk.table = TRUE, 
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#440154FF","#FDE725FF"), 
           legend.title = "",
           legend.labs = c("Synthetic Control", "Intervention"),
           risk.table.title = "Number at risk (the number of unvaccinated adults)",
#           risk.table = TRUE,
           fontsize = 7,
           tables.height = 0.2)

p1 = p$plot
p2 = p$table
class(p$table)
str(p$table)
plotp = cowplot::plot_grid(p1,p2,align = "v",ncol =1,rel_heights = c(4,1))


ggsave(normalizePath(file.path(a,"survival.png")), 
       survminer:::.buildggsurvplot(p), width = 16, height = 14,device="png",
       units = "in", dpi=700)


ggsave(normalizePath(file.path(a,"survival.png")), plot = plotp,
       width = 16, height = 14,device="png",
       units = "in", dpi=700)

ggsave(normalizePath(file.path(a,"survival.jpeg")),  plot = plotp,
       width = 16, height = 14,device="jpeg",
       units = "in", dpi=400)

ggsave(normalizePath(file.path(a,"survival.svg")),  plot = plotp,
       width = 16, height = 14,device="svg",
       units = "in", dpi=700)



cox <- coxph(Surv(week_from_vaccbus, outcome2) ~ pc_imd*vacc_bus+ethn_f*vacc_bus+agegrp*vacc_bus+
               sex+chronic+carer+social_care+mdist_min, data =surv_data,weights = main)
summary(cox)
str(cox)
str(summary(cox))

pvc <- round(coef(summary(cox))[,6],3)
hr_ci <- round(summary(cox)$conf.int,2)

class(hr_ci)
str(hr_ci)

res_cox <- as.data.table(hr_ci)
res_cox[,var:=rownames(hr_ci)]
res_cox[,`exp(-coef)`:=NULL]
setcolorder(res_cox, c(4,1,2,3))
res_cox[,p:=as.data.table(pvc)]

write.csv(res_cox, normalizePath(file.path(a,"appendix_cox.csv")))


####################### 22. calculate sample sizes for the subgroup analysis  #########################

prop.table(table(surv_data$ethn_f, surv_data$vacc_bus), 2)


summary(surv_data[,list(pc_imd, ethn_f, agegrp,vacc_bus)])
table(surv_data[,list(pc_imd, ethn_f, agegrp),by=.(vacc_bus)])

subgrp_sz <- as.data.table(surv_data[vacc_bus==1 & !is.na(agegrp) & !is.na(ethn_f),
                                     list(subgrp_sz=NROW(vacc_bus),vacc_bus), 
                                     by=.(pc_imd, ethn_f, agegrp)])


ggplot(data = subgrp_sz, mapping = aes(x = pc_imd,
                                      y = agegrp,
                                      fill = subgrp_sz)) +
  scale_fill_viridis(name="Sample size",
                     begin = 0.1, end = 1,
                     na.value="white") +
  geom_tile(aes(fill = subgrp_sz)) +
  geom_text(aes(label = subgrp_sz), fontface = "bold", size = 12,hjust="left") +
  facet_grid(ethn_f~vacc_bus)+
  theme_classic() +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"))+
  xlab(label = "Tercile of deprivation (from least to most deprived)") +
  ylab(label = "Age group")
#
#
ggsave(normalizePath(file.path(a,"subgrp_sz_heatmap.png")), width = 16, height = 14,device="png",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"subgrp_sz_heatmap.svg")), width = 16, height = 14,device="svg",
       units = "in", dpi=700)
ggsave(normalizePath(file.path(a,"subgrp_sz_heatmap.jpeg")), width = 16, height = 14,device="jpeg",
       units = "in", dpi=400)
