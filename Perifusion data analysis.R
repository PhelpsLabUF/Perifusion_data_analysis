############# Perifusion data for insulin secretion - nPOD ###############

# Written by: Alexandra Cuaycal
## This R script was adapted from live Ca imaging analysis by AC

## set working directory to folder containing perifusion data
#setwd()

## load packages
library(devtools)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggbreak) 
library(patchwork)
library(viridis)
library(DescTools)
library(tidyr)
library(agricolae)
library(ggpubr)
library(multcompView)


## load table 
data <- read.csv("Perifusion_data_FINAL.csv")

### loading donor metdata 
donor_info <- read.csv("Donor_info.csv")

# Cases to exclude:
## 6597, 6558 because there was issues in the perifusion experiment 
## 6535, 6539 and 6562 because values out of ranges not corrected

# Cases with <3 values out of range (values excluded, set to NA)
#  6552, 6582, 6512, 6505, 6518

# case 6516 had two values with low perfusate (excluded, set to NA) and one out of range (excluded,  set to NA)
# case 6530, one value with low perfusate (excluded,  set to NA)

# negative values (under the detection limit) are set to 0


excluded_cases <- c("Case_6597", "Case_6558", "Case_6535", "Case_6539", "Case_6562")

# Removing excluded cases
data <- data[, !( colnames(data) %in% excluded_cases) ]

# For Donor metdata table
excluded_cases_donor_info <- as.numeric(mapply(excluded_cases, FUN=function(x){
  a<-strsplit(x,split = "_")
  a<-unlist(a)
  a<- a[2]
  
} ))

#excluding cases in donor metdata
donor_info<- donor_info[!(donor_info$Case_ID %in% excluded_cases_donor_info ),]

# Export final metdata with donors analyzed
write.csv(donor_info, file = "donor_info_insulin_data.csv")

## Setting negative values to 0
data <-as.data.frame(apply(data, 2, function(x) if_else(x < 0, 0, x)) )

## setting up time in stimulations (total time = 80 min)

## LG is 10 min
## HG is 20 min
## LG is 30 min
## KCl is 5 min
## LG is 15 min

stimulations_times=c(10,20,30,5, 15)
# name of stimulations
stimulations<- c("LG_1", "HG", "LG_2", "KCl", "LG_3")


## The following function gets the index for each stimulation

## Time correction= 3 min delay (determined by nPOD)
# framerate = 1 because it is measuring one value per min
# stimulations_time_min = stimulation_times
# name_stimulations = stimulations

time_stimulations <- function(data, framerate=1, stimulations_time_min=stimulations_times, name_stimulations = stimulations , time_correction=3){
  
  # naming vector with stimulations times
  names(stimulations_time_min) <- name_stimulations
 
  
  # time vector with time per frame number in min
  time <- seq(from = 1, to = nrow(data), by=1 )
  
  # total time 
  time_min <- framerate*nrow(data)
  
  ## getting time points in sequence
  ## i.e. additive stimulation times
  
  if (length(stimulations_time_min) == 5){
    time_in_seq <- c(stimulations_time_min[1], stimulations_time_min[1]+stimulations_time_min[2], 
                     stimulations_time_min[1]+stimulations_time_min[2]+stimulations_time_min[3], 
                     stimulations_time_min[1]+stimulations_time_min[2]+stimulations_time_min[3]+ stimulations_time_min[4], last(time))
    
    names(time_in_seq)<- names(stimulations_time_min)
    
  } else{
    # if n stimulations = 4 i.e. LG, HG, LG, KCl
    time_in_seq <- c(stimulations_time_min[1], stimulations_time_min[1]+stimulations_time_min[2], 
                     stimulations_time_min[1]+stimulations_time_min[2]+stimulations_time_min[3], last(time))
    
    names(time_in_seq)<- names(stimulations_time_min)
    
  }
  
  
  ## correcting time in sequence with time correction (delay)
  
  if (is.null (time_correction) == F){
    time_in_seq <- time_in_seq + time_correction
    
    # last stimulation duration is subtracted the time correction
    time_in_seq [length(time_in_seq)] <- last(time_in_seq)-time_correction
  } 
  
  ## getting the frames/indexes for each stimulation
  stimulation_frames <- list()
  
  if (length(stimulations_time_min) == 5){
    
    for (i in 1:length(time_in_seq)){
      
      if (i==1){
        stimulation_frames[[i]] = time[time <= time_in_seq[i]]
      }
      else if (i==length(time_in_seq)){
        stimulation_frames[[i]] = time[time > time_in_seq[i-1]]
      }
      else{
        stimulation_frames[[i]] = time[time <= time_in_seq[i] & time > time_in_seq[i-1]]
      }
    }
  }else{
    # if n stimulations = 4
    for (i in 1:length(time_in_seq)){
      
      if (i==1){
        stimulation_frames[[i]] = time[time >= time_in_seq[i] & time < time_in_seq[i+1]]
      }
      else if (i==length(time_in_seq)){
        stimulation_frames[[i]] = time[time > time_in_seq[i]]
      }
      else stimulation_frames[[i]] = time[time >= time_in_seq[i] & time < time_in_seq[i+1]]
      
    }
    
  }
  
  names(stimulation_frames) <- names(stimulations_time_min)
  
  ## returns the time vector, the additive stimulation times, and stimulation time indexes/frames
  results<- list(time_axis = time, time_points = time_in_seq, stimulation_frames = stimulation_frames)
  return(results)
  
}


## using function to get the index for each stimulation
time_stimulations_data <- time_stimulations(data =data, stimulations_time_min = stimulations_times, name_stimulations = stimulations)

#stimulation time indexes/frames
stimulation_frames<- time_stimulations_data$stimulation_frames

#time vector
time<- time_stimulations_data$time_axis

#additive times
time_points <- time_stimulations_data$time_points



################# normalizing data to first LG stimulation  ###############################

# data to normalize to 3G (i.e. stimulation index (SI))
data_for_si_3G <- data[stimulation_frames$LG_1,]

## mean of selected stimulation/baseline for normalization
# it needs donor_info in case a case is excluded because baseline = 0

data_to_norm <- function(data_for_norm=data_for_si_3G, data=data,donor_info=donor_info){
  
  ## average 3G per donor
  for_norm_mean <- apply(data_for_norm, 2, function(x) mean(x, na.rm = TRUE))

  ## if we have a mean of zero, we have to exclude that case, otherwise it will
  # create infinite values when normalizing data
  
  if (any(for_norm_mean== 0) == T){
    case_excluded <- unlist(strsplit(names(for_norm_mean[ for_norm_mean==0 ]), split="_"))[2]
    
    ##excluding that value from the 3G averages
    for_norm_mean <- for_norm_mean[! for_norm_mean==0 ]
    
    # excluding case from the data and metdata
    data <- data[, names(for_norm_mean)]
    donor_info <- subset(donor_info, Case_ID != case_excluded )
    
    # print excluded donors
    print(case_excluded)
  }
  
  ## normalizing data to 3G stimulation
  data_norm <-  mapply('/', data,for_norm_mean)
  
  # returns norm data and updated metdata
  results<-list(data_norm=data_norm, donor_info=donor_info)
  return(results)
  
}

# using function to normalize data
data_norm_3G <- data_to_norm(data_for_norm = data_for_si_3G, data = data, donor_info = donor_info)$data_norm
donor_info_norm_3G <- data_to_norm(data_for_norm = data_for_si_3G, data = data, donor_info = donor_info)$donor_info


################## Function to calculate AUC given time intervals ###################

## function takes norm_data and the index for the time intervals from which AUC is calculated

## we need the index in a given time interval/stimulation to extract the data from norm_data:

# we are going to calculate AUC for LG_1, HG, and KCl stimulations
frames_for_auc <- stimulation_frames[ which( !(names(stimulation_frames) %in%  c("LG_2","LG_3"))) ]

# Getting the indexes for each time interval 
index_per_stimulation <- list()

for (y in 1:length(frames_for_auc)){
  
  index_per_stimulation[[y]] <- which(time %in% frames_for_auc[[y]])
  
}
names(index_per_stimulation)<- names(frames_for_auc)

## defining function
auc_roi <- function(data = data, index_interval = index_per_min){
  
  auc_per_interval <-NULL #will create a matrix where result will be saved (one element per ROI/Donor)
  
  for (r in 1:ncol(data)){
    auc_interval<-NULL
    for (k in 1:length(index_interval)){
      
      # calculating AUC for each time interval per ROI/donor
      auc_individual <- AUC(x= time[index_interval[[k]]],y= data[index_interval[[k]], r],
                            na.rm = TRUE, method = "trapezoid")
      # appending data in a vector )
      auc_interval <- c(auc_interval, auc_individual)
      
      
    }
    #assigning names to elements in auc_interval so we know to which stimulation it corresponds
    names(auc_interval) <- names(index_interval)
    auc_per_interval<- cbind(auc_per_interval, auc_interval)
    
  }
  
  ## naming elements in list auc_per_min based on ROI/Donor
  auc_per_interval <- as.data.frame(auc_per_interval)
  colnames(auc_per_interval) <- colnames(data)
  rownames(auc_per_interval) <- names(index_interval)
  
  ## calculating average AUC for all ROIs (only used for Ca imaging)
  
  average_auc <- apply(auc_per_interval, 1, function(x) mean(x, na.rm=TRUE))
  
  # returns auc per interval/stimulation and average auc (for ROIs in Ca imaging)
  results <- list(auc_per_interval = auc_per_interval, average_auc = average_auc)
  
  return(results)
  
  
}



############### Calculating AUC for each case for LG_1, HG, KCl stimulations with response correction ###

# using function
auc_data_norm_to_3G <- auc_roi(data = data_norm_3G, index_interval = index_per_stimulation)

# getting the AUC per interval table
auc_per_case_corr_3G <- auc_data_norm_to_3G$auc_per_interval


######## calculating AUC per min ##################

# first we determined the real time that each stimulation lasted (based on correction)

time_each_stimulation <- unlist(lapply(index_per_stimulation, length))


## we divide AUC by that time
auc_per_case_corr_3G <- auc_per_case_corr_3G / time_each_stimulation

############### converting to long data frame for plots #########################333

auc_per_case_f_3G <- data.frame(Stimulation = names(index_per_stimulation), auc_per_case_corr_3G)

### making long tables ##
auc_per_case_f_3G <- auc_per_case_f_3G %>% pivot_longer(cols = ! Stimulation, names_to ="Case_ID", values_to = "AUC")

## getting only case numbers for Case_ID column
case_number_3G <- NULL
for (y in 1:length(auc_per_case_f_3G$Case_ID)){
  
  case <- unlist(strsplit(auc_per_case_f_3G$Case_ID[y], split = "_"))
  case_number_3G <- c(case_number_3G, case[2])
} 

## replacing in Case_ID column

auc_per_case_f_3G$Case_ID <- factor(case_number_3G, levels = unique(case_number_3G)) 
auc_per_case_f_3G$Stimulation <- factor(auc_per_case_f_3G$Stimulation, levels = c("LG_1", "HG", "KCl"))


### adding clinical_phenotype column to AUC data
# we subset donor info for Case_Id and Clinical Phenotype only
donor_info_subset  <- donor_info_norm_3G[,1:2]

# merging to AUC/min table
auc_per_case_f_3G <- merge(auc_per_case_f_3G, donor_info_subset)

auc_per_case_f_3G$Clinical_phenotype <- factor(auc_per_case_f_3G$Clinical_phenotype, levels = c("ND", "sAAb+", "mAAb+", "T1D"))

auc_per_case_f_3G <- auc_per_case_f_3G %>% arrange(Clinical_phenotype, Stimulation)

## exporting AUC/min tables for graphpad
write.csv(auc_per_case_f_3G, file = "AUC_per_min_normalized_to_3G.csv")

