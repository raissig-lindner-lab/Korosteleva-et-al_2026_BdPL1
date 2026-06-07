#POME script for intensity plot
install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("paletteer")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(paletteer)

rm(list = ls()) #remove everything from the global environment

print(getwd()) #check that your working directory corresponds to the location of the source script
#install.packages("rstudioapi") # for RStudio users
#library(rstudioapi) # for RStudio users
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to the file location


Reconstruction <- function(x) # Move all NA measurements of each angel to the end of the row
{num.na <- sum(is.na(x))
x <- x[!is.na(x)]
x <- c(x, rep(NA, num.na))
return(x)}


intensity <- function(file_name) {
  Cell <-  readr::read_csv(file_name, na = c("NaN"))# Convert NaN from FIJI output to NA
  Cell$Angle <- c(1:63)*360/63
  Cell <- data.frame(t(apply(Cell, 1, Reconstruction))) # Move NA to the end of the row
  Cell <- dplyr::select(Cell,c(2,5,6,7,8)) 
  colnames(Cell) <- c("Angle",9:12)# Select only the first 4 pixels to represent the membrane region
  Cell$mean_val <- rowMeans(Cell[ , c(2:5)], na.rm=TRUE)
  
  cellmin <- min(Cell$mean_val)
  cellmax <- max(Cell$mean_val)
  percentage <- ((Cell$mean_val - cellmin) / (cellmax - cellmin) * 100)
  percentage <- as.data.frame(percentage)
  Cell <- cbind(Cell, percentage)
  
  Cell
}

files <- list.files(path = "./data/01_PL1-PL1", pattern = ".csv", full.names = TRUE) #set path to a folder
l <- length(files)

for(i in 1:l){
  x <- intensity(files[i])
  assign(paste0("Cell",i),x,.GlobalEnv)
}

# Dynamically create a list of all data frames whose names start with "Cell"
dfs <- mget(ls(pattern = "^(Cell)\\d*"))

# Extract the percentage column from each data frame in the list
columns <- lapply(dfs, function(df) df[, 7])

# Combine the extracted columns 
Cell_all <- cbind(Cell1[, c(1, 7)], do.call(cbind, columns))
Cell_all$percentage <- NULL

col_num <- ncol(Cell_all)
Cell_all$mean_all <- rowMeans(Cell_all[, c(2:col_num)], na.rm = TRUE) #for percentage, not for the raw data
Cell_all$median_all <- apply(Cell_all[, 2:col_num], 1, median, na.rm = TRUE)

get_mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

Cell_all$mode_all <- apply(Cell_all[, 2:col_num], 1, get_mode)

# Extract the raw mean
columns_r <- lapply(dfs, function(df) df[, 6])

# Combine the extracted columns with the first and 6th columns from Cell1
Cell_all_r <- cbind(Cell1[, c(1, 6)], do.call(cbind, columns_r))
Cell_all_r$mean_val <- NULL

Cell_all_long <- gather(data = Cell_all_r, key="Cell", value="Value", 2:col_num) # se and sd calculated from the raw mean

#to get sd ans se data must be in long format
tf <- group_by(Cell_all_long, Angle) %>%
  summarise(w=mean(Value), sd = sd(Value)) %>% 
  arrange(desc(w))

se <- group_by(Cell_all_long, Angle) %>%
  summarise(w=mean(Value), se = std.error(Value)) %>% 
  arrange(desc(w))


Cell_all <- cbind(Cell_all, tf[, 3], se[, 3]) #combine se and sd for the raw mean with percentage mean
y_placement <- c(rep.int(150, 63))
y_placement <- as.data.frame(y_placement)
Cell_all <- cbind(Cell_all, y_placement)

#################################################################################################################
#write.csv(Cell_all,file="./Cell_all_POLAR-PL1-polar.csv",na='')
#Cell_all <- read.csv("Cell_all_POLAR-PL1-polar.csv")
####################################

ggplot(Cell_all, aes(x = Angle, y = y_placement)) +
  geom_segment(aes(x = 315, xend = 315, y = 150, yend = 250), linewidth = 0.1, color = "lightgray")+
  geom_segment(aes(x = 135, xend = 135, y = 150, yend = 250), linewidth = 0.1, color = "lightgray")+
  
  geom_point(aes(fill = median_all), shape = 21, color = "black", size = 6, stroke = 0.5, na.rm = TRUE) +
  geom_point(aes(x = 0, y = 0), colour = "lightgray", shape = 4, size = 6) +
  theme_classic() +
  labs(title = "XXXX", fill = "Fluorescence intensity, %") +
  ylim(0, 250) +
  scale_fill_paletteer_c("grDevices::GnBu", 100,
    direction = 1,
    limits = c(0, 100),
    breaks = c(0, 50, 100)
  ) +
  scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270)) +
  theme(
    legend.position = "top",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_polar(theta = "x", start = -pi / 2) +
  geom_line(aes(x = Angle, y = 150 + sd), color = "black", linewidth = 0.75) +
  geom_point(aes(x=315, y=150), colour="darkgray", fill = "lightgray", shape=22,size=1) +
  geom_point(aes(x=315, y=200), colour="lightgray", fill = "lightgray", shape=20,size=3) +
  geom_point(aes(x=315, y=250), colour="lightgray", fill = "lightgray", shape=20,size=3) +
  geom_text(aes(x = 318, y = 250, label = "100"), color = "darkgray", size = 2, angle = 45, vjust = 0) +
  geom_text(aes(x = 318, y = 200, label = "50"), color = "darkgray", size = 2, angle = 45, vjust = 0) +
  geom_point(aes(x=135, y=150), colour="darkgray", fill = "lightgray", shape=22,size=1) +
  geom_point(aes(x=135, y=200), colour="lightgray", fill = "lightgray", shape=20,size=3) +
  geom_point(aes(x=135, y=250), colour="lightgray", fill = "lightgray", shape=20,size=3) 

