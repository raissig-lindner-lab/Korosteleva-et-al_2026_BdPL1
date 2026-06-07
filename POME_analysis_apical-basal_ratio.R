install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("multcompView")
install.packages("multcomp")
install.packages("qpcR")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(multcompView)
library(multcomp)
library(qpcR)

print(getwd()) #check that your working directory corresponds to the location of the source script
#install.packages("rstudioapi") # for RStudio users
#library(rstudioapi) # for RStudio users
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to the file location


rm(list = ls()) #remove everything from the global environment

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

  Cell
}
##################################################################################################

files <- list.files(path = "./data/01_PL1-PL1", pattern = ".csv", full.names = TRUE)
l <- length(files)

for(i in 1:9){
  x <- intensity(files[i])
  assign(paste0("Cell0",i),x,.GlobalEnv)
}
for(i in 10:l){
  x <- intensity(files[i])
  assign(paste0("Cell",i),x,.GlobalEnv)
}

# Dynamically create a list of all data frames whose names start with "Cell"
dfs <- mget(ls(pattern = "^(Cell)\\d*"))

columns <- lapply(dfs, function(df) df[, 6])

# Combine the extracted columns for individual cells
Cell_all <- cbind(Cell01[, c(1, 6)], do.call(cbind, columns))
Cell_all$mean_val <- NULL

Cell_ratio <- (Cell_all[16, c(2:21)]/Cell_all[47, c(2:21)])

PL1_PL1_ratio_long <- gather(data = Cell_ratio, key="Cell", value="Ratio", 1:20) 

PL1_PL1_ratio_long <- PL1_PL1_ratio_long[ ,2]
PL1_PL1_ratio_long <- as.data.frame(PL1_PL1_ratio_long)

rm(list = ls(pattern = "^Cell")) #remove data frames for individual cells

###################################################################################################
###################################################################################################

files <- list.files(path = "./data/02_POLAR-POLAR", pattern = ".csv", full.names = TRUE)
l <- length(files)

for(i in 1:9){
  x <- intensity(files[i])
  assign(paste0("Cell0",i),x,.GlobalEnv)
}
for(i in 10:l){
  x <- intensity(files[i])
  assign(paste0("Cell",i),x,.GlobalEnv)
}

dfs <- mget(ls(pattern = "^(Cell)\\d*"))

columns <- lapply(dfs, function(df) df[, 6])

# Combine the extracted columns 
Cell_all <- cbind(Cell01[, c(1, 6)], do.call(cbind, columns))
Cell_all$mean_val <- NULL

Cell_ratio <- (Cell_all[16, c(2:31)]/Cell_all[47, c(2:31)])

POLAR_POLAR_ratio_long <- gather(data = Cell_ratio, key="Cell", value="Ratio", 1:30) 

POLAR_POLAR_ratio_long <- POLAR_POLAR_ratio_long[ ,2]
POLAR_POLAR_ratio_long <- as.data.frame(POLAR_POLAR_ratio_long)

rm(list = ls(pattern = "^Cell"))

###################################################################################################
###################################################################################################

files <- list.files(path = "./data/05_PL1-POLAR-pl1", pattern = ".csv", full.names = TRUE)
l <- length(files)

for(i in 1:9){
  x <- intensity(files[i])
  assign(paste0("Cell0",i),x,.GlobalEnv)
}
for(i in 10:l){
  x <- intensity(files[i])
  assign(paste0("Cell",i),x,.GlobalEnv)
}

dfs <- mget(ls(pattern = "^(Cell)\\d*"))

columns <- lapply(dfs, function(df) df[, 6])

# Combine the extracted columns 
Cell_all <- cbind(Cell01[, c(1, 6)], do.call(cbind, columns))
Cell_all$mean_val <- NULL

Cell_ratio <- (Cell_all[16, c(2:24)]/Cell_all[47, c(2:24)])

PL1_POLAR_pl1_ratio_long <- gather(data = Cell_ratio, key="Cell", value="Ratio", 1:23) 

PL1_POLAR_pl1_ratio_long <- PL1_POLAR_pl1_ratio_long[ ,2]
PL1_POLAR_pl1_ratio_long <- as.data.frame(PL1_POLAR_pl1_ratio_long)

rm(list = ls(pattern = "^Cell"))

###################################################################################################
###################################################################################################

files <- list.files(path = "./data/03_PL1-POLAR-wt", pattern = ".csv", full.names = TRUE)
l <- length(files)

for(i in 1:9){
  x <- intensity(files[i])
  assign(paste0("Cell0",i),x,.GlobalEnv)
}
for(i in 10:l){
  x <- intensity(files[i])
  assign(paste0("Cell",i),x,.GlobalEnv)
}

dfs <- mget(ls(pattern = "^(Cell)\\d*"))

columns <- lapply(dfs, function(df) df[, 6])

# Combine the extracted columns 
Cell_all <- cbind(Cell01[, c(1, 6)], do.call(cbind, columns))
Cell_all$mean_val <- NULL

Cell_ratio <- (Cell_all[16, c(2:21)]/Cell_all[47, c(2:21)])

PL1_POLAR_wt_ratio_long <- gather(data = Cell_ratio, key="Cell", value="Ratio", 1:20) 

PL1_POLAR_wt_ratio_long <- PL1_POLAR_wt_ratio_long[ ,2]
PL1_POLAR_wt_ratio_long <- as.data.frame(PL1_POLAR_wt_ratio_long)

rm(list = ls(pattern = "^Cell"))

###################################################################################################
###################################################################################################

files <- list.files(path = "./data/06_POLAR-PL1-polar", pattern = ".csv", full.names = TRUE)
l <- length(files)

for(i in 1:9){
  x <- intensity(files[i])
  assign(paste0("Cell0",i),x,.GlobalEnv)
}
for(i in 10:l){
  x <- intensity(files[i])
  assign(paste0("Cell",i),x,.GlobalEnv)
}

dfs <- mget(ls(pattern = "^(Cell)\\d*"))

columns <- lapply(dfs, function(df) df[, 6])

# Combine the extracted columns 
Cell_all <- cbind(Cell01[, c(1, 6)], do.call(cbind, columns))
Cell_all$mean_val <- NULL

Cell_ratio <- (Cell_all[16, c(2:19)]/Cell_all[47, c(2:19)])

POLAR_PL1_polar_ratio_long <- gather(data = Cell_ratio, key="Cell", value="Ratio", 1:18)

POLAR_PL1_polar_ratio_long <- POLAR_PL1_polar_ratio_long[ ,2]
POLAR_PL1_polar_ratio_long <- as.data.frame(POLAR_PL1_polar_ratio_long)

rm(list = ls(pattern = "^Cell"))

###################################################################################################
###################################################################################################

files <- list.files(path = "./data/04_POLAR-PL1-wt", pattern = ".csv", full.names = TRUE)
l <- length(files)

for(i in 1:9){
  x <- intensity(files[i])
  assign(paste0("Cell0",i),x,.GlobalEnv)
}
for(i in 10:l){
  x <- intensity(files[i])
  assign(paste0("Cell",i),x,.GlobalEnv)
}

dfs <- mget(ls(pattern = "^(Cell)\\d*"))

columns <- lapply(dfs, function(df) df[, 6])

# Combine the extracted columns 
Cell_all <- cbind(Cell01[, c(1, 6)], do.call(cbind, columns))
Cell_all$mean_val <- NULL

Cell_ratio <- (Cell_all[16, c(2:21)]/Cell_all[47, c(2:21)])

POLAR_PL1_wt_ratio_long <- gather(data = Cell_ratio, key="Cell", value="Ratio", 1:20)

POLAR_PL1_wt_ratio_long <- POLAR_PL1_wt_ratio_long[ ,2]
POLAR_PL1_wt_ratio_long <- as.data.frame(POLAR_PL1_wt_ratio_long)

rm(list = ls(pattern = "^Cell"))

#################################################################################################


Ratio_all <- qpcR:::cbind.na(PL1_PL1_ratio_long, PL1_POLAR_pl1_ratio_long, PL1_POLAR_wt_ratio_long,
                   POLAR_POLAR_ratio_long, POLAR_PL1_polar_ratio_long, POLAR_PL1_wt_ratio_long) #combine all extracted columns

Ratio_all_long <- gather(data = Ratio_all, key = line, value = "Value", , na.rm=TRUE) #convert into long format

#write.csv(Ratio_all_long,"./ratio_all_long.csv", row.names = T) #save gathered data
#Ratio_all_long <- read.csv("ratio_all_long.csv",header = TRUE, sep = ",", check.names = F) #reload gathered data

########## statistics
anova <- aov(Value ~ line, Ratio_all_long)
tukey <- TukeyHSD(anova)
print(tukey)
# compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# table with factors and 3rd quantile
dt <- group_by(Ratio_all_long, line) %>%
  summarise(w=mean(Value), sd = sd(Value)) %>%
  arrange(desc(w))

# extracting the compact letter display and adding to the dt table
cld <- as.data.frame.list(cld$line)
dt$cld <- cld$Letters
print(dt)

#plot
ggplot(Ratio_all_long,aes(x = line, y = Value ))+
  geom_boxplot(aes(fill = Value), fill = "#A2DBB5", color = "black") +
  geom_jitter(width=0.1, height = 0, color = "black", show.legend = F, size = 1.5, alpha = 0.5)+
  scale_y_continuous(limits=c(0,3)) +
  labs(y = "Ratio value, apical/basal", x = "Lines")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position="bottom")+
  geom_text(data = dt, check_overlap = TRUE, 
            aes(x = line, y = 2.05, label = cld), size = 10, hjust =0.5) +
  theme( text=element_text(family="sans"),
         axis.text.x = element_text(color="black", size=15, angle=15, margin = margin(t = 12)),              
         axis.title.x = element_text(color="dimgrey", size=14, margin = margin(t = -10)),
         axis.title.y = element_text(color="dimgrey", size=14, margin = margin(r = 5)))

