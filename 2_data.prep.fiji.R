#prepare fiji generated data (2D) for syanpse count analysis

#install.packages("tidyverse")
library(tidyverse)


gen_start_log<-function(name="Debug"){
  dir.create("log", recursive=TRUE, showWarnings=FALSE) # make the folders before saving data
  local_log_file<-paste(format(Sys.time(),"%Y%m%d%H%M%S"), name, "log", sep=".")
  return(local_log_file)
}

gen_log<-function(text, log_file, datetime=FALSE){
  dir.create(file.path(".","log"), recursive=TRUE, showWarnings=FALSE)
  if (datetime) {
    write(paste(Sys.time(), text, sep=": "), file=file.path(".","log",log_file), append=TRUE)
    print(paste(Sys.time(), text, sep=": "))
  }
  else {
    write(text, file=file.path(".","log",log_file), append=TRUE)
    print(text)
  }
}

load.data <- function(data){
  if (file.exists(data)){
    suppressMessages(d <- read_csv(data, col_names=TRUE))
    name <- gsub(pattern=".csv", replacement="", data, ignore.case = TRUE)
    d$group <- name  #add group
    gen_log(paste("preparing data:", name), logfile)  
    return(d)
  }else{
    print(paste(data, "not found!!!"))
  }
}

#log results in file
logfile <- gen_start_log("Fiji.data.prep")
gen_log("prepare Fiji data", logfile, datetime=TRUE)

#delete pdf to start with clean one
if (file.exists("data2D.Rda")){
  gen_log("Deleted old data2D.Rda", logfile)
  unlink("data2D.Rda", recursive = TRUE)
}

#files to read
files <- c("DAPI.csv", 
           "Post.csv",
           "Pre.csv"
           )

#read files and combined then in to one tibble
data <- list()
for (i in seq_along(files)){
  if (file.exists(files[i])){
    data[[i]] <- load.data(files[i])
  } else {
      gen_log(paste0("File ", files[i], " is missing!!!"), logfile)
  }
}
data <- bind_rows(data)

colnames(data)[1] <- "ROI"            #fix first columns name
data$group <- factor(data$group)      #convert group to factor
#Perim. causes problems when splitting using . %Area is problematic too
colnames(data)[which(colnames(data) == "Perim.")] <- "Perimeter"
colnames(data)[which(colnames(data) == "%Area")] <- "PercentArea"
colnames(data)[which(colnames(data) == "Circ.")] <- "Circularity"
print(data)

save(data, file=paste0("data2D.Rda"))
gen_log("Fiji data collected in data2D.Rda", logfile, datetime=TRUE)
