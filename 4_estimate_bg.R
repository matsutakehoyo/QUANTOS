#2017Nov30	Improved bg estimation using image histogram (all pixels)
library(grid)					  #grid_arrange_shared_legend & vplayout
library(tidyverse)                #ggplot, tibble, tidyr, readr, purrr, dplyr
library("wesanderson")


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



#graph styles
basic_style = theme(plot.title = element_text(lineheight=.8, face="bold"),
				plot.background = element_rect(fill = "transparent",colour = NA),
				legend.background= element_rect(fill = "transparent",colour = NA),
				axis.title=element_text(face="bold"))
				

keynote = theme(plot.title = element_text(size=18),
                axis.title=element_text(size=14,face="bold"), 
                axis.title.x = element_text(size=16,face="bold"), 
                axis.title.y = element_text(size=16,face="bold"))

dark_background =  theme(plot.title = element_text(color="white"),
                axis.title=element_text(color="white"), 
                axis.title.x = element_text(color="white"), 
                axis.title.y = element_text(color="white"))

transparent = theme(
    panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    legend.background= element_rect(fill = "transparent",colour = NA))


logfile <- gen_start_log("Estimate Background intensity")
gen_log("Estimate image background intensity and normaize signals", logfile, datetime=TRUE)

background_th = 8  #igrnore values smaller than this. (Changed from 2 to 8 on 2018/3/30)
gen_log(paste("Background intensity threshold:", background_th), logfile)
#load data
data <- readRDS("data.rds")
data_pair <- readRDS("data_pair.rds")

#load text image
mat_dapi <- read.table("DAPIZproject.txt") %>% data.matrix()
mat_post <- read.table("PostZproject.txt") %>% data.matrix()
mat_pre <- read.table("PreZproject.txt") %>% data.matrix()

intensity <- tibble(group=rep("DAPI", dim(mat_dapi)[1] * dim(mat_dapi)[2]), value=as.vector(mat_dapi))
intensity <- intensity %>% bind_rows(., tibble(group=rep("Pre", dim(mat_pre)[1] * dim(mat_pre)[2]), 
											   value=as.vector(mat_pre)))
intensity <- intensity %>% bind_rows(., tibble(group=rep("Post", dim(mat_post)[1] * dim(mat_post)[2]),
											   value=as.vector(mat_post)))

background <- intensity %>% 
  filter(value > background_th) %>% #avoid completely empty pixels
  dplyr::group_by(group) %>%
  dplyr::summarise(intensity_median = median(value),
                   intensity_mean = mean(value), 
                   intensity_sd = sd(value)) %>%
  gather(center, value, -group) %>%
  mutate(center=gsub(pattern="intensity_", replace="", center))

gen_log(capture.output(print(background)), logfile)

fibonacci <- function(n) {
	if(n <= 1) {
		as.integer(n)
	} else {
		as.integer(fibonacci(n-1) + fibonacci(n-2))
	}
}

histo <- list()    #data to meke histogram of different bins
dens <- list()	   #estimate denity	
peaks <- list()	   #peak(max) values of the density
normal <- list()
i <- 1
for (group_str in c("Pre", "Post", "DAPI")){
	gen_log(paste("Calculating", group_str, " peak"), logfile, datetime=TRUE)
	for (n in 2:9){
		bin <- fibonacci(n)
		h <- intensity%>% filter(group==group_str) %>% 
			dplyr::select(value) %>% pull() %>%
			hist(breaks=seq(0, 255+bin, bin), plot = FALSE)
		histo[[i]] <- tibble(bins = bin, counts = h$counts, density = h$density, mids = h$mids, group = group_str)
		kde <- density(intensity %>% filter(group==group_str & value > background_th) %>% pull(), adjust=bin*2, from=0, to=255, n=2550)
		dens_kde <- approxfun(kde)
		dens[[i]] <- tibble(bins = bin, intensity = seq(0,255, 0.1), density = dens_kde(seq(0,255, 0.1)), group = group_str)
		peak_index <- which.max(dens_kde(seq(0, 255, 0.1)))
		peak_x <- seq(0, 255, 0.1)[peak_index]
		peak_y <- dens_kde(seq(0, 255, 0.1)[peak_index])
		peaks[[i]] <- tibble(bins = bin, peak = peak_x, value = peak_y, group = group_str)
		i <- i + 1
	}
	normal[[which(c("Pre", "Post", "DAPI") %in% group_str)]] <-
		tibble(intensity = seq(0,255),
			   value = dnorm(x = seq(0,255), 
			   				 mean = background %>% filter(center=="mean" & group==group_str) %>% dplyr::select(value) %>% pull(), 
							 sd = background %>% filter(center=="sd" & group==group_str) %>% dplyr::select(value) %>% pull()),
			   group = group_str)

}
histo <- bind_rows(histo)
dens <- bind_rows(dens)
peaks <- bind_rows(peaks)
normal <- bind_rows(normal)
gen_log(paste("Peak estimation done"), logfile, datetime=TRUE)
gen_log(capture.output(print(peaks, n=Inf)), logfile)


h1 <- ggplot(histo, aes(x = mids, y = density, fill=as.factor(bins))) + 
	geom_bar(stat = "identity", aes(width = bins)) +
	scale_fill_manual(values = wes_palette(names(wes_palettes)[10], 8, type = "continuous")) +
	geom_line(data = dens, aes(x = intensity, y = density), color = "gray50") +
	# geom_vline(data = peaks, aes(xintercept = peak), color = "gray40") +
	geom_point(data = peaks, aes(x = peak, y = value), color="red") +
	geom_vline(data = background %>% filter(center=="median"), 
			   aes(xintercept = value), color = "gray60") +
	geom_vline(data = background %>% filter(center=="mean"), 
			   aes(xintercept = value), color = "gray70") +
	geom_text(data = peaks, aes(x = peak, y = value, label = peak), 
			  hjust = "left", vjust = "bottom") +
	facet_grid(group~bins, scales="free_y", labeller=label_both) +
	xlab("intensity") + ylab("density") +
	ggtitle(paste0(group_str, " Intensity Histogram")) +
	basic_style +
	theme(legend.position="none")
ggsave("Intensity Histogram.pdf", h1, width = 15, height = 10)
	
h2 <- ggplot(histo %>% filter(bins==1), aes(x = mids, y = density, fill=as.factor(group))) + 
	geom_bar(stat = "identity", aes(width = bins)) +
	#scale_fill_manual(values = wes_palette("Zissou")) +
	geom_vline(data = background %>% filter(center=="median"), 
			   aes(xintercept = value), color = "gray60") +
	geom_vline(data = background %>% filter(center=="mean"), 
			   aes(xintercept = value), color = "gray70") +
	geom_line(data = normal, aes(x = intensity, y = value), color="gray50") +
	facet_grid(group~., scales="free_y", labeller=label_both) +
	xlab("intensity") + ylab("density") +
	ggtitle(paste0(group_str, " Intensity Histogram (fit Normal Curve)")) +
	basic_style +
	theme(legend.position="none")
ggsave("Intensity Histogram (normal).pdf", h2, width = 15, height = 10)

background_pre <- peaks %>% filter(group=="Pre" & bins==5) %>% dplyr::select(peak) %>% pull()
background_post <- peaks %>% filter(group=="Post" & bins==5) %>% dplyr::select(peak) %>% pull()
background_DAPI <- peaks %>% filter(group=="DAPI" & bins==5) %>% dplyr::select(peak) %>% pull()

gen_log(paste("Pre intensity bg:", background_pre), logfile)
gen_log(paste("Post intensity bg:", background_post), logfile)
gen_log(paste("DAPI intensity bg:", background_DAPI), logfile)


#remove markers with low signal
#data_bg <- data %>% filter(!(group=="Pre" & Mean < background_pre * bg_level_pre))
#data_bg <- data_bg %>% filter(!(group=="Post" & Mean < background_post * bg_level_post))


#normalize marker signal (mean) with background
normalize_data <- function(df, group_str){
	df <- filter(df, group==group_str)
	#bg <- peaks %>% filter(group==group_str & bins==5) %>% dplyr::select(peak) %>% pull()
	bg <- peaks %>% filter(group==group_str & bins==8) %>% dplyr::select(peak) %>% pull()
	df$bg <- bg
	df$Mean_bg <- df$Mean / bg
	df$Mode_bg <- df$Mode / bg
	df$Median_bg <- df$Mean / bg
	df$Min_bg <- df$Min / bg
	df$Max_bg <- df$Max / bg
	df$StdDev_bg <- df$StdDev / bg
	df$IntDen_bg <- df$IntDen / bg
	df$RawIntDen_bg <- df$RawIntDen / bg
	#bind_rows(filter(df, group!=group_str), df)
	df
}

if (length(data_pair)){
  data_pair_bg <- bind_rows(normalize_data(data_pair, "Pre"), 
	              					  normalize_data(data_pair, "Post"))
} else {
  data_pair_bg <-data_pair
}
data_bg <- bind_rows(normalize_data(data, "Pre"), 
					 normalize_data(data, "Post"), 
					 normalize_data(data, "DAPI"))

#normalize marker signal (mean) with background
standardise_data <- function(df, group_str){
	df <- filter(df, group==group_str)
	sd <- background %>% filter(group==group_str & center=="sd") %>% dplyr::select(value) %>% pull()
	mu <- background %>% filter(group==group_str & center=="mean") %>% dplyr::select(value) %>% pull()
	df$bg_sd <- sd
	df$bg_mu <- mu
	df$Mean_z <- df$Mean - mu / sd
	df$Mode_z <- df$Mode - mu / sd
	df$Median_z <- df$Mean - mu / sd
	df$Min_z <- df$Min - mu / sd
	df$Max_z <- df$Max - mu / sd
	df$StdDev_z <- df$StdDev - mu / sd
	df$IntDen_z <- df$IntDen - mu / sd
	df$RawIntDen_z <- df$RawIntDen - mu / sd
	#bind_rows(filter(df, group!=group_str), df)
	df
}

if (length(data_pair)){
  data_pair_bg <- bind_rows(standardise_data(data_pair_bg, "Pre"), 
						  standardise_data(data_pair_bg, "Post"))
}

data_bg <- bind_rows(standardise_data(data_bg, "Pre"), 
					 standardise_data(data_bg, "Post"), 
					 standardise_data(data_bg, "DAPI"))

if (file.exists("data_pair_bg.rds")){
  gen_log("Deleted old data_pair_bg.rds", logfile)
  unlink("data_pair_bg.rds", recursive = TRUE)
}
saveRDS(data_pair_bg, file = "data_pair_bg.rds")        

if (file.exists("data_bg.rds")){
  gen_log("Deleted old data_bg.rds", logfile)
  unlink("data_bg.rds", recursive = TRUE)
}
saveRDS(data_bg, file = "data_bg.rds")          
gen_log("data tibble created: data.rds", logfile)

gen_log("Done", logfile, datetime=TRUE)