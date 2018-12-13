#summarize synapse count
#2017Sep28 Estimate pdf with optimum distribution propagate:fitDistr
#2017Oct09 Parametric estimation results in too many NaNs, values out of 
#          estimated densities > implemented non parametric density estimation
#2017Oct22 Remove unecessary libraries
#          kde estimation with boundaries > renormalize
#          density is not calculated for noise data, need to calculate local density
#2018Jan03 Collect data from Synapse and Noise to estimate pdf together


library(MASS)                    #for for fitdistr()  
library(RColorBrewer)
# library(GGally)                  #for ggpairs()
# library(PerformanceAnalytics)
library(gridExtra)               #for viewport
library(grid)                    #for grid.newpage
library(ks)                      #kde
library(bde)                     #bde
library(tidyverse)


#Script ot log R output ot file: moified original file by Sunagawa gen-log.R
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
	} else {
		write(text, file=file.path(".","log",log_file), append=TRUE)
		print(text)
	}
}

#log results in file
logfile <- gen_start_log("SynapseCount")
gen_log("Ideel Marker pdf", logfile, datetime=TRUE)

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

## Convert x-ticks to fractional x-ticks with a symbol multiplier
fracAx <- function(p, symbol, width=0.5) {
	val <- tryCatch(eval(parse(text=symbol)), error=function(e) 1)
	info <- ggplot_build(p)
	xrange <- info[[2]]$panel_ranges[[1]]$x.range/val             # get the x-range of figure
	vec.breaks <- seq(floor(xrange[1]), ceiling(xrange[2]), by=width)
	fracs <- strsplit(attr(fractions(vec.breaks), "fracs"), "/")  # convert to fractions
	labels <- sapply(fracs, function(i)
	if (length(i) > 1) { paste(i[1], "*", symbol, "/", i[2]) }
	else { paste(i, "*", symbol) })
	invisible(p <- p + scale_x_continuous(breaks=vec.breaks*val, labels=parse(text=labels)))
}

#for multiple plots in pages
vplayout <- function(x,y) viewport(layout.pos.row = x, layout.pos.col =y)

#plots in a grid with shared leyend
grid_arrange_shared_legend <- function(plots, ncol = length(plots), nrow = 1, position = c("bottom", "right")) {  
	position <- match.arg(position)
	g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
	legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
	lheight <- sum(legend$height)
	lwidth <- sum(legend$width)
	gl <- lapply(plots, function(x) x + theme(legend.position="none"))
	gl <- c(gl, ncol = ncol, nrow = nrow)
	combined <- switch(position,
					 "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
											legend,
											ncol = 1,
											heights = unit.c(unit(1, "npc") - lheight, lheight)),
					 "right" = arrangeGrob(do.call(arrangeGrob, gl),
										   legend,
										   ncol = 2,
										   widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
	grid.newpage()
	grid.draw(combined)
	invisible(combined)                           # return gtable invisibly
}

#local density calculation
calc_density <- function(x, y, df, distance){    #df is tibble with X & Y coordinates
	df %>% filter((X > x - distance/2 & X < x + distance/2) & (Y > y - distance/2 & Y < y + distance/2)) %>%
		nrow() -> count
	(count)/(distance^2)
}


#read rds fiels with the name filename and combine them"
gen_log("Loading data...", logfile, datetime=TRUE)
files_synapse <- list.files(path = "./IdealSynapse", pattern = paste0("data_pair_bg.rds"), recursive=TRUE)
files_noise <- list.files(path = "./IdealNoise", pattern = paste0("data_bg.rds"), recursive=TRUE)
d <- lst()
for (i in seq_along(files_synapse)){
	d[[i]] <- readRDS(paste0("./IdealSynapse/", files_synapse[i])) %>% 
		filter(0.3 < pair_dist & pair_dist < 0.8) %>% #use crème de la crème
		mutate(summary = "synapse") 
	gen_log(paste0(files_synapse[i], "  >>  ", dim(d[[i]])[1], " markers loaded"), logfile)
}	
for (i in seq_along(files_noise) + length(files_synapse)){
	d[[i]] <- readRDS(paste0("./IdealNoise/", files_noise[i - length(files_synapse)])) %>% 
		filter(group!="DAPI") %>% 
		mutate(shape = sqrt((XM-(BX+Width/2))^2+(YM-(BY+Height/2))^2)) %>%
		mutate(summary = "noise")
	d[[i]]$density <- vector("double", nrow(d[[i]]))  
	for (r in 1:nrow(d[[i]])){
		d[[i]]$density[r] <- calc_density(d[[i]]$X[r], d[[i]]$Y[r], filter(d[[i]], group==as.character(d[[i]]$group[r])), 4)
	}
	gen_log(paste0(files_noise[i - length(files_synapse)], "  >>  ", dim(d[[i]])[1], " markers loaded"), logfile)
}	

data <- bind_rows(d) 
gen_log(capture.output(data %>% glimpse()), logfile)	
gen_log("Data loaded", logfile, datetime=TRUE)

if (file.exists("data_summary.rds")){
	gen_log("Deleted old data_summary.rds", logfile)
	unlink("data_summary.rds", recursive = TRUE)
}
saveRDS(data, file = "data_summary.rds")          
gen_log("data tibble created: data_summary.rds", logfile)


#all the features   
features <- colnames(data[sapply(data,is.double)])
#remove some unecessary features
features <- features[!features %in% c("X", "Y", "XM","YM", "BX", "BY", "bg", 
									  "X.Area", "MinThr", "MaxThr", "PercentArea", 
									  "FeretX", "FeretY")]

#delete pdf to start with clean one
if (file.exists("pdf")){
	unlink("pdf", recursive = TRUE)
}
dir.create("pdf", showWarnings=FALSE) # make the folders before saving data

#distance pdf
bw_dist <- data %>% 
		   filter(summary=="synapse") %>% 
		   dplyr::select(pair_dist) %>% 
		   pull() %>%
		   hpi()
h_pair_dist_row <- ggplot(data %>% filter(summary=="synapse"), aes(x=pair_dist)) + 
	geom_histogram(aes(y=..density..), binwidth = bw_dist, position="identity") + 
	basic_style

dist.pdf <- fitdistr(filter(data %>% filter(summary=="synapse"), pair_dist<0.8)$pair_dist,"normal") 
h_pair_dist <- h_pair_dist_row + 
	#use normal distribution for synapse
	stat_function(fun = dnorm, n = 200, 
			args = list(mean = dist.pdf$estimate[1], sd=dist.pdf$estimate[2]), color="darkgray", size=1) +
	stat_function(fun = function(x){2/3*x^2}, color="red", size=1)

saveRDS(dist.pdf, file=(paste0("pdf/dist.pdf.synapse.rds")))

#angle pdf
bw_angle <- data %>% 
			filter(summary=="synapse" & pair_dist < 0.8) %>% 
			dplyr::select(pair_angle) %>% 
			pull() %>%
			hpi()
h_pair_angle_row <- ggplot(data %>% filter(summary=="synapse" & pair_dist < 0.8), aes(x=pair_angle)) + 
	geom_histogram(aes(y=..density..), binwidth=bw_angle, position="identity")  + 
	basic_style
h_pair_angle <- h_pair_angle_row + 
	#use normal distribution for synapse
	stat_function(fun = dnorm, n = 200, args = list(mean = 0, sd=pi/4), color="darkgray", size=1)+
	stat_function(fun = propagate:::dlaplace, n = 200, args = list(mean = 0, sigma=pi/2), color="red", size=1)
h_pair_angle <- fracAx(h_pair_angle, "pi")    

#plot the pre pdf 
stat_function_kde <-function(pdf){
  stat_function(fun = dkde, n = 200, args = list(fhat = pdf), aes(color="kde"), alpha=1/2)
}

stat_function_bde <-function(pdf){
  stat_function(fun = bde::density, n = 200, args = list(x = pdf), aes(color="bde"), alpha=1/2)
}

robust_hpi <- function(data){
  tryCatch({
	  min <- data %>% na.omit() %>% pull() %>% min()
	  hpi <- hpi(data %>% na.omit() %>% pull())
	  hpi
	  }, warning = function(w){
		print(paste("Warning in hpi estimation"))
	  }, error = function(err){
		print(paste("Error in hpi estimation using min value insted:", min * 2))
		min
	  }
  )
}

#make pdf, plot and save > use bde to acount for boundaries
#adjust: higher number results in more smoothing
estimate_pdf <- function(feature, group_str, pdf, adjust_bde=1, adjust_kde=1, ...){
	print(paste("Estimating", group_str, feature, pdf, "density distribution"))
	bw <- data %>% 
			filter(group==group_str) %>% 
			dplyr::select(feature) %>% 
			robust_hpi()
	
	#calculate a reasonable range for posterior distributions excluding ontliers
	posterior_range <- function(posterior, feature, hdi=0.99){
	  posterior_range <- posterior %>%
	    summarise(hdi_lower = hdi(., credMass=hdi)[1], 
	              hdi_upper = hdi(., credMass=hdi)[2]) %>%
	    summarise(min = min(hdi_lower),
	              max = max(hdi_upper)) %>%
	    unlist()
	  return(posterior_range)
	}

	feature_max <- posterior_range(data %>% filter(group==group_str) %>% dplyr::select(feature), hdi=0.999)[2]*1.1
	feature_min <- posterior_range(data %>% filter(group==group_str) %>% dplyr::select(feature), hdi=0.999)[1]*0.9
	

	h_data <- ggplot() + 
		geom_histogram(data= data %>% filter(group==group_str), aes(y=..density.., x=eval(parse(text=paste0(feature))), fill=summary), binwidth = bw*4, alpha=1/2, position="identity") +
		scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
	  xlim(feature_min, feature_max)+
		xlab(paste(feature, "(", group_str, ")")) +
		basic_style + keynote

	bde <- data %>% 
			filter(group==group_str) %>% 
			dplyr::select(feature) %>% 
			na.omit() %>% pull() %>%
			bde(., estimator="vitale"
					, ...
				)

	bde_synapse <- data %>% 
			filter(summary=="synapse") %>% 
			filter(group==group_str) %>% 
			dplyr::select(feature) %>% 
			na.omit() %>% pull() %>%
			bde(., estimator="vitale"
					,b = 1/getm(bde) * adjust_bde
					, ... 	
				)

	bde_noise <- data %>% 
			filter(summary=="noise") %>% 
			filter(group==group_str) %>% 
			dplyr::select(feature) %>% 
			na.omit() %>% pull() %>%
			bde(., estimator="vitale"
					,b = 1/getm(bde) * adjust_bde
					, ...					
				)

	#bw for kde 
	bw <- bw * adjust_kde 
	kde_synapse <- data %>% 
		filter(summary=="synapse") %>% 
		filter(group==group_str) %>% 
		dplyr::select(feature) %>% 
		na.omit() %>% pull() %>%
		kde(h=bw, 
			xmin=min(bde_synapse@lower.limit, bde_noise@lower.limit), 
			xmax=min(bde_synapse@upper.limit, bde_noise@upper.limit), 
			)

	kde_noise <- data %>% 
		filter(summary=="noise") %>% 
		filter(group==group_str) %>% 
		dplyr::select(feature) %>% 
		na.omit() %>% pull() %>%
		kde(h=bw, 
			xmin=min(bde_synapse@lower.limit, bde_noise@lower.limit), 
			xmax=min(bde_synapse@upper.limit, bde_noise@upper.limit), 
			)

	pdf_points <- seq(min(bde_synapse@lower.limit, bde_noise@lower.limit), max(bde_noise@upper.limit, bde_synapse@upper.limit), length.out=200)
	pdf_data <- tibble(y=bde::density(bde_synapse, pdf_points), x=pdf_points, estimate="synapse", pdf="bde") %>%
				bind_rows(tibble(y=bde::density(bde_noise, pdf_points), x=pdf_points, estimate="noise", pdf="bde")) %>%
				bind_rows(tibble(y=dkde(fhat=kde_synapse, x=pdf_points), x=pdf_points, estimate="synapse", pdf="kde")) %>%
				bind_rows(tibble(y=dkde(fhat=kde_noise, x=pdf_points), x=pdf_points, estimate="noise", pdf="kde"))
	
	
	h <- h_data + 
		geom_line(data=pdf_data, mapping=aes(x=x, y=y, colour=estimate, linetype=pdf)) +
	  xlim(feature_min, feature_max)+
		scale_color_manual(values=c("#E69F00", "#56B4E9")) 

	if (pdf=="bde") {
		assign(paste0(feature, ".", group_str, ".bde.synapse"), bde_synapse, envir=globalenv())
		assign(paste0(feature, ".", group_str, ".bde.noise"), bde_noise, envir=globalenv())
		saveRDS(bde_synapse, file = paste0("pdf/", feature, ".", group_str, ".bde.synapse.rds"))
		saveRDS(bde_noise, file = paste0("pdf/", feature, ".", group_str, ".bde.noise.rds"))
		h <- h + 
			annotate("text", x=0, y=Inf, label = "bde", hjust = 0, vjust = 1)
	} else if (pdf=="kde") {
		assign(paste0(feature, ".", group_str, ".kde.synapse"), kde_synapse, envir=globalenv())
		assign(paste0(feature, ".", group_str, ".kde.noise"), kde_noise, envir=globalenv())
		saveRDS(kde_synapse, file = paste0("pdf/", feature, ".", group_str, ".kde.synapse.rds"))
		saveRDS(kde_noise, file = paste0("pdf/", feature, ".", group_str, ".kde.noise.rds"))
		h <- h + 
			annotate("text", x=0, y=Inf, label = "kde", hjust = 0, vjust = 1)
	} else {
		print("Unspecified pdf type!!!")
	}


	ggsave(paste0("pdf_", feature, "_", group_str, ".pdf"))
	assign(paste0("h_", feature, "_", group_str), h, envir=globalenv())  
}

gen_log("estimating pdf...", logfile, datetime=TRUE)

features_signal <- features[features %in% c("Mean", "Median", "Mode", "Max", "Min")]
for (feature in features_signal){
  estimate_pdf(feature, "Pre", "bde", adjust_bde=3, adjust_kde=10, lower.limit=0, upper.limit=255)
  estimate_pdf(feature, "Post", "bde", adjust_bde=3, adjust_kde=10, lower.limit=0, upper.limit=255)
}

features_signal_bg <- features[features %in% c("Mean_bg", "Median_bg", "Mode_bg", "Max_bg", "Min_bg")]
for (feature in features_signal_bg){
  bg_Pre <- data %>% filter(group=="Pre") %>% dplyr::select(bg) %>% unique %>% pull %>% mean
  estimate_pdf(feature, "Pre", "bde", adjust_bde=3, adjust_kde=10,  lower.limit=0, upper.limit=255/bg_Pre)
  bg_Post <- data %>% filter(group=="Post") %>% dplyr::select(bg) %>% unique %>% pull %>% mean
  estimate_pdf(feature, "Post", "bde", adjust_bde=3, adjust_kde=10, lower.limit=0, upper.limit=255/bg_Post)
}

features_signal_z <- features[features %in% c("Mean_z", "Median_z", "Mode_z", "Max_z", "Min_z")]
for (feature in features_signal_z){
  min_Pre <- data %>% filter(group=="Pre") %>% dplyr::select(feature) %>% min()
  estimate_pdf(feature, "Pre", "bde", adjust_bde=3, adjust_kde=10, lower.limit = min_Pre)
  min_Post <- data %>% filter(group=="Post") %>% dplyr::select(feature) %>% min()
  estimate_pdf(feature, "Post", "bde", adjust_bde=3, adjust_kde=10, lower.limit = min_Post)
}


estimate_pdf("Area", "Pre", "bde", adjust_bde=1/2, adjust_kde=4, lower.limit=0, upper.limit=5)
estimate_pdf("Area", "Post", "bde", adjust_bde=1/2, adjust_kde=4, lower.limit=0, upper.limit=5)

features_zero1 <- features[features %in% c("StdDev", "StdDev_bg", "StdDev_z",
										   "IntDen", "IntDen_bg", "IntDen_z", 
										  "RawIntDen", "RawIntDen_bg", "RawIntDen_z")]
for (feature in features_zero1){
  min_Pre <- data %>% filter(group=="Pre") %>% dplyr::select(feature) %>% min()
  estimate_pdf(feature, "Pre", "kde", adjust_bde=1, adjust_kde=8, lower.limit = min_Pre)
  min_Post <- data %>% filter(group=="Post") %>% dplyr::select(feature) %>% min()
  estimate_pdf(feature, "Post", "kde", adjust_bde=1, adjust_kde=8, lower.limit = min_Post)
}

estimate_pdf("density", "Pre", "bde", adjust_bde=1, adjust_kde=16, lower.limit=0)
estimate_pdf("density", "Post", "bde", adjust_bde=1, adjust_kde=16, lower.limit=0)

features_zero2 <- features[features %in% c("Minor", "Major", "Height", "Width", "Perimeter",
										   "shape", "Feret")]
for (feature in features_zero2){
  estimate_pdf(feature, "Pre", "kde", adjust_bde=1/2, adjust_kde=8, lower.limit=0)
  estimate_pdf(feature, "Post", "kde", adjust_bde=1/2, adjust_kde=8, lower.limit=0)
}

estimate_pdf("MinFeret", "Pre", "kde", adjust_bde=1/2, adjust_kde=16, lower.limit=0)
estimate_pdf("MinFeret", "Post", "kde", adjust_bde=1/2, adjust_kde=16, lower.limit=0)
estimate_pdf("Skew", "Pre", "kde", adjust_bde=1/4, adjust_kde=4, lower.limit=-2.7, upper.limit=5)
estimate_pdf("Skew", "Post", "kde", adjust_bde=1/4, adjust_kde=4, lower.limit=-2.7, upper.limit=5)
estimate_pdf("Kurt", "Pre", "kde", adjust_bde=1/4, adjust_kde=4,  lower.limit=-2, upper.limit=10)
estimate_pdf("Kurt", "Post", "kde", adjust_bde=1/4, adjust_kde=4,  lower.limit=-2, upper.limit=10)
estimate_pdf("AR", "Pre", "kde", adjust_bde=1/2, adjust_kde=6, lower.limit=1)
estimate_pdf("AR", "Post", "kde", adjust_bde=1/2, adjust_kde=6, lower.limit=1)
estimate_pdf("Round", "Pre", "bde", adjust_bde=1/2, adjust_kde=12, lower.limit=0, upper.limit=1)
estimate_pdf("Round", "Post", "bde", adjust_bde=1/2, adjust_kde=12, lower.limit=0, upper.limit=1)
estimate_pdf("Circularity", "Pre", "bde", adjust_bde=1/2, adjust_kde=4, lower.limit=0, upper.limit=1)
estimate_pdf("Circularity", "Post", "bde", adjust_bde=1/2, adjust_kde=4, lower.limit=0, upper.limit=1)
estimate_pdf("Solidity", "Pre", "bde", adjust_bde=1/4, adjust_kde=8, upper.limit=1)
estimate_pdf("Solidity", "Post", "bde", adjust_bde=1/4, adjust_kde=8, upper.limit=1)
estimate_pdf("FeretAngle", "Pre", "bde", adjust_bde=1, adjust_kde=8, lower.limit=0, upper.limit=180)
estimate_pdf("FeretAngle", "Post", "bde", adjust_bde=1, adjust_kde=8, lower.limit=0, upper.limit=180)
estimate_pdf("Angle", "Pre", "bde", adjust_bde=2, adjust_kde=8, lower.limit=0, upper.limit=180)
estimate_pdf("Angle", "Post", "bde", adjust_bde=2, adjust_kde=8, lower.limit=0, upper.limit=180)

#features[!features %in% dplyr::combine(features_signal, features_signal_bg, features_zero)]

gen_log("pdf done", logfile, datetime=TRUE)


geometry <- list(h_Area_Pre, h_Area_Post, 
				 h_IntDen_Pre, h_IntDen_Post, h_IntDen_bg_Pre, h_IntDen_bg_Post, 
				 h_RawIntDen_Pre, h_RawIntDen_Post, h_RawIntDen_bg_Pre, h_RawIntDen_bg_Post)  

h_geometry <- grid_arrange_shared_legend(geometry, ncol = 4, nrow = 3)
ggsave("Feature distribution (geometry).pdf", plot=h_geometry , width = 12, height = 6)

signal <- list(h_Mean_Pre, h_Mean_bg_Pre, h_Mean_z_Pre,
			   h_Mean_Post,  h_Mean_bg_Post, h_Mean_z_Post,
			   h_Median_Pre, h_Median_bg_Pre, h_Median_z_Pre,
			   h_Median_Post, h_Median_bg_Post, h_Median_z_Post,
			   h_Mode_Pre, h_Mode_bg_Pre, h_Mode_z_Pre,
			   h_Mode_Post, h_Mode_bg_Post, h_Mode_z_Post,
			   h_StdDev_Pre, h_StdDev_bg_Pre, h_StdDev_z_Pre,
			   h_StdDev_Post, h_StdDev_bg_Post, h_StdDev_z_Post,
			   h_Min_Pre, h_Min_bg_Pre, h_Min_z_Pre,
			   h_Min_Post, h_Min_bg_Post, h_Min_z_Post,
			   h_Max_Pre, h_Max_bg_Pre, h_Max_z_Pre,
			   h_Max_Post,  h_Max_bg_Post, h_Max_z_Post,
			   h_Skew_Pre, h_Skew_Post, h_Kurt_Pre, h_Kurt_Post)

h_signal1 <- grid_arrange_shared_legend(signal[1:18], ncol = 3, nrow = 6)
ggsave("Feature distribution (signal) 1.pdf", plot=h_signal1 , width = 12, height = 12)

h_signal2 <- grid_arrange_shared_legend(signal[19:36], ncol = 3, nrow = 6)
ggsave("Feature distribution (signal) 2.pdf", plot=h_signal2 , width = 12, height = 12)

h_signal3 <- grid_arrange_shared_legend(signal[37:40], ncol = 4, nrow = 1)
ggsave("Feature distribution (signal) 3.pdf", plot=h_signal3 , width = 12, height = 4)


morphology <- list(h_Perimeter_Pre, h_Perimeter_Post, h_shape_Pre, h_shape_Post, 
				   h_Width_Pre, h_Width_Post, h_Height_Pre, h_Height_Post, 
				   h_Major_Pre, h_Major_Post, h_Minor_Pre, h_Minor_Post, 
				   h_Angle_Pre, h_Angle_Post, h_AR_Pre, h_AR_Post,
				   h_Round_Pre, h_Round_Post, h_Circularity_Pre, h_Circularity_Post, 
				   h_Solidity_Pre, h_Solidity_Post, h_FeretAngle_Pre, h_FeretAngle_Post,
				   h_Feret_Pre, h_Feret_Post, h_MinFeret_Pre, h_MinFeret_Post)
h_morphology1 <- grid_arrange_shared_legend(morphology[1:14], ncol = 4, nrow = 4)
ggsave("Feature distribution (morphology) 1.pdf", plot=h_morphology1 , width = 12, height = 12)

h_morphology2 <- grid_arrange_shared_legend(morphology[15:28], ncol = 4, nrow = 4)
ggsave("Feature distribution (morphology) 2.pdf", plot=h_morphology2 , width = 12, height = 12)

gen_log("Ideal Marker pdf estmation done", logfile, datetime=TRUE)


ggplot(data, aes(x=X, y=Y, color=group))+geom_point(aes(size=Area, alpha=1/255*Mean))