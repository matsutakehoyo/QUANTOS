library(MASS)                     #for pi breaks 2dkdes
library(grid)					  #grid_arrange_shared_legend
library(gridExtra)                #grid_arrange_shared_legend
library(propagate)                #for distDistr
library(ks)
library(plyr)                     #for mapvalues to replace NaN
library(tidyverse)                #ggplot, tibble, tidyr, readr, purrr, dplyr
library(bde)

dist_cutoff_min = 0
dist_cutoff_max = 1.2
distance_density = 4                                 #distance for marker density calculation

#type of image 3D or 2D
image_3D = FALSE

#log R output ot file: moified original file by Sunagawa gen-log.R
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
vplayout <- function(x,y) viewport(layout.pos = x, layout.pos.col =y)

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

#type of image 3D or 2D
image_3D = FALSE

#log results in file
logfile <- gen_start_log("FindSynapse")
gen_log("Applying Niave Bayes Classifier to paired synapse markers", logfile, datetime=TRUE)


#load data
data <- readRDS("data_bg.rds")
data_pair <- readRDS("data_pair_bg.rds")

#remove na rows
na <- na.omit(data_pair) %>% na.action() %>% as.integer
na_id <- data_pair[na, ] %>% dplyr::select(pair_id) %>% pull()
data_pair <- data_pair %>% filter(!pair_id %in% na_id)
gen_log(paste("Removed", length(na_id), "NA/NaN pairs"), logfile)


if (length(data_pair)){
  background_pre <- data_pair %>% filter(group=="Pre") %>%dplyr::select(bg) %>% unique() %>% pull()
  background_post <- data_pair %>% filter(group=="Post") %>%dplyr::select(bg) %>% unique() %>% pull()
}

calc_global_density <- function(group_str){
  n <- data %>% filter(group==group_str) %>% nrow()
  gen_log(paste("Number of", group_str, "markers:", n), logfile)	
  x_min <- data %>% filter(group==group_str) %>%dplyr::select(X) %>% min()
  x_max <- data %>% filter(group==group_str) %>%dplyr::select(X) %>% max()
  y_min <- data %>% filter(group==group_str) %>%dplyr::select(Y) %>% min()
  y_max <- data %>% filter(group==group_str) %>%dplyr::select(Y) %>% max()
  
  n / ((x_max - x_min) * (y_max - y_min))
}

dPre <- calc_global_density("Pre")
gen_log(paste("Pre density:", dPre), logfile)
dPost <- calc_global_density("Post")
gen_log(paste("Post density:", dPost), logfile)

#read feature pdf
read_pdf <- function(filename){
  files <- list.files(path = "./pdf/", pattern=paste0(filename, "*.rds"), recursive=TRUE)
  #print(files) 
  for (i in seq_along(files)){
	print(paste("Loading file:", files[i] ))
	pdf <- readRDS(paste0("pdf/", files[i]))
	assign(gsub(pattern=".rds|pdf/", "", x=files[i]), pdf, envir=globalenv())
	print(gsub(pattern=".rds|pdf/", "", x=files[i]))
  }
  files <- gsub(pattern=".rds|pdf/", "", x=files)
}
read_pdf("pdf.synapse")
pdf_noise <- dplyr::combine(read_pdf("bde.noise"), read_pdf("kde.noise"))
pdf_synapse <- dplyr::combine(read_pdf("bde.synapse"), read_pdf("kde.synapse"))

if (length(pdf_noise)!=length(pdf_synapse)){
  gen_log("There are missing pdf's !!", logfile)	
  stop("The value is TRUE, so the script must end here")
}
#capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# stddev_bg_pre_pdf + geom_histogram(data=(data_pair %>% filter(group=="Pre")), aes(x=StdDev_bg, y=..density..), alpha=1/3)
# ggplot()+ geom_histogram(data=(data_pair %>% filter(group=="Pre")), aes(x=Mean_bg, y=..density..))

#plot synapse & noise pdf
dist_pdf <- ggplot(data = data.frame(x = c(dist_cutoff_min, dist_cutoff_max)), aes(x), fill=eval) +
  stat_function(fun = dnorm, n = 200, 
				args = list(mean = dist.pdf.synapse$estimate[1],
							sd=dist.pdf.synapse$estimate[2]),
				geom="area", alpha = 1/3, aes(fill="synapse")) +
  stat_function(fun = function(x){3/dist_cutoff_max^3*x^2}, geom="area", alpha = 1/3, aes(fill="noise")) +
  xlab("distance (um)") + ylab("probability") + ggtitle("Distance") +
  scale_fill_manual("Marker", values = c("synapse" ="#56B4E9","noise" = "#E69F00")) +
  basic_style
if (length(data_pair)){
  dist_pdf <- dist_pdf + 
	geom_histogram(data=(data_pair), aes(x = pair_dist, y = ..density..), alpha=1/3) 
}  

angle_pdf <- ggplot(data = data.frame(x = c(-pi, pi)), aes(x)) +
  stat_function(fun = propagate:::dlaplace, n = 200, args = list(mean = 0, sigma=pi/2), geom="area", alpha = 1/3, aes(fill="synapse")) +
  stat_function(fun = dunif, n = 200, args = list(min=-pi, max=pi), geom="area", alpha = 1/3, fill="#E69F00") + 
  xlab("angle") + ylab("probability") + ggtitle("Angle")+
  scale_fill_manual("Marker", values = c("synapse" ="#56B4E9","noise" = "#E69F00")) +
  basic_style 
if (length(data_pair)){
  angle_pdf <- angle_pdf +
	geom_histogram(data=(data_pair), aes(x = pair_angle, y = ..density..), alpha=1/3) 
}
angle_pdf <- fracAx(angle_pdf, "pi")

robust_hpi <- function(data){
  tryCatch({
	min <- data %>% na.omit() %>% pull() %>% min()
	hpi <- hpi(data %>% na.omit() %>% pull())
	hpi
  }, warning = function(w){
	print(paste("Warning in hpi estimation"))
  }, error = function(err){
	print(paste("Error in hpi estimation using min value insted:", min))
	min
  }
  )
}

#function to plot the pdf for noise and synapse pairs
plot_pdf <- function(pdf){
  title <- gsub(pattern=".bde|kde", replacement="", pdf) %>% 
	strsplit("[.]") %>% unlist() %>%
	firstup(.) %>% paste(collapse=" ")
  feature <- pdf %>% strsplit("[.]") %>% unlist() %>% .[1]
  group_name <- pdf %>% strsplit("[.]") %>% unlist() %>% .[2] 
  pdf_type <- pdf %>% strsplit("[.]") %>% unlist() %>% .[3]
  #print(paste(group_name, feature))
  
  if (pdf_type=="bde"){
	min <- eval(parse(text=paste0(pdf, "@lower.limit"))) %>% floor
	max <- eval(parse(text=paste0(pdf, "@upper.limit"))) %>% ceiling
  } else if (pdf_type=="kde"){
	min <- eval(parse(text=paste0(pdf, "$eval.points"))) %>% min %>% floor
	max <- eval(parse(text=paste0(pdf, "$eval.points"))) %>% max %>% ceiling
  }
  #print(paste("Min:", min, "Max:", max))
  p <- ggplot(data = data.frame(x = c(min, max)), aes(x)) +
	scale_fill_manual(values = c("synapse" ="#56B4E9","noise" = "#E69F00")) +
	xlim(min, max) +
	xlab(feature) + ylab("probability") + ggtitle(title) + 
	basic_style + theme(legend.position="none")

  if (length(data_pair)){
	bw <- robust_hpi(data_pair %>% filter(group==group_name) %>% dplyr::select(feature))
	p <- p +
	geom_histogram(data = (data_pair %>% filter(group == group_name)), 
				   aes(x = eval(parse(text = feature)), y = ..density..), 
				   binwidth = bw, alpha = 1/3)
  }
  
  if (pdf_type=="bde"){
	p <- p +
	  stat_function(fun = bde::density, n = 200, args = list(x = eval(parse(text=paste0(feature, ".", group_name, ".", pdf_type, ".synapse")))), 
					geom="area", alpha = 1/3, aes(fill="synapse")) +    
	  stat_function(fun = bde::density, n = 200, args = list(x = eval(parse(text=paste0(feature, ".", group_name, ".", pdf_type, ".noise")))), 
					geom="area", alpha = 1/3, aes(fill="noise")) 
  } else if (pdf_type=="kde"){
	p <- p +
	  stat_function(fun = dkde, n = 200, args = list(fhat = eval(parse(text=paste0(feature, ".", group_name, ".", pdf_type, ".synapse")))), 
					geom="area", alpha = 1/3, aes(fill="synapse")) +    
	  stat_function(fun = dkde, n = 200, args = list(fhat = eval(parse(text=paste0(feature, ".", group_name, ".", pdf_type, ".noise")))), 
					geom="area", alpha = 1/3, aes(fill="noise")) 
  }
}

for (pdf in pdf_synapse){
  feature <- pdf %>% strsplit("[.]") %>% unlist() %>% .[1]
  group_name <- pdf %>% strsplit("[.]") %>% unlist() %>% .[2] 
  plot <- plot_pdf(pdf)
  assign(paste0(feature, "_", group_name, "_pdf"), plot)
}



pdf_geometry <- list(dist_pdf, angle_pdf, 
					 Area_Pre_pdf, Area_Post_pdf, 
					 density_Pre_pdf, density_Post_pdf,
					 IntDen_Pre_pdf, IntDen_bg_Pre_pdf, IntDen_z_Pre_pdf,
					 IntDen_Post_pdf, IntDen_bg_Post_pdf, IntDen_z_Post_pdf,
					 RawIntDen_Pre_pdf, RawIntDen_bg_Pre_pdf, RawIntDen_z_Pre_pdf,
					 RawIntDen_Post_pdf, RawIntDen_bg_Post_pdf, RawIntDen_z_Post_pdf)

pdf_signal <- list(Mean_Pre_pdf, Mean_bg_Pre_pdf, Mean_z_Pre_pdf,
				   Mean_Post_pdf, Mean_bg_Post_pdf, Mean_z_Post_pdf, 
				   Median_Pre_pdf, Median_bg_Pre_pdf, Median_z_Pre_pdf, 
				   Median_Post_pdf, Median_bg_Post_pdf, Median_z_Post_pdf, 
				   Mode_Pre_pdf, Mode_bg_Pre_pdf, Mode_z_Pre_pdf,
				   Mode_Post_pdf, Mode_bg_Post_pdf, Mode_z_Post_pdf, 
				   Min_Pre_pdf, Min_bg_Pre_pdf,  Min_z_Pre_pdf,
				   Min_Post_pdf, Min_bg_Post_pdf, Min_z_Post_pdf, 
				   Max_Pre_pdf, Max_bg_Pre_pdf, Max_z_Pre_pdf,
				   Max_Post_pdf, Max_bg_Post_pdf, Max_z_Post_pdf, 
				   StdDev_Pre_pdf, StdDev_bg_Pre_pdf, StdDev_z_Pre_pdf, 
				   StdDev_Post_pdf, StdDev_bg_Post_pdf, StdDev_z_Post_pdf, 
				   Skew_Pre_pdf, Skew_Post_pdf, Kurt_Pre_pdf, Kurt_Post_pdf)

pdf_morphology <- list(Perimeter_Pre_pdf, Perimeter_Post_pdf, shape_Pre_pdf, shape_Post_pdf, 
					   Width_Pre_pdf, Width_Post_pdf, Height_Pre_pdf, Height_Post_pdf, 
					   Major_Pre_pdf, Major_Post_pdf, Minor_Pre_pdf, Minor_Post_pdf, 
					   Angle_Pre_pdf, Angle_Post_pdf, AR_Pre_pdf, AR_Post_pdf, 
					   Round_Pre_pdf, Round_Post_pdf, Circularity_Pre_pdf, Circularity_Post_pdf,
					   Solidity_Pre_pdf, Solidity_Post_pdf, FeretAngle_Pre_pdf, FeretAngle_Post_pdf, 
					   Feret_Pre_pdf, Feret_Post_pdf, MinFeret_Pre_pdf, MinFeret_Post_pdf)

h_pdf_geometry1 <- grid_arrange_shared_legend(pdf_geometry[1:6], ncol = 2, nrow = 3)
ggsave("ROI sats probability distribution (geometry1).pdf", plot=h_pdf_geometry1 , width = 15, height = 9)

h_pdf_geometry2 <- grid_arrange_shared_legend(pdf_geometry[7:18], ncol = 3, nrow = 4)
ggsave("ROI sats probability distribution (geometry2).pdf", plot=h_pdf_geometry2 , width = 15, height = 10)


h_pdf_signal1 <- grid_arrange_shared_legend(pdf_signal[1:18], ncol = 3, nrow = 6)
ggsave("ROI sats probability distribution (signal1).pdf", plot=h_pdf_signal1 , width = 15, height = 10)

h_pdf_signal2 <- grid_arrange_shared_legend(pdf_signal[19:36], ncol = 3, nrow = 6)
ggsave("ROI sats probability distribution (signal2).pdf", plot=h_pdf_signal2 , width = 15, height = 10)

h_pdf_signal3 <- grid_arrange_shared_legend(pdf_signal[37:40], ncol = 2, nrow = 2)
ggsave("ROI sats probability distribution (signal3).pdf", plot=h_pdf_signal3 , width = 15, height = 6)

h_pdf_morphology <- grid_arrange_shared_legend(pdf_morphology, ncol = 4, nrow = 7)
ggsave("ROI sats probability distribution (shape).pdf", plot=h_pdf_morphology , width = 15, height = 10)


#function to score each feature according to its pdf, subset the data to the group 
calc_score <- function(pdf){
  score_name <- gsub(pattern=".(Pre|Post).(bde|kde).(synapse|noise)", replacement="", pdf)
  group_name <- pdf %>% strsplit("[.]") %>% unlist() %>% .[2] %>% firstup(.) #pre or post
  pdf_type <- pdf %>% strsplit("[.]") %>% unlist() %>% .[3] #bde or kde
  marker_type <- pdf %>% strsplit("[.]") %>% unlist() %>% .[4] #noise or synapse
  
  if (pdf_type=="bde"){
	data_scored <- eval(parse(text=paste0("dplyr::mutate(filter(data_pair_scored, group==\"",group_name, "\"), 
										  score_", score_name, "_", marker_type, "=",  
										  "bde::density(values = ", score_name, ", x = ", pdf, "))")))
  } else if (pdf_type=="kde"){
	data_scored <- eval(parse(text=paste0("dplyr::mutate(filter(data_pair_scored, group==\"",group_name, "\"), 
										  score_", score_name, "_", marker_type, "=",  
										  "dkde(x = ", score_name, ", fhat = ", pdf, "))")))
  }
  #na <- data_scored %>% dplyr::select(data_scored %>% dim() %>% .[2]) %>% is.na() %>% sum()
  #zero <- data_scored %>% dplyr::select(starts_with("score")) %>% filter_all(any_vars(. == 0)) %>% dim() %>% .[1]
  #zero <- data_scored %>% dplyr::select(data_scored %>% dim() %>% .[2]) %>% filter(.==0) %>% dim() %>% .[1]
  
  # if (na){
  # 	gen_log(paste0("Scoring: ",  score_name, "   ", group_name, "   ", pdf_type), logfile)
  #   gen_log(paste0(na, "  NaNs generated"), logfile)  
  # }
  # 
  # if (zero){
  # 	#data_scored %>% dplyr::select(data_scored %>% dim() %>% .[2]) %>% filter(.==0) %>% print(width=Inf)
  # 	gen_log(paste0("Scoring: ",  score_name, "   ", group_name, "   ", pdf_type), logfile)
  #   gen_log(paste0(zero, "  zeros generated"), logfile)
  # }
  #print(data_scored)
  #print(data_scored %>% dplyr::select(data_pair_scored %>% dim() %>% .[2]))
  d <- dplyr::bind_rows(filter(data_pair_scored, group!=group_name), data_scored)
}

pdf_list <- dplyr::combine(pdf_synapse, pdf_noise)

if (!(nrow(data_pair))){
  gen_log("Empty data", logfile)
  data_pair <- tibble(X = as.numeric(), Y = as.numeric(), 
					  XM = as.numeric(), YM = as.numeric(), 
					  pair_dist = as.numeric(), 
					  pair_angle = as.numeric(), 
					  pair_id = as.integer(), 
					  marker_id = as.integer(), 
					  group = as.character())
  features <- pdf_list %>% strsplit("[.]") %>% lapply(., `[[`, 1) %>% unlist() %>% unique()
  for (feature in features){
	data_pair <- eval(parse(text=paste0("add_column(data_pair, ", feature, "= as.numeric())"))) 
  }
}

#calculate synapse probability
data_pair_scored <- data_pair %>% 
  #distance
  mutate(score_dist_synapse = dnorm(pair_dist, 
									mean = dist.pdf.synapse$estimate[1], 
									sd = dist.pdf.synapse$estimate[2])) %>%
  mutate(score_dist_noise = 3 / dist_cutoff_max^3 * pair_dist^2) %>%
  #angle
  #mutate(score_angle_synapse = dnorm(pair_angle, mean = 0, sd = pi/4)) %>%
  mutate(score_angle_synapse = propagate:::dlaplace(pair_angle, mean = 0, sigma = pi/2)) %>%
  mutate(score_angle_noise = dunif(pair_angle, min= -pi, max= pi))
#distance is calculated manually/normal distribution  


for (pdf in pdf_list){
  data_pair_scored <- calc_score(pdf)
  # gen_log(paste0("Scoring ", pdf), logfile)
}

#fix zero frequency problem
for (score in colnames(data_pair_scored)[colnames(data_pair_scored) %>% grep("score_", .)]){
	zero <- data_pair_scored %>%dplyr::select(score) %>% filter(.==0) %>% dim() %>% .[1]
	if (zero){
		feature_name <- gsub(pattern = "score_", "", score) %>% 
		gsub(pattern = "_synapse|_noise", "", .)
		signal_name <- score %>% strsplit("[_]") %>% unlist() %>% .[length(.)]
		gen_log(paste(zero, "zeros in score:", feature_name, signal_name), logfile)
		zero_rows <- which(data_pair_scored %>% dplyr::select(score) ==0) %>% print
		
		if (signal_name == "synapse"){
		  opposite_score <- paste0("score_", feature_name, "_noise")
		}else if (signal_name == "noise"){
		  opposite_score <- paste0("score_", feature_name, "_synapse")
		}
		data_pair_scored[zero_rows,score] <-0.5
		data_pair_scored[zero_rows,opposite_score] <-0.5
	}
	  
	na <- data_pair_scored %>% dplyr::select(score) %>% is.na() %>% sum()
	if (na){
		feature_name <- gsub(pattern = "score_", "", score) %>% 
		  gsub(pattern = "_synapse|_noise", "", .)
		signal_name <- score %>% strsplit("[_]") %>% unlist() %>% .[length(.)]
		gen_log(paste(na, "NA in score:", feature_name, signal_name), logfile)
		na_rows <- which(data_pair_scored %>% dplyr::select(score) %>% is.na()) %>% print
		
		if (signal_name == "synapse"){
		  opposite_score <- paste0("score_", feature_name, "_noise")
		}else if (signal_name == "noise"){
		  opposite_score <- paste0("score_", feature_name, "_synapse")
		}
		data_pair_scored[na_rows,score] <-0.5
		data_pair_scored[na_rows,opposite_score] <-0.5
	}
	  
	nan <- data_pair_scored %>% dplyr::select(score) %>% pull %>% is.nan() %>% sum()
	  
	if (nan){
		feature_name <- gsub(pattern = "score_", "", score) %>% 
		  gsub(pattern = "_synapse|_noise", "", .)
		signal_name <- score %>% strsplit("[_]") %>% unlist() %>% .[length(.)]
		gen_log(paste(na, "NaN in score:", feature_name, signal_name), logfile)
		nan_rows <- which(data_pair_scored %>% dplyr::select(score) %>% pull() %>% is.nan()) %>% print
		
		if (signal_name == "synapse"){
		  opposite_score <- paste0("score_", feature_name, "_noise")
		}else if (signal_name == "noise"){
		  opposite_score <- paste0("score_", feature_name, "_synapse")
		}
		data_pair_scored[nan_rows,score] <-0.5
		data_pair_scored[nan_rows,opposite_score] <-0.5
	}
}



#geometry score
data_pair_scored_pre <- data_pair_scored %>% 
	filter(group=="Pre")
data_pair_scored_pre <- data_pair_scored_pre %>%	
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_synapse")) %>% 
	dplyr::select(score_dist_synapse, score_angle_synapse
				,score_density_synapse
				,score_Area_synapse
				,score_IntDen_bg_synapse
				,score_RawIntDen_bg_synapse
				) %>%
	transmute(score_geometry_synapse = apply(., 1, prod)) %>%
	bind_cols(data_pair_scored_pre, .)

data_pair_scored_pre <- data_pair_scored_pre %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_noise")) %>%
	dplyr::select(score_dist_noise, score_angle_noise
				,score_density_noise
				,score_Area_noise
				,score_IntDen_bg_noise
				,score_RawIntDen_bg_noise
				) %>%
	transmute(score_geometry_noise = apply(., 1, prod)) %>%
	bind_cols(data_pair_scored_pre, .)

data_pair_scored_post <- data_pair_scored %>% 
	filter(group=="Post") 
data_pair_scored_post <- data_pair_scored_post %>%
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_synapse")) %>% 
	dplyr::select(score_dist_synapse, score_angle_synapse
				  ,score_density_synapse
				  ,score_Area_synapse
				  ,score_IntDen_bg_synapse
				  ,score_RawIntDen_bg_synapse
	) %>%
	transmute(score_geometry_synapse = apply(., 1, prod)) %>%
	bind_cols(data_pair_scored_post, .)

data_pair_scored_post <- data_pair_scored_post %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_noise")) %>%
	dplyr::select(score_dist_noise, score_angle_noise
				,score_density_noise
				,score_Area_noise
				,score_IntDen_bg_noise
				,score_RawIntDen_bg_noise
				) %>%
	transmute(score_geometry_noise = apply(., 1, prod)) %>%
	bind_cols(data_pair_scored_post, .)

#signal score
data_pair_scored_pre <- data_pair_scored_pre %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_synapse")) %>% 
	dplyr::select(
				  score_Mean_synapse
				  ,score_Median_synapse
				  ,score_Mode_synapse
				  ,score_Mean_bg_synapse
				  ,score_Median_bg_synapse
				  ,score_Mode_bg_synapse
				  ,score_Min_synapse
				  ,score_StdDev_bg_synapse
				  ,score_Max_synapse
				  ,score_Max_bg_synapse
				  ,score_Min_bg_synapse
				  ,score_Skew_synapse
				  ,score_Kurt_synapse 
  				) %>%
	transmute(score_signal_synapse = apply(., 1, prod)) %>%
	bind_cols(data_pair_scored_pre, .)

data_pair_scored_pre <- data_pair_scored_pre %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_noise")) %>%
	dplyr::select(
				  score_Mean_noise
				  ,score_Median_noise
				  ,score_Mode_noise
				  ,score_Mean_bg_noise
				  ,score_Median_bg_noise
				  ,score_Mode_bg_noise
				  ,score_Min_noise
				  ,score_StdDev_bg_noise
				  ,score_Max_noise
				  ,score_Max_bg_noise
				  ,score_Min_bg_noise
				  ,score_Skew_noise
				  ,score_Kurt_noise
  				 ) %>%
	transmute(score_signal_noise = apply(., 1, prod)) %>% 
	bind_cols(data_pair_scored_pre, .)

data_pair_scored_post <- data_pair_scored_post %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_synapse")) %>% 
	dplyr::select(
		score_Mean_synapse
		,score_Median_synapse
		,score_Mode_synapse
		,score_Mean_bg_synapse
		,score_Median_bg_synapse
		,score_Mode_bg_synapse
		,score_Max_synapse
		,score_Min_synapse
		,score_StdDev_bg_synapse
		,score_Max_bg_synapse
		,score_Min_bg_synapse
		,score_Skew_synapse
		,score_Kurt_synapse 
	) %>%
	transmute(score_signal_synapse = apply(., 1, prod)) %>%
	bind_cols(data_pair_scored_post, .)

data_pair_scored_post <- data_pair_scored_post %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_noise")) %>%
	dplyr::select(
		score_Mean_noise
		,score_Median_noise
		,score_Mode_noise
		,score_Mean_bg_noise
		,score_Median_bg_noise 
		,score_Mode_bg_noise
		,score_Max_noise
		,score_Min_noise
		,score_StdDev_bg_noise
		,score_Max_noise
		,score_Max_bg_noise
		,score_Min_bg_noise
		,score_Skew_noise
		,score_Kurt_noise
	) %>%
	transmute(score_signal_noise = apply(., 1, prod)) %>% 
	bind_cols(data_pair_scored_post, .)

#morphology score
data_pair_scored_pre <- data_pair_scored_pre %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_synapse")) %>% 
	dplyr::select(
				 	score_Perimeter_synapse
				 	,score_shape_synapse
					,score_Width_synapse
					,score_Height_synapse			#bounding box
					,score_Major_synapse
					,score_Minor_synapse
					,score_Angle_synapse
					,score_AR_synapse
					,score_Round_synapse
					,score_Circularity_synapse
					,score_Solidity_synapse
					,score_Feret_synapse
					,score_MinFeret_synapse
					,score_FeretAngle_synapse
				 ) %>%
  transmute(score_morph_synapse = apply(., 1, prod)) %>% 
  bind_cols(data_pair_scored_pre, .)

data_pair_scored_pre <- data_pair_scored_pre %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_noise")) %>%
	dplyr::select(
					score_Perimeter_noise
					,score_shape_noise
					,score_Width_noise
					,score_Height_noise				#bounding box
					,score_Major_noise
					,score_Minor_noise
					,score_Angle_noise
					,score_AR_noise
					,score_Round_noise
					,score_Circularity_noise
					,score_Solidity_noise
					,score_Feret_noise
					,score_MinFeret_noise
					,score_FeretAngle_noise
  				 ) %>%
	transmute(score_morph_noise = apply(., 1, prod)) %>% 
	bind_cols(data_pair_scored_pre, .)

data_pair_scored_post <- data_pair_scored_post %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_synapse")) %>% 
	dplyr::select(
		score_Perimeter_synapse
		,score_shape_synapse
		,score_Width_synapse
		,score_Height_synapse			#bounding box
		,score_Major_synapse
		,score_Minor_synapse
		,score_Angle_synapse
		,score_AR_synapse
		,score_Round_synapse
		,score_Circularity_synapse
		,score_Solidity_synapse
		,score_Feret_synapse
		,score_MinFeret_synapse
		,score_FeretAngle_synapse
	) %>%
	transmute(score_morph_synapse = apply(., 1, prod)) %>% 
	bind_cols(data_pair_scored_post, .)

data_pair_scored_post <- data_pair_scored_post %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_noise")) %>%
	dplyr::select(
		score_Perimeter_noise
		,score_shape_noise
		,score_Width_noise
		,score_Height_noise				#bounding box
		,score_Major_noise
		,score_Minor_noise
		,score_Angle_noise
		,score_AR_noise
		,score_Round_noise
		,score_Circularity_noise
		,score_Solidity_noise
		,score_Feret_noise
		,score_MinFeret_noise
		,score_FeretAngle_noise
	) %>%
	transmute(score_morph_noise = apply(., 1, prod)) %>% 
	bind_cols(data_pair_scored_post, .)

data_pair_scored <- bind_rows(data_pair_scored_pre, data_pair_scored_post)

#total score
data_pair_scored <- data_pair_scored %>%  
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_synapse")) %>% 
	dplyr::select(score_geometry_synapse, score_signal_synapse, score_morph_synapse) %>%
	transmute(score_all_synapse = apply(., 1, prod)) %>%
	bind_cols(data_pair_scored, .)

data_pair_scored <- data_pair_scored %>% 
	dplyr::select(starts_with("score_")) %>% 
	dplyr::select(ends_with("_noise")) %>%
	dplyr::select(score_geometry_noise, score_signal_noise, score_morph_noise) %>%
	transmute(score_all_noise = apply(., 1, prod)) %>%
	bind_cols(data_pair_scored, .)

fix_eval <- function(df, str){
	tryCatch(
		{
			mutate(df, eval=str)
		}, warning = function(w){
			print("Warning in mutate")
		}, error = function(err){
			print("Error in mutate:")
			df
		}
  	)
}

#marker features score
marker_features_syn <- data_pair_scored %>% 
	dplyr::select(pair_id, marker_id, group,  ends_with("_synapse")) %>%
	gather(feature, score, matches("score"), factor_key = TRUE) %>% 
	filter(feature!="score_all_synapse") %>%
	separate(feature, c("a", "feature", "eval"), sep = c(6,-8)) %>%  #cannnot separate at "_" because _bg
	dplyr::select(-a) %>% 
	mutate(eval="synapse")

marker_features_noi <- data_pair_scored %>% 
	dplyr::select(pair_id, marker_id, group,  ends_with("_noise")) %>%
	gather(feature, score, matches("score"), factor_key = TRUE) %>%
	filter(feature!="score_all_noise") %>%
	separate(feature, c("a", "feature", "eval"), sep = c(6, -6)) %>% 
	dplyr::select(-a) %>% 
	mutate(eval="noise")

marker_features <- bind_rows(marker_features_syn, marker_features_noi) %>%
	filter(score>1e-300) #this avoid errors when plotting

s_marker_features <- ggplot(marker_features) + basic_style
if (nrow(marker_features)){
	s_marker_features <- s_marker_features + 
		geom_point(data=marker_features %>% 
						spread(eval, score) %>% 
						filter(feature!="geometry" & feature!="signal" & feature!="morph"), 
				aes(x=synapse, y=noise, color=group), size=0.1) +
	geom_abline(intercept=0, slope=1, color="darkgray") +
	scale_x_log10()+
	scale_y_log10()+
	facet_wrap(~feature, scales="free", ncol=6)      
}
ggsave("Synapse-Random score scatter.pdf", width=15, height=10)

s_marker_features_summary <- ggplot() + basic_style
if (nrow(marker_features)){
	s_marker_features_summary <- s_marker_features_summary + 
	geom_point(data = marker_features %>% 
				spread(eval, score) %>% 
				filter(feature=="geometry" | feature=="signal" | feature=="morph"), 
			aes(x=synapse, y=noise, color=group), size=0.1) +
	geom_abline(intercept=0, slope=1, color="darkgray") +
	scale_x_log10()+
	scale_y_log10()+
	facet_wrap(~feature, scales="free", ncol=6) 
}
ggsave("Synapse-Random score scatter summary.pdf", width=15, height=10)


#calculate probability of random/synapse
# poly(dist.max, 2, raw = TRUE)1 poly(dist.max, 2, raw = TRUE)2 
#                        877.599                       7136.395 
rand.pair <- function(dist, dens1, dens2){
	(877.599 * dist + 7136.395 * dist^2) * dens1 * dens2
}

synapse_prior <- function(dist, dens1, dens2){
	ifelse(rand.pair(dist, dens1, dens2) < 2, 0.5, 1/rand.pair(dist, dens1, dens2))
}

#calculate Pre Post combined score
if (nrow(data_pair_scored)){
	data_pair_scored <- data_pair_scored %>% 
		dplyr::select(pair_id, group, score_all_synapse) %>%
		spread(group, score_all_synapse) %>%
		mutate(score_total_synapse = Pre * Post) %>%
		gather(group, score_marker, Pre:Post) %>% 
		arrange(pair_id) %>% dplyr::select(score_total_synapse) %>%
		bind_cols(data_pair_scored %>% arrange(pair_id), .)
  
	data_pair_scored <- data_pair_scored %>% 
		dplyr::select(pair_id, group, score_all_noise) %>%
		spread(group, score_all_noise) %>%
		mutate(score_total_noise = Pre * Post) %>%
		gather(group, score_marker, Pre:Post) %>% 
		arrange(pair_id) %>% dplyr::select(score_total_noise) %>%
		bind_cols(data_pair_scored %>% arrange(pair_id), .)
  
	#calculate Pre Post combined prior
	data_pair_scored <- data_pair_scored %>% 
		dplyr::select(pair_id, group, pair_dist, density) %>%
		spread(group, density) %>%
		#mutate(prior_synapse = synapse_prior(pair_dist, Pre, Post)) %>%
	  mutate(prior_synapse = synapse_prior(5, Pre, Post)) %>%
		gather(group, synapse_prior, Pre:Post) %>% 
		arrange(pair_id) %>% dplyr::select(prior_synapse) %>%
		bind_cols(data_pair_scored %>% arrange(pair_id), .)
  
  #calculate naive_bayes synapse
	data_pair_scored <- data_pair_scored %>% 
		dplyr::select(pair_id, group, score_total_synapse, prior_synapse) %>%
		mutate(naive_bayes_synapse = score_total_synapse * prior_synapse) %>%
		arrange(pair_id) %>% dplyr::select(naive_bayes_synapse) %>%
		bind_cols(data_pair_scored %>% arrange(pair_id), .)
  
	data_pair_scored <- data_pair_scored %>% 
		dplyr::select(pair_id, group, score_total_noise, prior_synapse) %>%
		mutate(naive_bayes_noise = score_total_noise * (1-prior_synapse)) %>%
		arrange(pair_id) %>% dplyr::select(naive_bayes_noise) %>%
		bind_cols(data_pair_scored %>% arrange(pair_id), .)

	gen_log(paste("All pairs:", length(data_pair_scored$pair_id %>% unique)), logfile)

	#remove duplicates first remove duplicates from one side then the other to avoid orphan markers
	data_pair_scored %>% 
		mutate(naive_bayes = naive_bayes_synapse - naive_bayes_noise) %>% 
		dplyr::select(pair_id, marker_id, group, naive_bayes) %>%
		spread(group, marker_id) %>% 
		group_by(Pre) %>% 
		filter(naive_bayes == max(naive_bayes)) %>%
		ungroup() %>% 
		dplyr::select(pair_id) %>% pull %>% unique() -> unique_pair
  
  	data_pair_scored <- data_pair_scored %>% filter(pair_id %in% unique_pair)
  
	data_pair_scored %>% 
		mutate(naive_bayes = naive_bayes_synapse - naive_bayes_noise) %>% 
		dplyr::select(pair_id, marker_id, group, naive_bayes) %>%
		spread(group, marker_id) %>% 
		group_by(Post) %>% 
		filter(naive_bayes == max(naive_bayes)) %>%
		ungroup() %>% 
		dplyr::select(pair_id) %>% pull %>% unique() -> unique_pair
  
  	data_pair_scored <- data_pair_scored %>% filter(pair_id %in% unique_pair)
} else {
	data_pair_scored <- data_pair_scored %>% 
		add_column(score_total_synapse = as.double(),
				   score_total_noise = as.double(),
				   prior_synapse = as.double(), 
				   naive_bayes_synapse = as.double(),
				   naive_bayes_noise = as.double(), 
				   naive_bayes = as.double())
  	unique_pair <- as.integer()
	gen_log(paste("All pairs:", length(data_pair_scored$pair_id %>% unique)), logfile)
}

gen_log(paste("Unique pairs:", length(unique_pair)), logfile)

p_synapse_prior <- ggplot(data_pair_scored, aes(x=prior_synapse)) + geom_histogram()

if (nrow(data_pair_scored)){
	data_scatter_plot <- data_pair_scored %>% 
		dplyr::select(marker_id, pair_id, group, starts_with("score_")) %>% 
		gather(score, value, matches("score_"), factor_key = TRUE) %>% 
		filter(value>1e-300) %>%
		spread(score, value)
  
	s_signal_scatter <- ggplot(data_scatter_plot) +
		# geom_point(aes(x=score_Mean_synapse, y=score_Mean_noise, color="Mean"), size=1, alpha=1/5) +
		# geom_point(aes(x=score_Mode_synapse, y=score_Mode_noise, color="Mode"), size=1, alpha=1/5) +
		# geom_point(aes(x=score_Median_synapse, y=score_Median_noise, color="Median"), size=1, alpha=1/5) +  
		geom_point(aes(x=score_Mean_bg_synapse, y=score_Mean_bg_noise, color="Mean_bg"), size=1, alpha=1/5) +
		geom_point(aes(x=score_Mode_bg_synapse, y=score_Mode_bg_noise, color="Mode_bg"), size=1, alpha=1/5) +
		geom_point(aes(x=score_Median_bg_synapse, y=score_Median_bg_noise, color="Median_bg"), size=1, alpha=1/5) +
		geom_point(aes(x=score_StdDev_bg_synapse, y=score_StdDev_bg_noise, color="StdDev_bg"), size=1, alpha=1/5) +
		#geom_point(aes(x=score_Min_bg_synapse, y=score_Min_bg_noise, color="Min_bg")) +
		geom_point(aes(x=score_Max_synapse, y=score_Max_noise, color="Max"), size=1, alpha=1/5) +
		geom_point(aes(x=score_Max_bg_synapse, y=score_Max_bg_noise, color="Max_bg"), size=1, alpha=1/5) +
		scale_x_log10()+
		scale_y_log10()+  
		facet_grid(.~group)+
		xlab("synapse") + ylab("noise") +
		geom_abline(intercept=0, slope=1, color="grey60")+
		basic_style
  	ggsave("Score scatter plot (signal).pdf")
  
	s_morph_scatter <- ggplot(data_scatter_plot) +
		geom_point(aes(x=score_Perimeter_synapse, y=score_Perimeter_noise, color="Perimeter"), size=1, alpha=1/5) +
		# geom_point(aes(x=score_AR_synapse, y=score_AR_noise, color="AR"), size=1, alpha=1/5) +
		geom_point(aes(x=score_shape_synapse, y=score_shape_noise, color="shape")) +
		# geom_point(aes(x=score_Width_synapse, y=score_Width_noise, color="Width")) +
		# geom_point(aes(x=score_Height_synapse, y=score_Height_noise, color="Height")) +
		# geom_point(aes(x=score_Major_synapse, y=score_Major_noise, color="Major")) +
		# geom_point(aes(x=score_Minor_synapse, y=score_Minor_noise, color="Minor")) +
		#geom_point(aes(x=score_Round_synapse, y=score_Round_noise, color="Round"), size=1, alpha=1/5) +
		# geom_point(aes(x=score_Circularity_synapse, y=score_Circularity_noise, color="Circularity"), size=1, alpha=1/5) +
		# geom_point(aes(x=score_Solidity_synapse, y=score_Solidity_noise, color="Solidity"), size=1, alpha=1/5) +
		geom_point(aes(x=score_Feret_synapse, y=score_Feret_noise, color="Feret"), size=1, alpha=1/5) +
		#geom_point(aes(x=score_MinFeret_synapse, y=score_MinFeret_noise, color="MinFeret"), size=1, alpha=1/5) +
		#geom_point(aes(x=score_FeretAngle_synapse, y=score_FeretAngle_noise, color="FeretAngle"), size=1, alpha=1/5) +
		scale_x_log10()+
		scale_y_log10()+
		xlab("synapse") + ylab("noise") +
		geom_abline(intercept=0, slope=1, color="grey60")+
		facet_grid(.~group)+
		basic_style
	ggsave("Score scatter plot (morph).pdf")
  
	s_geometry_scatter <- ggplot(data_scatter_plot) +
		geom_point(aes(x=score_dist_synapse, y=score_dist_noise, color="Distance"), size=1, alpha=1/5) +
		geom_point(aes(x=score_angle_synapse, y=score_angle_noise, color="Angle"), size=1, alpha=1/5) +
		geom_point(aes(x=score_Area_synapse, y=score_Area_noise, color="Area"), size=1, alpha=1/5) +
		#geom_point(aes(x=score_density_synapse, y=score_density_noise, color="density"), size=1, alpha=1/5) +
		#geom_point(aes(x=score_IntDen_synapse, y=score_IntDen_noise, color="IntDen")) +
		# geom_point(aes(x=score_IntDen_bg_synapse, y=score_IntDen_bg_noise, color="IntDen_bg"), size=1, alpha=1/5) +
		#geom_point(aes(x=score_RawIntDen_synapse, y=score_RawIntDen_noise, color="RawIntDen")) +
		geom_point(aes(x=score_RawIntDen_bg_synapse, y=score_RawIntDen_bg_noise, color="RawIntDen_bg"), size=1, alpha=1/5) +
		scale_x_log10()+
		scale_y_log10()+
		xlab("synapse") + ylab("noise") +
		geom_abline(intercept=0, slope=1, color="grey60")+
		facet_grid(.~group)
		basic_style
	ggsave("Score scatter plot (geometry).pdf")
  
	s_bayes_scatter <- ggplot(data_pair_scored %>% 
							  filter(naive_bayes_synapse > 1e-300 & naive_bayes_noise > 1e-300), 
							aes(x=naive_bayes_synapse, y=naive_bayes_noise, color="naive_bayes"))+
		geom_point()+
		geom_point(aes(x=score_total_synapse, y=score_total_noise, color="score"))+
		geom_segment(aes(x=naive_bayes_synapse, y=naive_bayes_noise, 
						 xend=score_total_synapse, yend=score_total_noise), color="gray")+
		#geom_text(aes(label=pair_id))+
		#geom_text(aes(x=score_total_synapse, y=score_total_noise, label=pair_id), color="gray") +
		scale_x_log10()+
		scale_y_log10()+
		geom_abline(intercept=0, slope=1, color="grey60")+
		scale_color_manual(values = c(score = "blue", naive_bayes = "red"))+
		basic_style
	ggsave("Naive Bayes scatter plot.pdf")
}

#map score to all data so it can be plotted together
data$all <- TRUE
# data$pair <- data$X %in% coord_subset1[1,] & data$Y %in% coord_subset1[2,] | 
#                     data$X %in% coord_subset1[4,] & data$Y %in% coord_subset1[5,]

data$pair <- data$marker_id %in% data_pair$marker_id

data$geometry <- data$marker_id %in% 
	(data_pair_scored %>% filter(score_geometry_synapse  > score_geometry_noise) %>% .$marker_id)

data$signal <- data$marker_id %in% 
	(data_pair_scored %>% filter(score_signal_synapse > score_signal_noise) %>% .$marker_id)

data$morph <- data$marker_id %in% 
	(data_pair_scored %>% filter(score_morph_synapse > score_morph_noise) %>% .$marker_id)

data$score_total <- data$marker_id %in% 
	(data_pair_scored %>% filter(score_total_synapse > score_total_noise) %>% .$marker_id)

data$naive_bayes <- data$marker_id %in% 
	(data_pair_scored %>% filter(naive_bayes_synapse > naive_bayes_noise) %>% .$marker_id)


#tidy up data
data_synapse <- data %>% 
  filter(group != "DAPI") %>%
  gather(all:naive_bayes, key="evaluation", value="logic", factor_key=TRUE) %>% 
  filter(logic) #seelct only the ones that are TRUE

gen_log(paste("Number of synapses (pair)", data$pair %>% sum/2), logfile)
gen_log(paste("Number of synapses (geometry)", data$geometry %>% sum/2), logfile)
gen_log(paste("Number of synapses (signal)", data$signal %>% sum/2), logfile)
gen_log(paste("Number of synapses (morph)", data$morph %>% sum/2), logfile)
gen_log(paste("Number of synapses (Score)", data$score_total %>% sum/2), logfile)
gen_log(paste("Number of synapses (Bayes)", data$naive_bayes %>% sum/2), logfile)


gen_log("Plot synapse coord", logfile)
if (nrow(data_pair_scored)){
	p_synapse <- ggplot(data_synapse %>% filter(evaluation != "all" & evaluation !="pair")) +
		geom_point(data=(data_synapse %>% filter(evaluation=="all") %>% dplyr::select(-evaluation)), aes(x=X, y=Y), color = "gray100", size=1) + #all marker points
		geom_point(data=(data_synapse %>% filter(evaluation=="pair") %>% dplyr::select(-evaluation)), aes(x=X, y=Y), color = "gray70", size=0.6) + #all marker points
		geom_point(aes(x=X, y=Y, color=evaluation), size=0.2) +
		facet_grid(evaluation~group)+
		{if (image_3D) labs(title = "Synapse Marker (z-project)") 
		  else labs(title = "Synapse Marker")} +
		scale_color_manual(values = c("blue", "green", "pink", "yellow", "red")) +
		xlab("x(um)") + ylab("y(um)") +
		scale_y_reverse() +
		basic_style 
  		#theme(legend.position="none")
	ggsave("synapse_coord.pdf")
  
	gen_log("Marker stats plots", logfile, datetime=TRUE)
	if (image_3D) {
		p_ROI_vs_Y <- ggplot(data=data_synapse, aes(x=Area, y=Center.of.Image.Mass.Y, color=evaluation))+
			geom_point(size=0.1)+facet_wrap(~group, scales = "free_x")+
			scale_y_reverse()+
			ylab("y(um)") +
			scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
			basic_style + theme(legend.position="none")
	}else{
		p_ROI_vs_Y <- ggplot(data=data_synapse, aes(x=Area, y=YM, color=evaluation))+
		geom_point(size=0.1)+facet_wrap(~group, scales = "free_x")+
		scale_y_reverse()+
		ylab("y(um)") +
		scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
		basic_style + theme(legend.position="none")
	}
  
	h_ROI <- ggplot(data=data_synapse, aes(x=Area, fill=evaluation))+
		geom_histogram(bins=60, position="identity", alpha=1) + #geom_density()+
		facet_wrap(~group, scales = "free")+ 
		scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
		basic_style + theme(legend.position="none")
}

# features <- colnames(data_pair[sapply(data_pair,is.double)])                         #all the stats   
# #remove some unecessary stats
# features <- features[!features %in% c("X", "Y", "XM","YM", "BX", "BY", "bg", 
#                                 "X.Area", "MinThr", "MaxThr", "PercentArea", 
#                                 "FeretX", "FeretY")]

# #pair_dis, pair_angle, density are not in data
# features <- features[!features %in% c("pair_dist", "pair_angle", "density", "shape")]
# print(features)
# #make the individual stats plots
# for (feature in features){
#   assign(paste0("h_ROI_", feature), h_ROI %+% aes_string(x=feature) + labs(title = paste("ROI", feature, "Histogram")))
#   assign(paste0("p_ROI_vs_Y_", feature), p_ROI_vs_Y %+% aes_string(x=feature) + labs(title = paste("ROI vs Y coord:", feature)))
# }

# #make multiple plots with common legedn using grid.arrange
# myplots <- list()  # new empty list
# for (i in 1:length(features)){
#   myplots[[i]] <- eval(parse(text=paste0("h_ROI_", features[i])))
# }

# if (image_3D) {
#   g <- grid_arrange_shared_legend(myplots, ncol = 2, nrow = 5)
#   ggsave("ROI sats histograms.pdf", plot=g, width = 15, height = 15)
# }else {
#   g1 <- grid_arrange_shared_legend(myplots[1:20], ncol = 4, nrow = 5)
#   ggsave("ROI sats histograms1.pdf", plot=g1, width = 15, height = 15)
#   g2 <- grid_arrange_shared_legend(myplots[21:length(myplots)], ncol = 4, nrow = 5)
#   ggsave("ROI sats histograms2.pdf", plot=g2, width = 15, height = 15)
# }


# #make multiple plots with common legned using grid.arrange
# myplots <- list()  # new empty list
# for (i in 1:length(features)){
#   myplots[[i]] <- eval(parse(text=paste0("p_ROI_vs_Y_", features[i])))
# }
# if (image_3D) {
#   g <- grid_arrange_shared_legend(myplots, ncol = 2, nrow = 5)
#   ggsave("ROI vs Y coord sats.pdf", plot=g, width = 15, height = 15)
# }else {
#   g1 <- grid_arrange_shared_legend(myplots[1:20], ncol = 4, nrow = 5)
#   ggsave("ROI vs Y coord sats1.pdf", plot=g1, width = 15, height = 15)
#   g2 <- grid_arrange_shared_legend(myplots[21:length(myplots)], ncol = 4, nrow = 5)
#   ggsave("ROI vs Y coord sats2.pdf", plot=g2, width = 15, height = 15)
# }


gen_log("Calculaing synapse coordinates", logfile)
synapse <- data_pair_scored %>% 
	filter(score_total_synapse * prior_synapse > score_total_noise * (1 - prior_synapse)) 
if (dim(synapse)[1]){
	synapse <- synapse %>%
		dplyr::select(group, XM, YM, FeretAngle, pair_id) %>%
		unite("X_Y_FeretAngle", XM:FeretAngle) %>% 
		spread(group, X_Y_FeretAngle) %>%
		separate(Pre, c("Pre.X", "Pre.Y", "Pre.FeretAngle"), sep="_", convert=TRUE) %>%
		separate(Post, c("Post.X", "Post.Y", "Post.FeretAngle"), sep="_", convert=TRUE) %>%
		mutate(Syn.X = (Pre.X + Post.X)/2) %>% 
		mutate(Syn.Y = (Pre.Y + Post.Y)/2) %>% 
		unite("Pre", Pre.X:Pre.FeretAngle) %>% 
		unite("Post", Post.X:Post.FeretAngle) %>% 
		unite("synapse", Syn.X:Syn.Y) %>% 
		gather(key=group, value=X_Y_FeretAngle, Pre:synapse) %>% #print(width=Inf, n=20) %>%
		separate(X_Y_FeretAngle, c("X", "Y", "FeretAngle"), sep="_", convert=TRUE) 
}
if (image_3D){
  #3D image
} else {
	#2d image
	write.csv(synapse %>% filter(group=="Pre")  %>% dplyr::select(X:FeretAngle), file = "synapse_pre.csv", quote = FALSE)
	write.csv(synapse %>% filter(group=="Post") %>% dplyr::select(X:FeretAngle), file = "synapse_post.csv", quote = FALSE)
	write.csv(synapse %>% filter(group=="synapse") %>% dplyr::select(X:Y), file = "synapse.csv", quote = FALSE)
}


saveRDS(data_pair_scored, file = "data_pair_scored.rds")  #subset of data with pairs within synapse distance
saveRDS(data_synapse, file = "data_synapse.rds")                #tibble to compare stats
saveRDS(synapse, file = "synapse.rds")

#plot synapse coord

gen_log("Done", logfile, datetime=TRUE)