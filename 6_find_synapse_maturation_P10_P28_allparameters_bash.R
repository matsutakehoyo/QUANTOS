
library(MASS)                     #for pi breaks 2dkdes
library(grid)					  #grid_arrange_shared_legend
library(gridExtra)                #grid_arrange_shared_legend
library(propagate)                #for distDistr
library(ks)
library(plyr)                     #for mapvalues to replace NaN
library(tidyverse)                #ggplot, tibble, tidyr, readr, purrr, dplyr
library(bde)
library(BEST)
library(RColorBrewer)

dist_cutoff_min = 0
dist_cutoff_max = 1.2
distance_density = 5                                 #distance for marker density calculation

#type of image 3D or 2D
image_3D = FALSE


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


#----------------------Loading files----------------------
data <- readRDS("data_bg.rds")
data_pair <- readRDS("data_pair_scored.rds") %>% dplyr::filter(naive_bayes_synapse > naive_bayes_noise)

#split dirname and filename
split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

#Function for reading files including subdirectories
read_files <- function(filename){
  #Getting list of files
  files <- list.files(pattern=paste0(filename, ".rds"), recursive=TRUE)
  print(files)
  #Reading rds file and Adding first line to tibble
  d <- readRDS(files[1])
  name <- split_path(files[1])[2]
  d <- d %>% add_column(name = "")
  if (nrow(d)){
    d$name <- name  
  }
  
  #assigning tibble to dataframe
  assign(filename, d, envir=globalenv())
  #Adding 2nd and other lines to dataframe
  if (length(files)>1){
    for (i in seq(2, length(files))){
      d <- readRDS(files[i])
      name <- split_path(files[i])[2]
      d <- d %>% add_column(name = "")
      if (nrow(d)){
        d$name <- name  
      }
      assign(filename, bind_rows(get(filename), d), envir=globalenv())
    } 
  }
}

  
  #log R output ot file: moified original file by Sunagawa gen-log.R
  gen_start_log<-function(name="Debug"){
    dir.create("log_mature", recursive=TRUE, showWarnings=FALSE) # make the folders before saving data
    local_log_file<-paste(format(Sys.time(),"%Y%m%d%H%M%S"), name, "log", sep=".")
    return(local_log_file)
  }
  
  gen_log<-function(text, log_file, datetime=FALSE){
    dir.create(file.path(".", "log_mature"), recursive=TRUE, showWarnings=FALSE)
    if (datetime) {
      write(paste(Sys.time(), text, sep=": "), file=file.path(".","log",log_file), append=TRUE)
      print(paste(Sys.time(), text, sep=": "))
    }
    else {
      write(text, file=file.path(".", "log_mature",log_file), append=TRUE)
      print(text)
    }
  }
  
  #log results in file
  logfile <- gen_start_log("FindMatureSynapses")
  gen_log("Applying Niave Bayes Classifier to paired mature synapse markers", logfile, datetime=TRUE)
    
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
      n <- data %>% dplyr::filter(group==group_str) %>% nrow()
      gen_log(paste("Number of", group_str, "markers:", n), logfile)	
      x_min <- data %>% dplyr::filter(group==group_str) %>%dplyr::select(X) %>% min()
      x_max <- data %>% dplyr::filter(group==group_str) %>%dplyr::select(X) %>% max()
      y_min <- data %>% dplyr::filter(group==group_str) %>%dplyr::select(Y) %>% min()
      y_max <- data %>% dplyr::filter(group==group_str) %>%dplyr::select(Y) %>% max()
      
      n / ((x_max - x_min) * (y_max - y_min))
    }
    
    dPre <- calc_global_density("Pre")
    gen_log(paste("Pre density:", dPre), logfile)
    dPost <- calc_global_density("Post")
    gen_log(paste("Post density:", dPost), logfile)
    
    #read feature pdf
    read_pdf <- function(filename){
      files <- list.files(path = "./pdf_mature/", pattern=paste0(filename, "*.rds"), recursive=TRUE)
      #print(files) 
      for (i in seq_along(files)){
    	print(paste("Loading file:", files[i] ))
    	pdf <- readRDS(paste0("pdf_mature/", files[i]))
    	assign(gsub(pattern=".rds|pdf_mature/", "", x=files[i]), pdf, envir=globalenv())
    	print(gsub(pattern=".rds|pdf_matur/", "", x=files[i]))
      }
      files <- gsub(pattern=".rds|pdf_mature/", "", x=files)
    }
    read_pdf("pdf.mature")
    pdf_immature <- dplyr::combine(read_pdf("bde.immature"), read_pdf("kde.immature"))
    pdf_mature <- dplyr::combine(read_pdf("bde.mature"), read_pdf("kde.mature"))
    
    if (length(pdf_immature)!=length(pdf_mature)){
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
    
    #plot mature & immature pdf
    dist_pdf <- ggplot(data = data.frame(x = c(dist_cutoff_min, dist_cutoff_max)), aes(x), fill=eval) +
      stat_function(fun = dnorm, n = 200, 
    				args = list(mean = dist.pdf.mature$estimate[1],
    							sd=dist.pdf.mature$estimate[2]),
    				geom="area", alpha = 1/3, aes(fill="mature")) +
      stat_function(fun = function(x){3/dist_cutoff_max^3*x^2}, geom="area", alpha = 1/3, aes(fill="immature")) +
      xlab("distance (um)") + ylab("probability") + ggtitle("Distance") +
      scale_fill_manual("Marker", values = c("mature" ="#56B4E9","immature" = "#E69F00")) +
      basic_style
    if (length(data_pair)){
      dist_pdf <- dist_pdf + 
    	geom_histogram(data=(data_pair), aes(x = pair_dist, y = ..density..), alpha=1/3) 
    }  
    
    angle_pdf <- ggplot(data = data.frame(x = c(-pi, pi)), aes(x)) +
      stat_function(fun = propagate:::dlaplace, n = 200, args = list(mean = 0, sigma=pi/2), geom="area", alpha = 1/3, aes(fill="mature")) +
      stat_function(fun = dunif, n = 200, args = list(min=-pi, max=pi), geom="area", alpha = 1/3, fill="#E69F00") + 
      xlab("angle") + ylab("probability") + ggtitle("Angle")+
      scale_fill_manual("Marker", values = c("mature" ="#56B4E9","immature" = "#E69F00")) +
      basic_style 
    if (length(data_pair)){
      angle_pdf <- angle_pdf +
    	geom_histogram(data=(data_pair), aes(x = pair_angle, y = ..density..), alpha=1/3) 
    }
    #angle_pdf <- fracAx(angle_pdf, "pi")
    
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
    
    #function to plot the pdf for immature and mature pairs
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
    	scale_fill_manual(values = c("mature" ="#56B4E9","immature" = "#E69F00")) +
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
    	  stat_function(fun = bde::density, n = 200, args = list(x = eval(parse(text=paste0(feature, ".", group_name, ".", pdf_type, ".mature")))), 
    					geom="area", alpha = 1/3, aes(fill="mature")) +    
    	  stat_function(fun = bde::density, n = 200, args = list(x = eval(parse(text=paste0(feature, ".", group_name, ".", pdf_type, ".immature")))), 
    					geom="area", alpha = 1/3, aes(fill="immature")) 
      } else if (pdf_type=="kde"){
    	p <- p +
    	  stat_function(fun = dkde, n = 200, args = list(fhat = eval(parse(text=paste0(feature, ".", group_name, ".", pdf_type, ".mature")))), 
    					geom="area", alpha = 1/3, aes(fill="mature")) +    
    	  stat_function(fun = dkde, n = 200, args = list(fhat = eval(parse(text=paste0(feature, ".", group_name, ".", pdf_type, ".immature")))), 
    					geom="area", alpha = 1/3, aes(fill="immature")) 
      }
    }
    
    for (pdf in pdf_mature){
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
    ggsave(filename ="ROI sats probability distribution (geometry1)_Mature.pdf", plot=h_pdf_geometry1 , width = 15, height = 9)
    
    h_pdf_geometry2 <- grid_arrange_shared_legend(pdf_geometry[7:18], ncol = 3, nrow = 4)
    ggsave(filename = "ROI sats probability distribution (geometry2)_Mature.pdf", plot=h_pdf_geometry2 , width = 15, height = 10)
    
    
    h_pdf_signal1 <- grid_arrange_shared_legend(pdf_signal[1:18], ncol = 3, nrow = 6)
    ggsave(filename = "ROI sats probability distribution (signal1)_Mature.pdf", plot=h_pdf_signal1 , width = 15, height = 10)
    
    h_pdf_signal2 <- grid_arrange_shared_legend(pdf_signal[19:36], ncol = 3, nrow = 6)
    ggsave(filename = "ROI sats probability distribution (signal2)_Mature.pdf", plot=h_pdf_signal2 , width = 15, height = 10)
    
    h_pdf_signal3 <- grid_arrange_shared_legend(pdf_signal[37:40], ncol = 2, nrow = 2)
    ggsave(filename = "ROI sats probability distribution (signal3)_Mature.pdf", plot=h_pdf_signal3 , width = 15, height = 6)
    
    h_pdf_morphology <- grid_arrange_shared_legend(pdf_morphology, ncol = 4, nrow = 7)
    ggsave(filename = "ROI sats probability distribution (shape)_Mature.pdf", plot=h_pdf_morphology , width = 15, height = 10)
    
    
    #function to score each feature according to its pdf, subset the data to the group 
    calc_score <- function(pdf){
      score_name <- gsub(pattern=".(Pre|Post).(bde|kde).(mature|immature)", replacement="", pdf)
      group_name <- pdf %>% strsplit("[.]") %>% unlist() %>% .[2] %>% firstup(.) #pre or post
      pdf_type <- pdf %>% strsplit("[.]") %>% unlist() %>% .[3] #bde or kde
      marker_type <- pdf %>% strsplit("[.]") %>% unlist() %>% .[4] #immature or mature
      
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
    
    pdf_list <- dplyr::combine(pdf_mature, pdf_immature)
    
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
    
    #calculate mature probability
    data_pair_scored <- data_pair %>% 
      #distance
      #mutate(score_dist_mature = dnorm(pair_dist, 
    	#								mean = dist.pdf.mature$estimate[1], 
    	#								sd = dist.pdf.mature$estimate[2])) %>%
      #mutate(score_dist_immature = 3 / dist_cutoff_max^3 * pair_dist^2) %>%
      mutate(score_dist_mature = 0.5) %>%
      mutate(score_dist_immature = 0.5) %>%
      #angle
      #mutate(score_angle_mature = dnorm(pair_angle, mean = 0, sd = pi/4)) %>%
      #mutate(score_angle_mature = propagate:::dlaplace(pair_angle, mean = 0, sigma = pi/2)) %>%
      #mutate(score_angle_immature = dunif(pair_angle, min= -pi, max= pi)) 
      mutate(score_angle_mature = 0.5) %>%
      mutate(score_angle_immature = 0.5) 
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
    		gsub(pattern = "_mature|_immature", "", .)
    		signal_name <- score %>% strsplit("[_]") %>% unlist() %>% .[length(.)]
    		gen_log(paste(zero, "zeros in score:", feature_name, signal_name), logfile)
    		zero_rows <- which(data_pair_scored %>% dplyr::select(score) ==0) %>% print
    		
    		if (signal_name == "mature"){
    		  opposite_score <- paste0("score_", feature_name, "_immature")
    		}else if (signal_name == "immature"){
    		  opposite_score <- paste0("score_", feature_name, "_mature")
    		}
    		data_pair_scored[zero_rows,score] <-0.5
    		data_pair_scored[zero_rows,opposite_score] <-0.5
    	}
    	  
    	na <- data_pair_scored %>% dplyr::select(score) %>% is.na() %>% sum()
    	if (na){
    		feature_name <- gsub(pattern = "score_", "", score) %>% 
    		  gsub(pattern = "_mature|_immature", "", .)
    		signal_name <- score %>% strsplit("[_]") %>% unlist() %>% .[length(.)]
    		gen_log(paste(na, "NA in score:", feature_name, signal_name), logfile)
    		na_rows <- which(data_pair_scored %>% dplyr::select(score) %>% is.na()) %>% print
    		
    		if (signal_name == "mature"){
    		  opposite_score <- paste0("score_", feature_name, "_immature")
    		}else if (signal_name == "immature"){
    		  opposite_score <- paste0("score_", feature_name, "_mature")
    		}
    		data_pair_scored[na_rows,score] <-0.5
    		data_pair_scored[na_rows,opposite_score] <-0.5
    	}
    	  
    	nan <- data_pair_scored %>% dplyr::select(score) %>% pull %>% is.nan() %>% sum()
    	  
    	if (nan){
    		feature_name <- gsub(pattern = "score_", "", score) %>% 
    		  gsub(pattern = "_mature|_immature", "", .)
    		signal_name <- score %>% strsplit("[_]") %>% unlist() %>% .[length(.)]
    		gen_log(paste(na, "NaN in score:", feature_name, signal_name), logfile)
    		nan_rows <- which(data_pair_scored %>% dplyr::select(score) %>% pull() %>% is.nan()) %>% print
    		
    		if (signal_name == "mature"){
    		  opposite_score <- paste0("score_", feature_name, "_immature")
    		}else if (signal_name == "immature"){
    		  opposite_score <- paste0("score_", feature_name, "_mature")
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
    	dplyr::select(ends_with("_mature")) %>% 
    	dplyr::select(
    	  score_dist_mature, 
    	  score_angle_mature, 
    	  score_density_mature,
    	  score_Area_mature, 
    	  score_IntDen_bg_mature,
    	  score_RawIntDen_bg_mature
    				) %>%
    	transmute(score_geometry_mature = apply(., 1, prod)) %>% 
    	bind_cols(data_pair_scored_pre, .) 
    
    data_pair_scored_pre <- data_pair_scored_pre %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_immature")) %>%
    	dplyr::select(
    	  score_dist_immature, 
    	  score_angle_immature,
    	  score_density_immature,
    	  score_Area_immature,
    	  score_IntDen_bg_immature,
    	  score_RawIntDen_bg_immature
    				) %>%
    	transmute(score_geometry_immature = apply(., 1, prod)) %>%
    	bind_cols(data_pair_scored_pre, .) 
    
    data_pair_scored_post <- data_pair_scored %>% 
    	filter(group=="Post") 
    data_pair_scored_post <- data_pair_scored_post %>%
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_mature")) %>% 
    	dplyr::select(
    	  score_dist_mature,
    	  score_angle_mature,
    	  score_density_mature,
    	  score_Area_mature,
    	  score_IntDen_bg_mature,
    	  score_RawIntDen_bg_mature
    	) %>%
    	transmute(score_geometry_mature = apply(., 1, prod)) %>%
    	bind_cols(data_pair_scored_post, .) 
    
    data_pair_scored_post <- data_pair_scored_post %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_immature")) %>%
    	dplyr::select(
    	  score_dist_immature, 
    	  score_angle_immature,
    	  score_density_immature,
    	  score_Area_immature,
    	  score_IntDen_bg_immature,
    	  score_RawIntDen_bg_immature
    				) %>%
    	transmute(score_geometry_immature = apply(., 1, prod)) %>%
    	bind_cols(data_pair_scored_post, .)
    
    #signal score
    data_pair_scored_pre <- data_pair_scored_pre %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_mature")) %>% 
    	dplyr::select(
    				   score_Mean_mature,
    				   score_Median_mature,
    				   score_Mode_mature,
    				   score_Mean_bg_mature,
    				   score_Median_bg_mature,
    				   score_Mode_bg_mature,
    				   score_Min_mature,
    				   score_StdDev_bg_mature,
    				   score_Max_mature,
    				   score_Max_bg_mature,
    				   score_Min_bg_mature,
    				   score_Skew_mature,
    				   score_Kurt_mature 
      				) %>%
    	transmute(score_signal_mature = apply(., 1, prod)) %>% 
    	bind_cols(data_pair_scored_pre, .) 
    
    data_pair_scored_pre <- data_pair_scored_pre %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_immature")) %>%
    	dplyr::select(
          	  score_Mean_immature,
          	  score_Median_immature,
          	  score_Mode_immature,
          	  score_Mean_bg_immature,
          	  score_Median_bg_immature,
          	  score_Mode_bg_immature,
          	  score_Min_immature,
          	  score_StdDev_bg_immature,
          	  score_Max_immature,
          	  score_Max_bg_immature,
          	  score_Min_bg_immature,
          	  score_Skew_immature,
          	  score_Kurt_immature 
      				 ) %>%
    	transmute(score_signal_immature = apply(., 1, prod)) %>% 
    	bind_cols(data_pair_scored_pre, .)
    
    data_pair_scored_post <- data_pair_scored_post %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_mature")) %>% 
    	dplyr::select(
    	        score_Mean_mature,
          	  score_Median_mature,
          	  score_Mode_mature,
          	  score_Mean_bg_mature,
          	  score_Median_bg_mature,
          	  score_Mode_bg_mature,
          	  score_Min_mature,
          	  score_StdDev_bg_mature,
          	  score_Max_mature,
          	  score_Max_bg_mature,
          	  score_Min_bg_mature,
          	  score_Skew_mature,
          	  score_Kurt_mature
    	) %>%
    	transmute(score_signal_mature = apply(., 1, prod)) %>%
    	bind_cols(data_pair_scored_post, .)
    
    data_pair_scored_post <- data_pair_scored_post %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_immature")) %>%
    	dplyr::select(
          	  score_Mean_immature,
          	  score_Median_immature,
          	  score_Mode_immature,
          	  score_Mean_bg_immature,
          	  score_Median_bg_immature,
          	  score_Mode_bg_immature,
          	  score_Min_immature,
          	  score_StdDev_bg_immature,
          	  score_Max_immature,
          	  score_Max_bg_immature,
          	  score_Min_bg_immature,
          	  score_Skew_immature,
          	  score_Kurt_immature 
    	 ) %>%
    	transmute(score_signal_immature = apply(., 1, prod)) %>% 
    	bind_cols(data_pair_scored_post, .)
    
    #morphology score
    data_pair_scored_pre <- data_pair_scored_pre %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_mature")) %>% 
    	dplyr::select(
    				 score_Perimeter_mature,
    				 score_shape_mature,
    				 score_Width_mature,
    				 score_Height_mature,			#bounding box
    				 score_Major_mature,
    				 score_Minor_mature,
    				 score_Angle_mature,
    				 score_AR_mature,
    				 score_Round_mature,
    				 score_Circularity_mature,
    				 score_Solidity_mature,
    				 score_Feret_mature,
    				 score_MinFeret_mature,
    		     	 score_FeretAngle_mature
    				 ) %>%
      transmute(score_morph_mature = apply(., 1, prod)) %>% 
      bind_cols(data_pair_scored_pre, .)
    
    data_pair_scored_pre <- data_pair_scored_pre %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_immature")) %>%
    	dplyr::select(
        	  score_Perimeter_immature,
        	  score_shape_immature,
        	  score_Width_immature,
        	  score_Height_immature,			#bounding box
        	  score_Major_immature,
        	  score_Minor_immature,
        	  score_Angle_immature,
        	  score_AR_immature,
        	  score_Round_immature,
        	  score_Circularity_immature,
        	  score_Solidity_immature,
        	  score_Feret_immature,
        	  score_MinFeret_immature,
        	  score_FeretAngle_immature
      				 ) %>%
    	transmute(score_morph_immature = apply(., 1, prod)) %>% 
    	bind_cols(data_pair_scored_pre, .)
    
    data_pair_scored_post <- data_pair_scored_post %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_mature")) %>% 
    	dplyr::select(
        	  score_Perimeter_mature,
        	  score_shape_mature,
        	  score_Width_mature,
        	  score_Height_mature,			#bounding box
        	  score_Major_mature,
        	  score_Minor_mature,
        	  score_Angle_mature,
        	  score_AR_mature,
        	  score_Round_mature,
        	  score_Circularity_mature,
        	  score_Solidity_mature,
        	  score_Feret_mature,
        	  score_MinFeret_mature,
        	  score_FeretAngle_mature
    	) %>%
    	transmute(score_morph_mature = apply(., 1, prod)) %>% 
    	bind_cols(data_pair_scored_post, .)
    
    data_pair_scored_post <- data_pair_scored_post %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_immature")) %>%
    	dplyr::select(
          	  score_Perimeter_immature,
          	  score_shape_immature,
          	  score_Width_immature,
          	  score_Height_immature,			#bounding box
          	  score_Major_immature,
          	  score_Minor_immature,
          	  score_Angle_immature,
          	  score_AR_immature,
          	  score_Round_immature,
          	  score_Circularity_immature,
          	  score_Solidity_immature,
          	  score_Feret_immature,
          	  score_MinFeret_immature,
          	  score_FeretAngle_immature
    	) %>%
    	transmute(score_morph_immature = apply(., 1, prod)) %>% 
    	bind_cols(data_pair_scored_post, .)
    
    data_pair_scored <- bind_rows(data_pair_scored_pre, data_pair_scored_post)
    
    #total score
    data_pair_scored <- data_pair_scored %>%  
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_mature")) %>% 
    	dplyr::select(score_geometry_mature, score_signal_mature, score_morph_mature) %>%
    	transmute(score_all_mature = apply(., 1, prod)) %>%
    	bind_cols(data_pair_scored, .)
    
    data_pair_scored <- data_pair_scored %>% 
    	dplyr::select(starts_with("score_")) %>% 
    	dplyr::select(ends_with("_immature")) %>%
    	dplyr::select(score_geometry_immature, score_signal_immature, score_morph_immature) %>%
    	transmute(score_all_immature = apply(., 1, prod)) %>%
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
    	dplyr::select(pair_id, marker_id, group,  ends_with("_mature")) %>%
    	gather(feature, score, matches("score"), factor_key = TRUE) %>% 
    	filter(feature!="score_all_mature") %>%
    	separate(feature, c("a", "feature", "eval"), sep = c(6,-7)) %>%  #cannnot separate at "_" because _bg
    	dplyr::select(-a) %>% 
    	mutate(eval="mature")
    
    marker_features_noi <- data_pair_scored %>% 
    	dplyr::select(pair_id, marker_id, group,  ends_with("_immature")) %>%
    	gather(feature, score, matches("score"), factor_key = TRUE) %>%
    	filter(feature!="score_all_immature") %>%
    	separate(feature, c("a", "feature", "eval"), sep = c(6, -9)) %>% 
    	dplyr::select(-a) %>% 
    	mutate(eval="immature")
    
    marker_features <- bind_rows(marker_features_syn, marker_features_noi) %>%
    	filter(score>1e-300) #this avoid errors when plotting
    
    s_marker_features <- ggplot(marker_features) + basic_style
    if (nrow(marker_features)){
    	s_marker_features <- s_marker_features + 
    		geom_point(data=marker_features %>% 
    						spread(eval, score) %>% 
    						filter(feature!="geometry" & feature!="signal" & feature!="morph"), 
    				aes(x=mature, y=immature, color=group), size=0.1) +
    	geom_abline(intercept=0, slope=1, color="darkgray") +
    	scale_x_log10()+
    	scale_y_log10()+
    	  facet_wrap(~feature, scales="free", ncol=6) 
    }
    ggsave(filename = "mature-immature score scatter.pdf", width=15, height=10)
    
    s_marker_features_summary <- ggplot() + basic_style
    if (nrow(marker_features)){
    	s_marker_features_summary <- s_marker_features_summary + 
    	geom_point(data = marker_features %>% 
    				spread(eval, score) %>% 
    				filter(feature=="geometry" | feature=="signal" | feature=="morph"), 
    			aes(x=mature, y=immature, color=group), size=0.1) +
    	geom_abline(intercept=0, slope=1, color="darkgray") +
    	scale_x_log10()+
    	scale_y_log10()+
    	facet_wrap(~feature, scales="free", ncol=6) 
    }
    ggsave(filename =  "mature-immature score scatter summary.pdf", width=15, height=10)
    
    #calculate probability of random/mature
    # poly(dist.max, 2, raw = TRUE)1 poly(dist.max, 2, raw = TRUE)2 
    #                        877.599                       7136.395 
    rand.pair <- function(dist, dens1, dens2){
    	(877.599 * dist + 7136.395 * dist^2) * dens1 * dens2
    }
    
    mature_prior <- function(dist, dens1, dens2){
    	ifelse(rand.pair(dist, dens1, dens2) < 1, 0.5, 1/rand.pair(dist, dens1, dens2))
    }
    
    #calculate Pre Post combined score
    if (nrow(data_pair_scored)){
    	data_pair_scored <- data_pair_scored %>% 
    		dplyr::select(pair_id, group, score_all_mature) %>%
    		spread(group, score_all_mature) %>%
    		mutate(score_total_mature = Pre * Post) %>%
    		gather(group, score_marker, Pre:Post) %>% 
    		arrange(pair_id) %>% dplyr::select(score_total_mature) %>%
    		bind_cols(data_pair_scored %>% arrange(pair_id), .)
      
    	data_pair_scored <- data_pair_scored %>% 
    		dplyr::select(pair_id, group, score_all_immature) %>%
    		spread(group, score_all_immature) %>%
    		mutate(score_total_immature = Pre * Post) %>%
    		gather(group, score_marker, Pre:Post) %>% 
    		arrange(pair_id) %>% dplyr::select(score_total_immature) %>%
    		bind_cols(data_pair_scored %>% arrange(pair_id), .)
      
    	#calculate Pre Post combined prior
    	data_pair_scored <- data_pair_scored %>% 
    		dplyr::select(pair_id, group, pair_dist, density) %>%
    		spread(group, density) %>%
    		#mutate(prior_mature = mature_prior(pair_dist, Pre, Post)) %>%
    	  mutate(prior_mature = 0.5) %>%
    		gather(group, mature_prior, Pre:Post) %>% 
    		arrange(pair_id) %>% dplyr::select(prior_mature) %>%
    		bind_cols(data_pair_scored %>% arrange(pair_id), .)
      
      #calculate naive_bayes mature
    	data_pair_scored <- data_pair_scored %>% 
    		dplyr::select(pair_id, group, score_total_mature, prior_mature) %>%
    		mutate(naive_bayes_mature = score_total_mature * prior_mature) %>%
    		arrange(pair_id) %>% dplyr::select(naive_bayes_mature) %>%
    		bind_cols(data_pair_scored %>% arrange(pair_id), .)
      
    	data_pair_scored <- data_pair_scored %>% 
    		dplyr::select(pair_id, group, score_total_immature, prior_mature) %>%
    		mutate(naive_bayes_immature = score_total_immature * (1-prior_mature)) %>%
    		arrange(pair_id) %>% dplyr::select(naive_bayes_immature) %>%
    		bind_cols(data_pair_scored %>% arrange(pair_id), .)
    
    	gen_log(paste("All pairs:", length(data_pair_scored$pair_id %>% unique)), logfile)
    
    	#remove duplicates first remove duplicates from one side then the other to avoid orphan markers
    	data_pair_scored %>% 
    		mutate(naive_bayes = naive_bayes_mature - naive_bayes_immature) %>% 
    		dplyr::select(pair_id, marker_id, group, naive_bayes) %>%
    		spread(group, marker_id) %>% 
    		group_by(Pre) %>% 
    		filter(naive_bayes == max(naive_bayes)) %>%
    		ungroup() %>% 
    		dplyr::select(pair_id) %>% pull %>% unique() -> unique_pair
      
      	data_pair_scored <- data_pair_scored %>% filter(pair_id %in% unique_pair)
      
    	data_pair_scored %>% 
    		mutate(naive_bayes = naive_bayes_mature - naive_bayes_immature) %>% 
    		dplyr::select(pair_id, marker_id, group, naive_bayes) %>%
    		spread(group, marker_id) %>% 
    		group_by(Post) %>% 
    		filter(naive_bayes == max(naive_bayes)) %>%
    		ungroup() %>% 
    		dplyr::select(pair_id) %>% pull %>% unique() -> unique_pair
      
      	data_pair_scored <- data_pair_scored %>% filter(pair_id %in% unique_pair)
    } else {
    	data_pair_scored <- data_pair_scored %>% 
    		add_column(score_total_mature = as.double(),
    				   score_total_immature = as.double(),
    				   prior_mature = as.double(), 
    				   naive_bayes_mature = as.double(),
    				   naive_bayes_immature = as.double(), 
    				   naive_bayes = as.double())
      	unique_pair <- as.integer()
    	gen_log(paste("All pairs:", length(data_pair_scored$pair_id %>% unique)), logfile)
    }
    
    gen_log(paste("Unique pairs:", length(unique_pair)), logfile)
    
    p_mature_prior <- ggplot(data_pair_scored, aes(x=prior_mature)) + geom_histogram()
    
    if (nrow(data_pair_scored)){
    	data_scatter_plot <- data_pair_scored %>% 
    		dplyr::select(marker_id, pair_id, group, starts_with("score_")) %>% 
    		gather(score, value, matches("score_"), factor_key = TRUE) %>% 
    		filter(value>1e-300) %>%
    		spread(score, value)
      
    	s_signal_scatter <- ggplot(data_scatter_plot) +
    		# geom_point(aes(x=score_Mean_mature, y=score_Mean_immature, color="Mean"), size=1, alpha=1/5) +
    		# geom_point(aes(x=score_Mode_mature, y=score_Mode_immature, color="Mode"), size=1, alpha=1/5) +
    		# geom_point(aes(x=score_Median_mature, y=score_Median_immature, color="Median"), size=1, alpha=1/5) +  
    		geom_point(aes(x=score_Mean_bg_mature, y=score_Mean_bg_immature, color="Mean_bg"), size=1, alpha=1/5) +
    		geom_point(aes(x=score_Mode_bg_mature, y=score_Mode_bg_immature, color="Mode_bg"), size=1, alpha=1/5) +
    		geom_point(aes(x=score_Median_bg_mature, y=score_Median_bg_immature, color="Median_bg"), size=1, alpha=1/5) +
    		geom_point(aes(x=score_StdDev_bg_mature, y=score_StdDev_bg_immature, color="StdDev_bg"), size=1, alpha=1/5) +
    		#geom_point(aes(x=score_Min_bg_mature, y=score_Min_bg_immature, color="Min_bg")) +
    		geom_point(aes(x=score_Max_mature, y=score_Max_immature, color="Max"), size=1, alpha=1/5) +
    		geom_point(aes(x=score_Max_bg_mature, y=score_Max_bg_immature, color="Max_bg"), size=1, alpha=1/5) +
    		scale_x_log10()+
    		scale_y_log10()+  
    		facet_grid(.~group)+
    		xlab("mature") + ylab("immature") +
    		geom_abline(intercept=0, slope=1, color="grey60")+
    		basic_style
      	ggsave(filename = "Score scatter plot (signal)_Mature.pdf")
     
    	s_morph_scatter <- ggplot(data_scatter_plot) +
    		geom_point(aes(x=score_Perimeter_mature, y=score_Perimeter_immature, color="Perimeter"), size=1, alpha=1/5) +
    		# geom_point(aes(x=score_AR_mature, y=score_AR_immature, color="AR"), size=1, alpha=1/5) +
    		geom_point(aes(x=score_shape_mature, y=score_shape_immature, color="shape")) +
    		# geom_point(aes(x=score_Width_mature, y=score_Width_immature, color="Width")) +
    		# geom_point(aes(x=score_Height_mature, y=score_Height_immature, color="Height")) +
    		# geom_point(aes(x=score_Major_mature, y=score_Major_immature, color="Major")) +
    		# geom_point(aes(x=score_Minor_mature, y=score_Minor_immature, color="Minor")) +
    		#geom_point(aes(x=score_Round_mature, y=score_Round_immature, color="Round"), size=1, alpha=1/5) +
    		# geom_point(aes(x=score_Circularity_mature, y=score_Circularity_immature, color="Circularity"), size=1, alpha=1/5) +
    		# geom_point(aes(x=score_Solidity_mature, y=score_Solidity_immature, color="Solidity"), size=1, alpha=1/5) +
    		geom_point(aes(x=score_Feret_mature, y=score_Feret_immature, color="Feret"), size=1, alpha=1/5) +
    		#geom_point(aes(x=score_MinFeret_mature, y=score_MinFeret_immature, color="MinFeret"), size=1, alpha=1/5) +
    		#geom_point(aes(x=score_FeretAngle_mature, y=score_FeretAngle_immature, color="FeretAngle"), size=1, alpha=1/5) +
    		scale_x_log10()+
    		scale_y_log10()+
    		xlab("mature") + ylab("immature") +
    		geom_abline(intercept=0, slope=1, color="grey60")+
    		facet_grid(.~group)+
    		basic_style
    	ggsave(filename = "Score scatter plot (morph)_Mature.pdf")
    
    	s_geometry_scatter <- ggplot(data_scatter_plot) +
    		geom_point(aes(x=score_dist_mature, y=score_dist_immature, color="Distance"), size=1, alpha=1/5) +
    		geom_point(aes(x=score_angle_mature, y=score_angle_immature, color="Angle"), size=1, alpha=1/5) +
    		geom_point(aes(x=score_Area_mature, y=score_Area_immature, color="Area"), size=1, alpha=1/5) +
    		#geom_point(aes(x=score_density_mature, y=score_density_immature, color="density"), size=1, alpha=1/5) +
    		#geom_point(aes(x=score_IntDen_mature, y=score_IntDen_immature, color="IntDen")) +
    		# geom_point(aes(x=score_IntDen_bg_mature, y=score_IntDen_bg_immature, color="IntDen_bg"), size=1, alpha=1/5) +
    		#geom_point(aes(x=score_RawIntDen_mature, y=score_RawIntDen_immature, color="RawIntDen")) +
    		geom_point(aes(x=score_RawIntDen_bg_mature, y=score_RawIntDen_bg_immature, color="RawIntDen_bg"), size=1, alpha=1/5) +
    		scale_x_log10()+
    		scale_y_log10()+
    		xlab("mature") + ylab("immature") +
    		geom_abline(intercept=0, slope=1, color="grey60")+
    		facet_grid(.~group)
    		basic_style
    	ggsave(filename = "Score scatter plot (geometry)_Mature.pdf")
  
    	s_bayes_scatter <- ggplot(data_pair_scored %>% 
    							  filter(naive_bayes_mature > 1e-300 & naive_bayes_immature > 1e-300), 
    							aes(x=naive_bayes_mature, y=naive_bayes_immature, color="naive_bayes"))+
    		geom_point()+
    		geom_point(aes(x=score_total_mature, y=score_total_immature, color="score"))+
    		geom_segment(aes(x=naive_bayes_mature, y=naive_bayes_immature, 
    						 xend=score_total_mature, yend=score_total_immature), color="gray")+
    		#geom_text(aes(label=pair_id))+
    		#geom_text(aes(x=score_total_mature, y=score_total_immature, label=pair_id), color="gray") +
    		scale_x_log10()+
    		scale_y_log10()+
    		geom_abline(intercept=0, slope=1, color="grey60")+
    		scale_color_manual(values = c(score = "blue", naive_bayes = "red"))+
    		basic_style
    	ggsave(filename = "Naive Bayes scatter plot_Mature.pdf")
    }
    
    #map score to all data so it can be plotted together
    data$all <- TRUE
    # data$pair <- data$X %in% coord_subset1[1,] & data$Y %in% coord_subset1[2,] | 
    #                     data$X %in% coord_subset1[4,] & data$Y %in% coord_subset1[5,]
    
    data$pair <- data$marker_id %in% data_pair$marker_id
    
    data$geometry <- data$marker_id %in% 
    	(data_pair_scored %>% filter(score_geometry_mature  > score_geometry_immature) %>% .$marker_id)
    
    data$signal <- data$marker_id %in% 
    	(data_pair_scored %>% filter(score_signal_mature > score_signal_immature) %>% .$marker_id)
    
    data$morph <- data$marker_id %in% 
    	(data_pair_scored %>% filter(score_morph_mature > score_morph_immature) %>% .$marker_id)
    
    data$score_total <- data$marker_id %in% 
    	(data_pair_scored %>% filter(score_total_mature > score_total_immature) %>% .$marker_id)
    
    data$naive_bayes <- data$marker_id %in% 
    	(data_pair_scored %>% filter(naive_bayes_mature > naive_bayes_immature) %>% .$marker_id)
    
    
    #tidy up data
    data_mature <- data %>% 
      filter(group != "DAPI") %>%
      gather(all:naive_bayes, key="evaluation", value="logic", factor_key=TRUE) %>% 
      filter(logic)  #seelct only the ones that are TRUE
  
    
    gen_log(paste("Number of matures (pair)", data$pair %>% sum/2), logfile)
    gen_log(paste("Number of matures (geometry)", data$geometry %>% sum/2), logfile)
    gen_log(paste("Number of matures (signal)", data$signal %>% sum/2), logfile)
    gen_log(paste("Number of matures (morph)", data$morph %>% sum/2), logfile)
    gen_log(paste("Number of matures (Score)", data$score_total %>% sum/2), logfile)
    gen_log(paste("Number of matures (Bayes)", data$naive_bayes %>% sum/2), logfile)
    
    
    gen_log("Plot mature coord", logfile)
    if (nrow(data_pair_scored)){
    	p_mature <- ggplot(data_mature %>% filter(evaluation != "all" & evaluation !="pair")) +
    		geom_point(data=(data_mature %>% filter(evaluation=="all") %>% dplyr::select(-evaluation)), aes(x=X, y=Y), color = "gray100", size=1) + #all marker points
    		geom_point(data=(data_mature %>% filter(evaluation=="pair") %>% dplyr::select(-evaluation)), aes(x=X, y=Y), color = "gray70", size=0.6) + #all marker points
    		geom_point(aes(x=X, y=Y, color=evaluation), size=0.2) +
    		facet_grid(evaluation~group)+
    		{if (image_3D) labs(title = "mature Marker (z-project)") 
    		  else labs(title = "mature Marker")} +
    		scale_color_manual(values = c("blue", "green", "pink", "yellow", "red")) +
    		xlab("x(um)") + ylab("y(um)") +
    		scale_y_reverse() +
    		basic_style 
      		#theme(legend.position="none")
    	ggsave(filename = "mature_coord.pdf")
    	
    	gen_log("Marker stats plots", logfile, datetime=TRUE)
    	if (image_3D) {
    		p_ROI_vs_Y <- ggplot(data=data_mature, aes(x=Area, y=Center.of.Image.Mass.Y, color=evaluation))+
    			geom_point(size=0.1)+facet_wrap(~group, scales = "free_x")+
    			scale_y_reverse()+
    			ylab("y(um)") +
    			scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
    			basic_style + theme(legend.position="none")
    	}else{
    		p_ROI_vs_Y <- ggplot(data=data_mature, aes(x=Area, y=YM, color=evaluation))+
    		geom_point(size=0.1)+facet_wrap(~group, scales = "free_x")+
    		scale_y_reverse()+
    		ylab("y(um)") +
    		scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
    		basic_style + theme(legend.position="none")
    	}
      
    	h_ROI <- ggplot(data=data_mature, aes(x=Area, fill=evaluation))+
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
    
    
    gen_log("Calculaing mature coordinates", logfile)
    mature <- data_pair_scored %>% 
    	filter(score_total_mature * prior_mature > score_total_immature * (1 - prior_mature)) 
    if (dim(mature)[1]){
    	mature <- mature %>%
    		dplyr::select(group, XM, YM, FeretAngle, pair_id) %>%
    		unite("X_Y_FeretAngle", XM:FeretAngle) %>% 
    		spread(group, X_Y_FeretAngle) %>%
    		separate(Pre, c("Pre.X", "Pre.Y", "Pre.FeretAngle"), sep="_", convert=TRUE) %>%
    		separate(Post, c("Post.X", "Post.Y", "Post.FeretAngle"), sep="_", convert=TRUE) %>%
    		mutate(Syn.X = (Pre.X + Post.X)/2) %>% 
    		mutate(Syn.Y = (Pre.Y + Post.Y)/2) %>% 
    		unite("Pre", Pre.X:Pre.FeretAngle) %>% 
    		unite("Post", Post.X:Post.FeretAngle) %>% 
    		unite("mature", Syn.X:Syn.Y) %>% 
    		gather(key=group, value=X_Y_FeretAngle, Pre:mature) %>% #print(width=Inf, n=20) %>%
    		separate(X_Y_FeretAngle, c("X", "Y", "FeretAngle"), sep="_", convert=TRUE) 
    }
    if (image_3D){
      #3D image
    } else {
    	#2d image
    	write.csv(mature %>% filter(group=="Pre")  %>% dplyr::select(X:FeretAngle), file = "mature_pre.csv", quote = FALSE)
    	write.csv(mature %>% filter(group=="Post") %>% dplyr::select(X:FeretAngle), file = "mature_post.csv", quote = FALSE)
    	write.csv(mature %>% filter(group=="mature") %>% dplyr::select(X:Y), file = "mature.csv", quote = FALSE)
    }
    
    
    saveRDS(data_pair_scored, file = "data_pair_scored_mature.rds")  #subset of data with pairs within mature distance
    saveRDS(data_mature, file = "data_mature.rds")                #tibble to compare stats
    saveRDS(mature, file = "mature.rds")
    
    #plot mature coord
    
    gen_log("Done", logfile, datetime=TRUE)
