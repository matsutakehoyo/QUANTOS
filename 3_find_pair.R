library(MASS)                     #for pi breaks 2dkdes
library(grid)                     #grid_arrange_shared_legend
library(gridExtra)                #grid_arrange_shared_legend
library(ks)                       #hpi
library(plyr)                     #for mapvalues to replace NaN
library(tidyverse)                #ggplot, tibble, tidyr, readr, purrr, dplyr


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

#functions to calculate distance & angles
if (image_3D) {
  calcDist <- function(x){(x[1]-x[4])^2+(x[2]-x[5])^2+(x[3]-x[6])^2}
}else {
  calcDist <- function(x){(x[1]-x[4])^2+(x[2]-x[5])^2}
}
#calculate angle Y
calcAngle.y <- function(x){atan2(x[1]-x[4], x[2]-x[5])}
#calculate angle Z
if (image_3D) {
  calcAngle.z <- function(x){atan2(x[1]-x[4], x[3]-x[6])}
}else {
  calcAngle.z <- function(x){NA}
}

#log results in file
logfile <- gen_start_log("FindPair")
gen_log("Analyzing synapse markers within pair distance", logfile, datetime=TRUE)


#Distance & Angle parameters: synapse distance in µm
dist_cutoff_min = 0
dist_cutoff_max = 1.2
#distance_density = 4                                 #distance for marker density calculation
distance_density = 5

gen_log(paste("Min distance:", dist_cutoff_min, "um"), logfile)
gen_log(paste("Max distance:", dist_cutoff_max, "um"), logfile)

#read pre/post synaptic marker XY coordinates
if (image_3D) load("data3D.Rda") else load("data2D.Rda")

#reorder markers:"DAPI" "Pre"  "Post" > "Pre", "Post", "DAPI"
data$group <- factor(data$group, levels = c("Pre", "Post", "DAPI"))
data$marker_id <- seq.int(nrow(data))

print(data)

saveRDS(data, file = "data.rds")  #subset of data_bg with pairs within synapse distance
gen_log("Pair tibble created: data.rds", logfile)


p_XY_dapi <- ggplot(data=filter(data, group=="DAPI"), aes(x=X, y=Y))+
  geom_point(aes(size=Area, alpha=1/255*Mean), color="royalblue3")+
  scale_y_reverse()+
  ylab("y(um)") +
  xlab("x(um)") +
  basic_style + theme(legend.position="none")
ggsave("DAPI distribution.pdf", width = 12, height = 6)


data_pre <- filter(data, group=="Pre")
data_post <- filter(data, group=="Post")
data_dapi <- filter(data, group=="DAPI")

nPre <- nrow(data_pre)                                               #number of pre markers
nPost <- nrow(data_post)                                             #number of post markers
nDAPI <- nrow(data_dapi)

gen_log(paste("Number of Pre marker:", nPre), logfile)
gen_log(paste("Number of Post marker:", nPost), logfile)
gen_log(paste("Number of DAPI marker:", nDAPI), logfile)


gen_log("Plot markers", logfile, datetime=TRUE)
p_marker_coord_facet <- ggplot(data, aes(x = X, y = Y, colour = group)) +
  #geom_point(size=0.5)+
  geom_point(size=0.5, aes(alpha=1/255*Mean))+
  facet_grid(.~group)+
  {if (image_3D) labs(title = "Marker (z-project)") 
  else labs(title = "Marker")} +
  xlab("x(um)") +
  ylab("y(um)") +
  scale_y_reverse()+
  scale_color_manual(values=c("chartreuse3", "red", "royalblue3"))+
  basic_style + theme(legend.position="none")

p_marker_coord <- ggplot(data, aes(x = X, y = Y, colour = group)) +
  geom_point(size=0.5, aes(alpha=1/10))+
  #geom_point(size=2, aes(alpha=1/255*Mean))+
  #facet_grid(.~group, margins = TRUE)+
  #scale_color_manual(values=c("chartreuse3", "red", "royalblue3"))+
  {if (image_3D) labs(title = "Merged Markers (z-project)") 
  else labs(title = "Merged Markers")} +
  xlab("x(um)") +
  ylab("y(um)") +
  scale_y_reverse()+
  scale_color_manual(values=c("chartreuse3", "red", "royalblue3"))+
  basic_style + theme(legend.position="none")

pdf("marker_coord.pdf", width = 10, height = 3)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,3)))
  print(p_marker_coord_facet, vp = vplayout(1,1:2))
  print(p_marker_coord, vp = vplayout(1,3))
dev.off()
gen_log("Done", logfile, datetime=TRUE)

gen_log("Make 2D histogram of coordinates", logfile)
#make hitogram of X  coord distribution
h_X <- ggplot(data, aes(x=X, y=..density.., fill=group)) +
  geom_histogram(bins=60) +
  geom_density(alpha=0.3, size=0.3, adjust=1/4)+
  facet_grid(group ~ ., scales = "free")+
  scale_fill_manual(values=c("chartreuse3", "red", "royalblue3"))+
  labs(title = "X coord distribution") +
  xlab("coordinate (um)") +
  basic_style+
  theme(legend.position="none")
#ggsave("X_coord_histo.pdf", width = 16, height = 12, units = "cm")

#make histogram of Y coordinate distribution
h_Y <- ggplot(data, aes(x=Y, y=..density.., fill=group)) +
  geom_histogram(bins = 60) +
  geom_density(alpha=0.3, size=0.3, adjust=1/4)+
  facet_grid(group ~ ., scales = "free") +
  scale_fill_manual(values=c("chartreuse3", "red", "royalblue3"))+
  labs(title = "Y coord distribution") +
  xlab("coordinate (um)") +
  basic_style +
  theme(legend.position="none")
#ggsave("Y_coord_histo.pdf", width = 16, height = 12, units = "cm")

#make histogram of Z coordinate distribution
if (image_3D) {
  h_Z <- ggplot(data, aes(x=Z, y=..density.., fill=group)) +
    geom_histogram(bins = 60) +
    geom_density(alpha=0.3, size=0.3)+
    facet_grid(group ~ ., scales = "free") +
    scale_fill_manual(values=c("chartreuse3", "red", "royalblue3"))+
    labs(title = "Z coord distribution") +
    xlab("coordinate (um)") +
    basic_style +
    theme(legend.position="none")
#ggsave("Z_coord_histo.pdf", width = 16, height = 12, units = "cm")
}

if (image_3D) {
  pdf("XYZ distributio hitogram.pdf", width = 15, height = 10)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1,3)))
    print(h_X, vp = vplayout(1,1))
    print(h_Y, vp = vplayout(1,2))
    print(h_Z, vp = vplayout(1,3))
  dev.off()
}else {
  pdf("XY distributio hitogram.pdf", width = 15, height = 10)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1,2)))
    print(h_X, vp = vplayout(1,1))
    print(h_Y, vp = vplayout(1,2))
  dev.off()
}

#2d histogram of XY cordinate distribution
h2d_XY <- ggplot(data, aes(x=X, y=Y)) + 
  #stat_bin2d(bins=40)+
  stat_binhex(bins=50)+
  #scale_fill_gradientn(colours=r, trans="log")+
  facet_grid(group ~ .)+
  xlab("X coordinate (um)") +
  ylab("Y coordinate (um)") +
  scale_y_reverse()+
  labs(title = "XY coord distribution") 
#ggsave("XY_coord_histo.pdf", width = 16, height = 12, units = "cm")

if (image_3D) {
  h2d_XY <- h2d_XY + theme(plot.title = element_text(lineheight=.8, face="bold"), legend.position="bottom")
  #2d histogram of XZ cordinate distribution
  h2d_XZ <- ggplot(data, aes(x=X, y=Z)) + 
    #stat_bin2d(bins=20)+
    stat_binhex(bins=50)+
    #scale_fill_gradientn(colours=r, trans="log")+
    facet_grid(group ~ .)+
    xlab("X coordinate (um)") +
    ylab("Z coordinate (um)") +
    labs(title = "XZ coord distribution") +
    theme(plot.title = element_text(lineheight=.8, face="bold"), legend.position="bottom")
  #ggsave("XZ_coord_histo.pdf", width = 16, height = 12, units = "cm")

  #2d histogram of YZ cordinate distribution
  h2d_ZY <- ggplot(data, aes(x=Z, y=Y)) + 
    #stat_bin2d(bins=20)+
    stat_binhex(bins=50)+
    #scale_fill_gradientn(colours=r, trans="log")+
    facet_grid(group ~ .)+
    xlab("Z coordinate (um)") +
    ylab("Y coordinate (um)") +
    scale_y_reverse()+
    labs(title = "ZY coord distribution") +
    theme(plot.title = element_text(lineheight=.8, face="bold"), legend.position="bottom")
  #ggsave("ZY_coord_histo.pdf", width = 16, height = 12, units = "cm")

  pdf("2D distributio hitogram.pdf", width = 15, height = 10)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1,3)))
    print(h2d_XY, vp = vplayout(1,1))
    print(h2d_XZ, vp = vplayout(1,2))
    print(h2d_ZY, vp = vplayout(1,3))
  dev.off()
} else {
  ggsave("2D distributio hitogram.pdf", width = 8, height = 10)
}

#calculate distance between pre / post points

#store marker XYZ values in coord (matrix)
coord <- matrix(0, nrow=6, ncol=nPost * nPre)
row.names(coord) <- c("PostX", "PostY", "PostZ", "PreX", "PreY", "PreZ")
temp <- rep(data_pre$X, nPost)                                       #X coord of pre marker
dim(temp) <- c(nPre, nPost)
coord[1,] <- c(t(temp))
temp <- rep(data_pre$Y, nPost)                                       #Y coord of pre marker
dim(temp) <- c(nPre, nPost)
coord[2,] <- c(t(temp))
if (image_3D){
  temp <- rep(data_pre$Z, nPost)
} else {
  temp <- rep(NA, nPre*nPost)
}                                          #Z coord of pre marker
dim(temp) <- c(nPre, nPost)
coord[3,] <- c(t(temp))
coord[4,] <- rep(data_post$X, nPre)                                  #X coord of post marker
coord[5,] <- rep(data_post$Y, nPre)                                  #Y coord of post marker
if (image_3D){
  coord[6,] <- rep(data_post$Z, nPre)
} else {
  coord[6,] <- rep(NA, nPre)
}                                      #Z coord of pre marker

#exclud marker that are apart more than 5
coord <- coord[,((coord["PostX",]-coord["PreX",])<5)&((coord["PostY",]-coord["PreY",])<5)]
row.names(coord) <- c("PostX", "PostY", "PostZ", "PreX", "PreY", "PreZ")


gen_log("Calculating distance of experimental data...", logfile, datetime=TRUE)
dist2 <- apply(coord, 2, calcDist)
dist <- sqrt(dist2)
coord <- rbind(coord[1:6, ], dist)
row.names(coord)[7]="distance"
distances <- tibble(distance=double(), group=character())
if (length(coord[7,])) 
  distances <- bind_rows(distances, data.frame(distance=as.vector(coord[7,]), group="experimental"))

#subset by distance parameters
coord_subset1 <- as.matrix(coord[, dist_cutoff_min < coord[7,] & coord[7,] < dist_cutoff_max])

#calculate  Y angles between Pre/Post Pairs
gen_log("Calculating angles: experimental/Y", logfile, datetime=TRUE)
angle.y <- apply(coord_subset1, 2, calcAngle.y)
coord_subset1 <- rbind(coord_subset1, angle.y)
row.names(coord_subset1)[8]="angleY"
angles <- tibble(angle=double(), group=character(), axis=character())
if (length(coord_subset1[8,])) 
  angles <- bind_rows(angles, data.frame(angle=as.vector(coord_subset1[8,]), group="experimental", axis="y"))


#calculate  Z angles between Pre/Post Pairs
gen_log("Calculating angles: experimental/Z", logfile, datetime=TRUE)
angle.z <- apply(coord_subset1, 2, calcAngle.z)
coord_subset1 <- rbind(coord_subset1, angle.z)
row.names(coord_subset1)[9]="angleZ"
if (length(coord_subset1[9,])) 
  angles <- bind_rows(angles, data.frame(angle=as.vector(coord_subset1[9,]), group="experimental", axis="z"))


#make histogram of distances
h_dist <- ggplot(distances, aes(x=distance)) +
  labs(title = "Pre/Postsynaptic marker distance histogram") +
  xlab("distance (um)") +
  basic_style    
if (nrow(distances)){
  h_dist <- h_dist +
  geom_histogram(stat = "bin", binwidth = 0.05, alpha=0.4, position="identity", aes(y = ..density..)) +
  geom_rect(data=data.frame(xmin=dist_cutoff_min, xmax=dist_cutoff_max, ymin=-Inf, ymax=Inf), 
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="deepskyblue", alpha=0.1, inherit.aes = FALSE) +
  xlim(0,5) 
  # geom_density(alpha=0.4)+
  # stat_density(adjust = 0.01)+
}

#ggsave("dist_histo.pdf", width = 16, height = 12, units = "cm")

#make histogram of angles
gen_log("Plot histogram of angles", logfile)
h_angle <- ggplot(angles, aes(x = angle)) + 
  labs(title = "Pre/Postsynaptic marker angle histogram") +
  xlab("angle (rad)") +
  basic_style
if (nrow(angles)){
  h_angle <- h_angle +
    geom_histogram(stat = "bin", binwidth=0.2, alpha=0.5, position="identity", aes(y = ..density..)) +
    {if (image_3D) facet_grid(. ~ axis)}
  h_angle <- fracAx(h_angle, "pi")
}
#ggsave("Y angle_histo.pdf", width = 16, height = 12, units = "cm")

pdf("distance & angle hitogram.pdf", width = 15, height = 10)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,2)))
  print(h_dist, vp = vplayout(1,1))
  print(h_angle, vp = vplayout(1,2))                                              # Make the graph with pi axis
dev.off()

gen_log("Mapping coordinates within dist & angle parameters", logfile)
#markers within parameters
map_coord <- function(d, m){                                                      #d is data.frame, m is matrix
  l <- list()
  for (i in 1:dim(m)[2]){
    pre <- d %>% filter(group == "Pre" & d$X %in% m[1,i] & d$Y %in% m[2,i]) %>%
                 bind_cols(., tibble(pair_dist = m[7,i], pair_angle = m[8,i], pair_id = i ))              #add distance and id
    post <- d %>% filter(group == "Post" & d$X %in% m[4,i] & d$Y %in% m[5,i]) %>%
                   bind_cols(., tibble(pair_dist = m[7,i], pair_angle = m[8,i], pair_id = i ))              #add distance and id    
    l[[i]] <- bind_rows(pre, post)
  }
  bind_rows(l)
}

#calculate marker density of a square of size distance^2 
calc_density <- function(x, y, df, distance){    #df is tibble with X & Y coordinates
  df %>% filter((X > x - distance/2 & X < x + distance/2) & (Y > y - distance/2 & Y < y + distance/2)) %>%
    nrow() -> count
  (count)/(distance^2)
}

#if there are no pairs return empty tibble to avoid errors
data_pair <- tibble()
if (length(coord_subset1)){
  data_pair <- map_coord(data, coord_subset1)

  #qunatify shape of synapse (ribeye horseshoe shape)
  data_pair <- data_pair %>% 
                  mutate(shape = sqrt((XM-(BX+Width/2))^2+(YM-(BY+Height/2))^2))

  gen_log("Calculating local marker densities", logfile, datetime=TRUE)
  gen_log(paste("Distance for density calculation", distance_density, "um"), logfile)
  data_pair$density <- vector("double", nrow(data_pair))  
  for (r in 1:nrow(data_pair)){
    data_pair$density[r] <- calc_density(data_pair$X[r], data_pair$Y[r], 
                                         filter(data, group==as.character(data_pair$group[r])), distance_density)
  }
}

# bw <- data_pair$density %>% min
# ggplot(data_pair, aes(x=density), color=group)+
#   geom_histogram(binwidth=bw)+
#   facet_grid(.~group)


#export data_pair data_bg once 
if (file.exists("data_pair.rds")){
  gen_log("Deleted old data_pair.rds", logfile)
  unlink("data_pair.rds", recursive = TRUE)
}

saveRDS(data_pair, file = "data_pair.rds")  #subset of data_bg with pairs within synapse distance
gen_log("Pair tibble created: data_pair.rds", logfile)

gen_log("Plotting marker stats", logfile, datetime=TRUE)

plot_feature <- function(feature, ...){
  h <- ggplot(data_pair) + xlab(feature) + basic_style
  if (nrow(data_pair)){
    h <- h + 
      geom_histogram(aes(x=eval(parse(text=feature)), fill=group), position="identity", alpha=1/3, ...) 
      
  }
  assign(paste0("h_", feature), h, envir=globalenv())
}

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

if (nrow(data_pair)){
  features <- colnames(data_pair[sapply(data_pair,is.double)])                         #all the stats     
  #remove some unecessary stats
  features <- features[!features %in% c("X", "Y", "XM","YM", "BX", "BY", "bg", 
                                        "X.Area", "MinThr", "MaxThr", "PercentArea", 
                                        "FeretX", "FeretY")]
  print(features)
  
  for (feature in features){
    bw <- robust_hpi(data_pair %>% dplyr::select(feature))
    plot_feature(feature, binwidth=bw)
  }

  geometry <- list(h_pair_dist, h_pair_angle, h_Area, 
                     h_IntDen, h_RawIntDen)
  h_geometry <- grid_arrange_shared_legend(geometry, ncol = 3, nrow = 2)
  ggsave("Feature distribution (geometry).pdf", plot=h_geometry , width = 12, height = 6)
  
  signal <- list(h_Mean, h_Median, h_Mode, h_StdDev, 
                 h_Min, h_Max, h_Skew, h_Kurt)
  h_signal <- grid_arrange_shared_legend(signal, ncol = 4, nrow = 2)
  ggsave("Feature distribution (signal).pdf", plot=h_signal , width = 12, height = 12)
  
  morphology <- list(h_Perimeter, h_Width, h_Height, h_shape, 
                     h_Feret, h_MinFeret, h_FeretAngle, h_density,
                     h_Major, h_Minor, h_Angle, h_AR,
                     h_Round, h_Circularity, h_Solidity)
  h_morphology <- grid_arrange_shared_legend(morphology, ncol = 4, nrow = 4)
  ggsave("Feature distribution (morphology).pdf", plot=h_morphology , width = 12, height = 12)
}

gen_log("find pair done", logfile, datetime = TRUE)
