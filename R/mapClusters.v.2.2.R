##' Map clusters
##' 
##' This function will map the clusters created by percolation. 
##' 
##' shape is a SpatialPolygonDataFrame (sp), which will be used for plotting as a background map, and determining the coordinate system.
##' 
##' map_name is a string, which will be used to title the maps.
##' 
##' source_file_name is a string, which will be used in a subtitle within the maps and should describe the source data used.
##' 
##' dpi is the output file resolution. The default is set to 300. Other resolutions may be needed for publication, for example.
##' 
##' Exported maps are in the folder "/maps": As many png map files as there are different radii values used in the percolation analysis. Also a shp-file with the point data.
##' 
##' For more information and a code and data example please check the vignette.
##'
#' @author Simon Maddison
#' @author Sophie C. Schmidt
#' @import sp
#' @import stats
#' @import igraph
#' @import maptools
#' @import stringr
#' @import grDevices
#' @import graphics
#' @import rgdal
#' @import raster
#' @import utils
#' @param shape needs input of SpatialPolygonDataFrame from 'sp' for background map (same ESPG!)
#' @param map_name needs string for how the maps shall be named
#' @param source_file_name string for the name of data set used, to include on maps
#' @param dpi dots per inch for the images. set to 300 as standard
#' @return as many maps as png as step_values have been calculated
#' 
#' @export mapClusters
#'
mapClusters <- function(shape, map_name, source_file_name, dpi = 300) {

  # set up path results and path_working to better seperate working data and analysis results data  

  # file.path creates paths according to platform used on user's computer
path_working <- file.path(getwd(), "working_data")
path_results <- file.path(getwd(), "analysis_results")

# create maps directory
path_maps <- file.path(getwd(), "maps")
dir.create(path_maps, showWarnings = FALSE)

# read in nodes list and coordinates
data_file <- paste(file.path(path_working,"nodes_list_d.txt"))

# Read in distance thresholds - this ensures same values used as in clustering script, where it was saved
file_name <- paste(file.path(path_working,"working_data.csv"))
radius_values <- read.csv(file_name,header=TRUE, sep = ",")
upper_radius <- radius_values$upper_radius
lower_radius <- radius_values$lower_radius
step_value <- radius_values$step_value
radius_unit <- radius_values$radius_unit

map_outline <- shape

# Read in nodes and grid coordinates
data_filename <- paste(file.path(path_working,"data.csv"))
xy_data <- read.csv2(data_filename, sep = ",")

xy_data$Easting <- as.numeric(as.character(xy_data$Easting))
xy_data$Northing <- as.numeric(as.character(xy_data$Northing))

# Convert data to SpatialPointsDataFrame, using project4string from selected map file
# get projection string for current map
projection_string <- proj4string(map_outline)
crs_projection <- CRS(projection_string)
# turn data into SpatialPointsDataFrame
coordinates(xy_data) <- c("Easting", "Northing")
proj4string(xy_data) <- crs_projection

# Read in data file with cluster id for each site and radius
file_name <- paste(file.path(path_results,"member_cluster_by_radius.csv"))
ranked_mem_clust_by_r <- read.csv(file_name, header=TRUE)


#check on decimal places if step_value non integer
dec_places <- if ((step_value %% 1) != 0)
	{nchar(strsplit(sub('0+$', '', as.character(step_value)), ".", fixed=TRUE)[[1]][[2]])
	} else {
	0}


if (radius_unit == 1)
{unit_text <- "m"
} else if (radius_unit == 1000)
{unit_text <- "km"
} else {
unit_text <- paste(radius_unit, "m", sep="")}

radius_values <- seq(upper_radius,lower_radius,by=-step_value)
# this accomodates non-integer radius values
radii_count <- length(radius_values)
loop_count <- seq(radii_count,1,by=-1)

# Define colours for clusters; top 15 and then gray for the rest
top15_colours <- colors()[c(553,29,258,654,91,115,456,122,48,8,149,86,102,40,12)]
the_rest_colour <- "#AEAEAE"

# this code can be used to plot all outputs as pdf as an alternative
# this code would need to be substiuted for the png plot code below, and in the for loop
# Maps plotted on multipage pdf
#file_map_pdf <- paste(path_maps,"/","percolation_plots_",map_name,".pdf",sep="")
#pdf(file=file_map_pdf, paper="a4", width=21/2.54, height=29.7/2.54)

# Plot as png, default. png is deemed a more practical solution
# page format default is A4
file_map_png <- paste(file.path(path_maps, "percolation_plots_"),map_name, "_all.png",sep="")
png(file=file_map_png, units="cm", width=21, height=29.7, res=dpi)
# plot the map outline
plot(map_outline, col="white",border=TRUE)

# point size adjusted for plots of more or less than 1000 points
if((nrow(xy_data) < 1000))
{point_dia <- 0.8
} else {
point_dia <- 0.4}

points(xy_data$Easting, xy_data$Northing, col='red', pch=20, cex=point_dia)
number_of_sites <- paste("Number of sites: ",nrow(xy_data))
mtext(number_of_sites)
plot_title <- str_to_title(map_name,locale="")
title(paste("Clusters of ", plot_title), sub=paste("Sources: ",source_file_name,"; ",map_name))
dev.off()

# Generates maps for each of the percolation radii in range of radius values
for(i in loop_count)
{
	radius <- radius_values[i]
	radius_name <- format(radius,nsmall=dec_places)
	file_map_png <- paste(file.path(path_maps, "percolation_plots_"),map_name, "_rad_", radius_name, ".png",sep="")
	png(file=file_map_png, units="cm", width=21, height=29.7, res=dpi)
	# Uses data for percolation radius computed, steps through each column for radius value defined by i
	ClstRad_col <- paste("ClstRad",radius,sep="")
	# extract data for this radius from cluster member by radius
	member_cluster <- ranked_mem_clust_by_r[c("PlcIndex", ClstRad_col)]
	colnames(member_cluster) <- c("node","cluster")

	# Ranks clusters by size, i.e. number of nodes for each cluster
	#  creates the number of instances for each cluster id ...
	frequency_clusters <- data.frame(table(member_cluster$cluster))
	names(frequency_clusters) <- c('cluster','number_of_points')
	ranked_clusters <- frequency_clusters[order(frequency_clusters$number_of_points, decreasing=TRUE),]

	# Number of nodes total for this radius
	total_nodes <- sum(ranked_clusters$number_of_points)
	# The number of clusters for this radius
	total_clusters <- nrow(ranked_clusters)
	print(paste("For radius: ",radius," there are: ",total_clusters," clusters and ",total_nodes," nodes"))

	# Pick top 15 (if there are as many as this)
	number_colours <- 15
	if(total_clusters < 15)
	{
		number_colours <- total_clusters
	}
	# Giving colours to the top 15 or the top ones when less than 15
	ranked_clusters$col = ''
	ranked_clusters$col[1:number_colours] <- top15_colours[1:number_colours]
	if(total_clusters>15)
	{
		# add colours of grey for the rest
		ranked_clusters$col[16:total_clusters] <- the_rest_colour
	}

    	#get sub-data from original xy_data
    	xy_at_d <- xy_data[xy_data$PlcIndex %in% member_cluster$node,]

    	#first add colour column to member_cluster then to xy_at_d
   	  member_cluster$col <- ranked_clusters$col[match(member_cluster$cluster, ranked_clusters$cluster)]
    	xy_at_d$col <- member_cluster$col[match(xy_at_d$PlcIndex, member_cluster$node)]
    	
    	# get extent of map for legend plotting
    	the_plot_extent <- extent(map_outline)
    	
    	# grab the upper right hand corner coordinates
    	furthest_pt_east <- the_plot_extent@xmax
    	furthest_pt_north <- the_plot_extent@ymax

  # margins to make space for plot, heading and sub-heading
  par(mar = c(5, 4, 4, 7))  	
	plot(map_outline, col="white",border=TRUE)

	# plot all points in grey, as background
	points(xy_data$Easting, xy_data$Northing, col='grey85', pch=4, cex=.3)
	# plot
	points(xy_at_d$Easting, xy_at_d$Northing, col=xy_at_d$col, pch=20, cex=point_dia)
	
	# legend
	legend_size <- 0.8
	legend(x = furthest_pt_east+20, y = furthest_pt_north+20, title="Rank Cluster#",legend=ranked_clusters$cluster[1:number_colours], cex=legend_size,
		fill=ranked_clusters$col[1:number_colours],  xpd = TRUE)
	
	title(main=plot_title, sub=paste("Source File: ",source_file_name, "; percolation distance:",radius_name," ",unit_text))
	
	# Create data for output in csv file
	# create ranking index for cluster numbers
	ranked_clusters$rank <- 1:nrow(ranked_clusters)
	ranked_clusters$number_of_points <- NULL
	ranked_clusters$col <- NULL
	new_col_name <- paste("RankRad", radius, sep="")
	colnames(ranked_clusters) <- c("cluster",new_col_name)
	# add cluster ranking as column into data frame with clusters and PlcIndex
	ranked_mem_clust_by_r <- merge(ranked_mem_clust_by_r, ranked_clusters, by.x=ClstRad_col, by.y="cluster",
		 all.x=TRUE)
	# Remove original cluster column for this radius
	ranked_mem_clust_by_r[,ClstRad_col] <- NULL
	dev.off()
}

# convert NA's in ranked_mem_clust_by_r to zero
ranked_mem_clust_by_r[is.na(ranked_mem_clust_by_r)] <- 0
# Write file with PlcIndex and cluster rank for each radius
# Output file name and location writes to a text file
file_name <- paste(file.path(path_results,"site_cluster_rank_by_radius.csv"))
# need to write this WITHOUT the row number
write.csv(ranked_mem_clust_by_r[order(ranked_mem_clust_by_r$PlcIndex),], file_name, row.names=FALSE)


# Write shape file with all data points
# First convert data to SpatialPointsDataFrame, using project4string from selected map file
# get projection string for current map
projection_string <- proj4string(map_outline)
crs_projection <- CRS(projection_string)

# name of output shapefile
output_shape_file <- paste(file.path(path_maps,source_file_name),".shp", sep = "")
layer_name <- "percolation"

xy_data <- sp::merge(xy_data, ranked_mem_clust_by_r, by = "PlcIndex")

writeOGR(xy_data,output_shape_file,layer=layer_name, driver="ESRI Shapefile",overwrite_layer=TRUE)

# output source file name for plotClustFreq
source_file_name_out <- paste(file.path(path_working,"source_file_name_out.csv"))
write.table(source_file_name, source_file_name_out)

}
