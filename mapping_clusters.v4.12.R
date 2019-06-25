# CLUSTER MAPPING SCRIPT
# Original code from Elsa Arcaute, CASA, UCL 28.2.15; extensively modified by Simon Maddison.
# Cite as:	arXiv:1504.08318 [physics.soc-ph]
# Version 3a: 10/5/15 & 2/10/15 - with additions to include county boundary outlines,
# This code maps the clusters overlaid on an outline of the British Isles or Ireland including county outlines
# England and Wales shape files added from Portsmouth University files 31/1/16
# Version 4: 24/3/16. Changes to plots made to:
# 	Include identity of clusters (i.e.cluster index) in legend - done
# 	Plot the source of the original data (e.g. Atlas, Hogg etc.) (taken from source file name)- done
# 	Parameterise the shape file to be used for plotting ...
#	Take out graded grey colours for clusters > 15
# Version 4.1: 29/3/16 changed outputs to be multi-page pdf
#	NOTE: Ireland not plotting points TO BE RESOLVED
# Version 4.2: 6/4/16 - error in plot of all site corrected; removed 'xy'; all refs now to xy_data
# Version 4.3: 10/4/16 - output shape file with all data points, with same projection as input map
#     NOTE: esri shape file driver truncates names to 10 characters but this should not be a
#     problem at the moment. Prob of e.g. PlacesIndex, ClusterNumber etc.
# 	shape file name is name of source file without .csv etc. extension
# Version 4.4: 29/4/16 - paging upgraded to export to A4 and use the full page
# Version 4.5: 6/5/16 - changed Index name to fit within 8 chars -> PlcIndex
# Version 4.6: 8/5/16 - working on Ireland plot to fix. Problem seemed to be an old Ireland map with
#	an obsolete coordinate system
#	Add in external input to country and shape file for source,
#	file is added to source_file.txt - map name and shape file name, headers also added
#	source_file, map_name, shape_file
#	REQUIRES CHANGE TO ALL OTHER PROGRAMS (done)
# Version 4.7: 18/5/16 -  Cluster data now comes in a single file (member_cluster_by_radius.csv)
#   	output cluster ranks as new file site_cluster_rank_by_radius, lists by PlcIndex the cluster rank
#	indexed by radius
# Version 4.8: 30/7/16 - Plot as png file output. Earlier issues with resolution resolved.
#	Easier to include in word document than pdf. However each plot is a separate file,
#	with suitable name increment
# Version 4.9: 03/08/2016 - minor change to put legend on right for Domesday counties
# Version 4.10: 7/11/2017 - working directory amended for new pc
#	Tweaks to legend location, size and plot title to allow for other data sets
# Version 4.11: 24/06/2018 - only worked for integer increments of radius. Modified to allow fractional
#	values of radius, e.g. step size of 0.1km
# Version 4.12: 04/12/18: - add in unit parameter to select km or m, in "radius_values.txt"
#	NOTE: values are assumed to be in metres, unit will factor this in calculations and charts by the unit value
#	so that a unit value of 1 shows through as metres, and a value of 1000 as km. Other values accomodated
#	Also modified so that when incrementing in decimal step values, all map names are formatted consistently

library(maptools)
library(rgdal)
library(stringr)

#setwd("/home/sophie/Dokumente/Konferenzen/percolatransect/") #path on laptop
# paths for results and output maps
path_source <- paste(getwd(),"/source_data",sep="")
path_results <- paste(getwd(),"/working_data",sep="")
path_maps <- paste(getwd(),"/maps",sep="")
path_shape <- paste(getwd(),"/shape_files",sep="")

# File is the list of nodes and x y grid references
# Read in source file name defined in file source_file.txt
# Edit this file for different input file name
file_name <- paste(path_results,"/","source_file.txt",sep="")
source_files <- read.csv(file_name,header=TRUE,stringsAsFactors=FALSE)
source_file_name <- source_files$source_file
source_file <- paste(path_source,"/",source_file_name,sep="")
map_name <- source_files$map_name
shape_file_name <- source_files$shape_file

# Load shape file
file_name <- paste(path_shape,"/",shape_file_name,sep="")
# Truncate shape file name to remove extension
layer_name <- substr(shape_file_name,1,(nchar(shape_file_name)-4))
print(map_name)
map_outline <- readOGR(dsn=file_name,layer=layer_name)

# Read in nodes and grid coordinates
xy_data <- read.csv(source_file, as.is = FALSE, sep = "\t") #added tab-seperation
# as.is is redundant as far as I can see as just confirms default behaviour of read.csv
# fields of interest:
# PlcIndex Easting Northing

# Convert data to SpatialPointsDataFrame, using project4string from selected map file
# get projection string for current map
projection_string <- proj4string(map_outline)
crs_projection <- CRS(projection_string)
# turn data into SpatialPointsDataFrame
coordinates(xy_data) <- c("Easting", "Northing")
proj4string(xy_data) <- crs_projection

# Read in data file with cluster id for each site and radius
file_name <- paste(path_results,"/","member_cluster_by_radius.csv",sep="")
ranked_mem_clust_by_r <- read.csv(file_name, header=TRUE)

# Read in distance thresholds - this ensures same values used as in clustering script
file_name <- paste(path_results,"/","radius_values.txt",sep="")
radius_values <- read.csv(file_name,header=TRUE)
upper_radius <- radius_values$upper_radius
lower_radius <- radius_values$lower_radius
step_value <- radius_values$step_value

#check on decimal places if step_value non integer
dec_places <- if ((step_value %% 1) != 0)
	{nchar(strsplit(sub('0+$', '', as.character(step_value)), ".", fixed=TRUE)[[1]][[2]])
	} else {
	0}

radius_unit <- radius_values$radius_unit

if (radius_unit == 1)
{unit_text <- "m"
} else if (radius_unit == 1000)
{unit_text <- "km"
} else {
unit_text <- paste(radius_unit, "m", sep="")}

radius_values <- seq(upper_radius,lower_radius,by=-step_value)
# Changes to accomodate non-integer radius values
radii_count <- length(radius_values)
loop_count <- seq(radii_count,1,by=-1)

# Define colours for clusters, top 15 and then gray for the rest
top15_colours <- colors()[c(553,29,258,654,91,115,456,122,48,8,149,86,102,40,12)]
the_rest_colour <- "#AEAEAE"

# plot all sites as pdf
# Maps plotted on multipage pdf
#file_map_pdf <- paste(path_maps,"/","percolation_plots_",map_name,".pdf",sep="")
#pdf(file=file_map_pdf, paper="a4", width=21/2.54, height=29.7/2.54)

# Plot as png, earlier issues with png now resolved

file_map_png <- paste(path_maps,"/","percolation_plots_",map_name, "_all.png",sep="")
png(file=file_map_png, units="cm", width=21, height=29.7, res=300)

plot(map_outline, col="white",border=TRUE)

# get from original file xy_data coords and create table
# Easting | Northing | ID | cluster | colour
if((nrow(xy_data) < 1000)| (map_name == 'Domesday'))
{point_dia <- 0.8
} else {
point_dia <- 0.4}

points(xy_data$Easting, xy_data$Northing, col='red', pch=20, cex=point_dia)
number_of_sites <- paste("Number of sites: ",nrow(xy_data))
mtext(number_of_sites)
plot_title <- str_to_title(map_name,locale="")
title(paste("Iron Age Hillforts in ", plot_title), sub=paste("Sources: ",source_file_name,"; ",shape_file_name))
dev.off()

# Generates maps for each of percolation radii in range of radius values
for(i in loop_count)
{
	radius <- radius_values[i]
	radius_name <- format(radius,nsmall=dec_places)
	file_map_png <- paste(path_maps,"/","percolation_plots_",map_name, "_rad_", radius_name, ".png",sep="")
	png(file=file_map_png, units="cm", width=21, height=29.7, res=300)
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

	plot(map_outline, col="white",border=TRUE)
	# get from original file xy_data coords and create table
	# Easting | Northing | ID | cluster | colour
	# this following to check all points included
	# plot all points in grey, as background
	points(xy_data$Easting, xy_data$Northing, col='grey85', pch=4, cex=.3)
	# plot
	points(xy_at_d$Easting, xy_at_d$Northing, col=xy_at_d$col, pch=20, cex=point_dia)
	if(plot_title=="Ireland" | plot_title=="Domesday" | plot_title=="Domesday Vills")
	{	legend_loc <- "bottomright"
		legend_size <- 0.6
	} else {
		legend_loc <- "bottomleft"
		legend_size <- 0.8}
	legend(legend_loc, inset= .01, title="Rank Cluster#",legend=ranked_clusters$cluster[1:number_colours], cex=legend_size,
		fill=ranked_clusters$col[1:number_colours])

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
dev.off()

# convert NA's in ranked_mem_clust_by_r to zero
ranked_mem_clust_by_r[is.na(ranked_mem_clust_by_r)] <- 0
# Write file with PlcIndex and cluster rank for each radius
# Output file name and location writes to a text file
file_name <- paste(path_results,"/","site_cluster_rank_by_radius.csv",sep="")
# need to write this WITHOUT the row number
write.csv(ranked_mem_clust_by_r[order(ranked_mem_clust_by_r$PlcIndex),], file_name, row.names=FALSE)

# Write shape file with all data points
# First convert data to SpatialPointsDataFrame, using project4string from selected map file
# get projection string for current map
projection_string <- proj4string(map_outline)
crs_projection <- CRS(projection_string)

# name of output shapefile
# Truncate source file name to remove extension
source_file_name <- substr(source_file_name,1,(nchar(source_file_name)-4))
output_shape_file <- paste(path_maps,"/",source_file_name,".shp",sep="")
writeOGR(xy_data,output_shape_file,layer=source_file_name,driver="ESRI Shapefile",overwrite_layer=TRUE)

