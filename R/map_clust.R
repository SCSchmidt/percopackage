
mapClusters <- function(data, shape, map_name) {

  # changed path results and path_working to better seperate working data and analysis results data  
path_results <- paste(getwd(),"/analysis_results",sep="")
path_working <- paste(getwd(),"/working_data",sep="")
path_maps <- paste(getwd(),"/maps",sep="")
data_file <- paste(path_working,"/","nodes_list_d.txt",sep="")

# Read in distance thresholds - this ensures same values used as in clustering script, where it was saved
file_name <- paste(path_working,"/","working_data.csv",sep="")
radius_values <- read.csv(file_name,header=TRUE, sep = ",")
upper_radius <- radius_values$upper_radius
lower_radius <- radius_values$lower_radius
step_value <- radius_values$step_value
radius_unit <- radius_values$radius_unit

map_outline <- shape
# Truncate shape file name to remove extension
#layer_name <- substr(shape,1,(nchar(shape)-4))
#print(map_name)

# Read in nodes and grid coordinates
xy_data <- data

# Convert data to SpatialPointsDataFrame, using project4string from selected map file
# get projection string for current map
projection_string <- proj4string(map_outline)
crs_projection <- CRS(projection_string)
# turn data into SpatialPointsDataFrame
coordinates(xy_data) <- c("Easting", "Northing")
proj4string(xy_data) <- crs_projection

# Read in data file with cluster id for each site and radius
file_name <- paste(path_working,"/","member_cluster_by_radius.csv",sep="")
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
title(paste("Clusters of ", plot_title), sub=paste("Sources: ",data,"; ",map_name))
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
	
	title(main=plot_title, sub=paste("; percolation distance:",radius_name," ",unit_text))

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

# name of output shapefile (renamed by Sophie)
# Truncate source file name to remove extension
output_shape_file <- paste(path_maps,"/analysis_results.shp",sep="")
layer_name <- "percolation"

xy_data <- sp::merge(xy_data, ranked_mem_clust_by_r, by = "PlcIndex")

writeOGR(xy_data,output_shape_file,layer=layer_name, driver="ESRI Shapefile",overwrite_layer=TRUE)

}
