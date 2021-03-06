##' Percolation
##' 
##' This is a function built based on the scripts by Simon Maddison, who based it on code developed by Elsa Arcaute, which is published under:  
##' E.Arcaute, C.Molinero, E.Hatna, R.Murcio, C.Vargas-Ruiz, A.P.Masucci and M.Batty. 
##'  Cities and Regions in Britain through hierarchical percolation. J. R. Soc. Open Science, 3 (4),150691. 
##'  doi:10.1098/rsos.150691, 2016. 
##' 
##' 
##' data = A data.frame, with the format and example data below:
##' 
##' PlcIndex,Easting,Northing
##' 
##' 1,350350,233050
##' 
##' 2,354700,260200
##' 
##' 3,358700,238900
##' 
##' PlcIndex is an assigned index number to distinguish the points. This does not need to be sequential and can for example be drawn from a source database.
##' 
##' Easting and Northing are self-explanatory, and are to the metre. (Note that additional columns, if any, will be ignored).
##' 
##' distance_table is set to be NULL as a default, which means, using the data input (coordinates) a distance matrix is calculated using Euclidean distance. If not NULL, the data frame given here needs to adhere to this format and example values:
##' 
##' ID1, ID2, d12
##' 
##' 1,2,12.1
##' 
##' 1,3,14.2
##' 
##' 2,3,2.9
##' 
##' ID1 and ID2 are point IDs (i.e. the PlcIndex values) and 
##' 
##' d12 denotes the weighted distance between them, if a non-euclidian distance is to be used. 
##' 
##' upper_radius is an integer, the upper value of the radius range to be used.
##' 
##' lower_radius is an integer, the lower value of the radius range to be used.
##' 
##' step_value is numeric or integer and the step value to be used between lower and upper radius.
##' 
##' limit is an integer, the value above which distances will not be calculated between sites.
##' 
##' radius_unit  is either 1 for meter or 1000 for km for all input radii.
##' 
##' Exported are in the folder "working_data": tables (csv) with the input-data, the list of PlcIndex, null_entries and duplicate_entries to be used by following functions.
##' 
##' In the folder "analysis_results" the tables (csv) analysis_by_radius and member_cluster_by_radius are exported. 
##' They will be used by the following functions as well, but may also be useful for other applications, therefore they are considered analysis results.

##' For more information and a code and data example please check the vignette "percolation".

#' @author Simon Maddison
#' @author Sophie C. Schmidt
#' @author Elsa Arcaute
#' 
#' @import stats
#' @import Hmisc
#' @import utils
#' @import igraph

#' @param data needs input of data.frame as specified below
#' @param distance_table to give a distance matrix instead of calculating one. Default is NULL, and Euclidean distance will be calculated from the given point data.
#' @param limit needs input of integer: is the value above which distances between sites will not be stored
#' @param lower_radius  needs input of integer, is the lower value of the radius range to be used
#' @param step_value integer or numeric, is the step value to be used between these two values
#' @return the function generates the data files analysis_by_radius.csv and member_cluster_by_radius.csv
#' @param radius_unit is either 1 for meter or 1000 for km; applied to all input values
#' @param upper_radius needs input of integer, is the upper value of the radius range to be used
#' @export
percolate <- function(data, distance_table = NULL, upper_radius, lower_radius, step_value, limit, radius_unit) {

# create directories
# file.path creates paths according to working directory used on user's computer
    
path_working <- file.path(getwd(), "working_data")
path_results <- file.path(getwd(), "analysis_results")

dir.create(path_working, showWarnings = FALSE)
dir.create(path_results, showWarnings = FALSE)

# this will save the input values in a csv and as w_data in the working environment 
w_data_colnames <- c("limit", "radius_unit", "upper_radius", "lower_radius", "step_value")
w_data_values <- c(limit, radius_unit, upper_radius, lower_radius, step_value)
w_data <- rbind(w_data_colnames,w_data_values)

file_name <- paste(file.path(path_working,"working_data.csv"))
write.table(w_data, file_name, row.names=FALSE, col.names = FALSE, sep = ",")

# NOTE: The coordinate data is still required, even if distance table is provided!
# Generate this in the data.csv file, for the mapping program
file_name <-paste(file.path(path_working,"data.csv"))
write.table(data, file_name, row.names=FALSE, col.names = TRUE, sep = ",")

ptm <- proc.time()
# This displays computation time, for information

# List any entries with null values and write to file
data_NA <- data[rowSums(is.na(data)) > 0,]
print(paste("Entries in source file with one or more null values: "))
data_NA
file_name <- paste(file.path(path_working,"null_entries.txt"))
write.table(data_NA, file_name, row.names=FALSE)
# Remove rows with null values from data
data <- na.omit(data)

# conditional depending on what kind of data type you are entering
# if the internode distance table needs to be created, compute and store it
if (is.null(distance_table)) {
	print("Computing inter-node distance table")
	# No distance table provided, so compute the distance table
  
	# Remove points that are superimposed and keep only the first ID
	# - This removes one of two sites that are very close to each other
	# Determined on basis of x y coordinates
	# Write list to file
	duplicate_xy <- data[duplicated(data[,2:3]),]
	data_unique <- data[!duplicated(data[,2:3]),]
	print(paste("Number of removed superimposed points: ",(nrow(data)-nrow(data_unique))))
	duplicate_xy
	file_name <- paste(file.path(path_working,"duplicate_entries.txt"))
	write.table(duplicate_xy, file_name, row.names=FALSE)

	x_vec <- data_unique$Easting
	y_vec <- data_unique$Northing
	ID <- data_unique$PlcIndex

	# Write file with list of PlcIndex values used
	PlcIndex_list <- matrix(ID,ncol=1)
	file_name <- paste(file.path(path_working,"PlcIndex.csv"))
	write.table(data_unique$PlcIndex, file_name, row.names=FALSE, col.names="PlcIndex")

	# Number of points/nodes in file with no duplicates
	n <- length(ID)
	print(paste('number of points: ',n))

	# Create matrix of internodal distances
	# Columns: NodeId1, NodeId2, distance 1-2
	col_list <- cbind('ID1','ID2','d12')
	ni <- n-1
	nj <- n
	n_rows <- n*(n-1)/2
	nodes_list <- matrix(, nrow = n_rows, ncol = 3)
	colnames(nodes_list) <- col_list

	# compute the internode distance table
	row <- 0
	i <- 1
	for (i in 1:ni)
	{
		j1 <- i+1
		for (j in j1:nj)
		{
			# Compute distance between nodes - using Pythagoras
			#  this is in units of 1m
			d <- sqrt((abs(x_vec[i]-x_vec[j]))^2+(abs(y_vec[i]-y_vec[j]))^2)
			# to give distance in m*unit, rounded to 2 decimal places; this also reduces file size
			d <- d/radius_unit # this factors the values by the unit. 1 gives metres, 1000 gives km
			d <- round(d,2)
			# Record this value only if less than limit of distances to be included in the internode distance table
			if(d < limit)
				{row <- row+1
				nodes_list[row,] <- cbind(ID[i],ID[j],as.numeric(d))
			}
		}
	}
	# give processing time, for information
	t1 <- proc.time() - ptm
	print('loop computed')
	
	# Output file name and location of internode distance table to a text file
	file_name <- paste(file.path(path_working,"nodes_list_d.txt"))
	# Remove the unused rows in the matrix
	nodes_list <- nodes_list[-(row+1:n_rows),]
	m_nodes <- as.matrix(nodes_list)
	# need to write this WITHOUT the row number
	write.table(m_nodes, file_name, row.names=FALSE)

	t2 <- proc.time() - ptm
	print('matrix copied')
	
	}
	else {
	  # Distance table is provided, no need to compute it
	  # However the unique IDs are still needed, for mapping
	  print("distance table given, not computed")
	  data_unique <- data.frame(stack(data[,1:2]))
	  data_unique <- unique(data_unique$values)
	  data_unique <- as.data.frame(data_unique)
	  colnames(data_unique) <- "PlcIndex"
	  
	  # Write file with list of PlcIndex values used
	  file_name <- paste(file.path(path_working,"PlcIndex.csv"))
	  write.table(data_unique$PlcIndex, file_name, row.names=FALSE, col.names="PlcIndex")

	  # Save distance_table as nodes_list_d.text so it corresponds to the work flow
	  file_name <- paste(file.path(path_working,"nodes_list_d.txt"))
	  write.table(distance_table, file_name, row.names=FALSE)
	  
}

## original clustering script

mem_clust_by_r <- as.data.frame(data_unique$PlcIndex)
names(mem_clust_by_r)[1] <- "PlcIndex"

# conditional depending on the input unit value

if (radius_unit == 1)
{unit_text <- "m"
} else if (radius_unit == 1000)
{unit_text <- "km"
} else {
  unit_text <- paste(radius_unit, "m", sep="")}

print(paste("Radii used for cluster analysis, ",unit_text,": upper ",upper_radius,
            "; lower ",lower_radius,"; step ", step_value))

# To compute and display computational time
ptm <- proc.time()

data_file <- paste(file.path(path_working,"nodes_list_d.txt"))

# The data table of nodes and internode distances is a Text file, with headers
matrix_IDs_distance <- read.table(data_file,header=TRUE)
# Columns are: node ID1, node ID2, d12 is distance between them. Note that this is generated or 
#  has to be specified beforehand
# If generated, it is with a limit to the maximum distance to reduce overall matrix size, and hence
# creates a partial matrix

# processing time for information
t <- as.vector(proc.time() - ptm)[3]
print(paste('time to get matrix in mins',t/60))
ptm2 <- proc.time()[3]

# Define range of percolation radius and step value to progressively reduce it
# Radius in defined unit value, as per radius_unit in radius_values.txt
radius_values <- seq(upper_radius,lower_radius,by=-step_value)

# To accomodate non-integer radius values
radii_count <- length(radius_values)
loop_count <- seq(radii_count,1,by=-1)

# Repeat cluster computation for each radius value

for (i in loop_count)
{
  radius <- radius_values[i]
  # Create sub-matrix such that all internode distances d12<=radius
  matrix_radius <- matrix_IDs_distance[matrix_IDs_distance[,3]<=radius,]
  # Create graph (note that characters and numerics are treated differently)
  matrix2 <- matrix_radius[,1:2]
  matrix2[,1] <- as.character(matrix2[,1])
  matrix2[,2] <- as.character(matrix2[,2])
  # This creates a graph from the sub-matrix matrix2
  # Directed means that the edges will only be 'counted' once
  # In order to interpret matrix2 as a matrix added 'as.matrix'
  g <- graph.edgelist(as.matrix(matrix2), directed=TRUE)
  
  # take subcomponents - description of how this works
  # http://stackoverflow.com/questions/20725411/finding-strong-and-weak-clusters-and-their-membership-in-r
  
  # Identifies clusters in the graph; creates list of nodes and associated cluster id
  # - weak refers to the mechanism used to generate the clusters and relates to compuational efficiency
  # Note that this does not include clusters of 1 node, so the counts are not meaningful at the lower limit
  member_of_cluster_id <- clusters(g, mode="weak")$membership
  # V labels vertices of the graph
  m <- cbind(V(g)$name,member_of_cluster_id)
  # add the values of the cluster identity for each node to mem_clust_by_r 
  new_col_name <- paste("ClstRad",radius,sep="")
  colnames(m) <- c("PlcIndex",new_col_name)
  mem_clust_by_r <-  merge(mem_clust_by_r, m, by="PlcIndex", all.x=TRUE)
  
  ptm1 <- proc.time()[3]
  t_print <-  as.vector(ptm1-ptm2)/60
  print(paste("for radius=",radius,"it took to compute",sprintf("%4.2f",t_print),"mins"))
  
}
# write the cluster identity file as a csv
file_name <- paste(file.path(path_results,"member_cluster_by_radius.csv"))
# Write file without row names
write.csv(mem_clust_by_r,file_name,quote=FALSE,row.names=FALSE)

### original cluster frequency script

## Generate analysis value data for cluster frequency plotting

# Read source file to establish number of nodes (entries in the file)
# subtract 1 for header line
total_nodes <- nrow(data) -1

# Create matrix of number of clusters and number of nodes (max, mean, median), for each radius
# Columns:  Radius, number of clusters, max cluster size, mean cluster size, median cluster size
col_list <- cbind('radius','num_clust','max_clust_size','max_normalized','mean_clust_size','median_clust_size')
n_cols <- length(col_list)

n_rows <- length(radius_values)
analysis_by_radius <- matrix(, nrow = n_rows, ncol = n_cols)
colnames(analysis_by_radius) <- col_list
row <- 1

for(i in loop_count)
{
  radius <- radius_values[i]
  # Reads in data for each percolation radius that has been computed, and if appropriate by minimum cluster size
  # This lists each node and the id of the cluster to which it is a member, for this radius
  ClstRad_col <- paste("ClstRad",radius,sep="")
  member_cluster <- mem_clust_by_r[c("PlcIndex", ClstRad_col)]
  names(member_cluster) <- c('node','cluster')
  member_cluster <- na.omit(member_cluster)
  
  # number of nodes in clusters, i.e. excluding single nodes
  total_clustered_nodes <- nrow(member_cluster)
  # Ranks clusters by size, i.e. number of nodes for each cluster
  frequency_clusters <- data.frame(table(member_cluster$cluster))
  names(frequency_clusters) <- c('cluster','number_of_nodes')
  ranked_clusters <- frequency_clusters[order(frequency_clusters$number_of_nodes, decreasing=TRUE),]
  # Convert cluster column from factor to numeric, else problems later
  ranked_clusters$cluster <- as.numeric(as.character(ranked_clusters$cluster))
  number_of_clusters <- nrow(ranked_clusters)
  
  # The number of nodes in the first ranked cluster is maximum
  max_nodes <- ranked_clusters$number_of_nodes[1]
  max_normalized <- max_nodes/total_nodes
  mean_nodes <- mean(ranked_clusters$number_of_nodes)
  median_nodes <- median(ranked_clusters$number_of_nodes)
  
  analysis_by_radius[row,] <- cbind(radius,number_of_clusters,max_nodes,max_normalized,mean_nodes,median_nodes)
  row <- row + 1
  
}

analysis_by_radius <<- as.data.frame(analysis_by_radius)
analysis_by_radius
# Save the analysis data as a csv file

file_name <- paste(file.path(path_results,"analysis_by_radius.csv"))
write.table(analysis_by_radius,file=file_name,col.names=TRUE,sep=",",row.names=FALSE,quote=FALSE)

}
