#' CREATE NODES LIST SCRIPT
#' Original code from Elsa Arcaute, CASA, UCL 28.2.15; extensively modified by Simon Maddison.
#' Cite as:	arXiv:1504.08318 [physics.soc-ph]
#' This is now version 3, incorporating corrections, performance and style improvements 10/5/15
#' This code creates the file nodes_list_d.txt containing the list for all the distances between pairs,
#' that is to say a partial matrix of distances between nodes.
#' Distances are computed for each node with each successive node, working through the list.
#' For ca. 3500 nodes this is no longer practicable, as computation become too long and files too big.
#' This uses an approach that creates file with only distances below a limit, e.g. 100km.
#' This will then not be a fully interconnected network;
#' 	need to ensure that this does not adversely affect clustering computation.
#' Note original used the bigmemory library but this is not available for windows, so all refs removed.
#' Version 3.1: 29/3/16 - add check for null values, display and strip out, also update on duplicate coordinates.
#' 	This data written out to text files.
#' Version 4.1: 4/4/2016 -
#' NOTE: problems with a matrix with #'REF entries in .csv -
#'	No check for this, and caused problems with matrix manipulation
#' Removed conversions to and from matrix data format for data, from original source. Not necessary
#' Version 4.2: 10/4/16 - Minor error in duplicate computation corrected. Old ref to matrix_xy; note also
#'     that input data now has additional Lat Long fields added, which are not used here
#' Version 4.3: 6/5/16 - Change index name in source file to bring within character limit in ArcGIS -> PlcIndex
#' Version 4.4: 8/5/16 - Change of source text file format to include map data
#' Version 4.5: 18/5/16 - Change to include output PlcIndex as csv, cleaned up list with only those processed
#' Version 4.6: 07/11/17 - Change working directory for updated machine
#' Version 4.7: 03/12/18 - add in unit parameter to select km or m, in "radius_values.txt"
#'	NOTE: values are assumed to be in metres, unit will factor this in calculations and charts by the unit value
#'	so that a unit value of 1 shows through as metres, and a value of 1000 as km. Other values accomodated

## Sophie's trial to make functions here.

percolate2 <- function(data, radius_values, limit, radius_unit, upper_radius, lower_radius, step_value) {
path_source <- paste(getwd(),"/source_data",sep="")
path_results <- paste(getwd(),"/working_data",sep="")

ptm <- proc.time()
#' For computation time

#' List any entries with null values and write to file
data_NA <- data[rowSums(is.na(data)) > 0,]
print(paste("Entries in source file with one or more null values: "))
data_NA
file_name <- paste(path_results,"/","null_entries.txt",sep="")
write.table(data_NA, file_name, row.names=FALSE)
#' Remove rows with null values from data
data <- na.omit(data)

#' Remove points that are superimposed and keep only the first ID
#' - This removes one of two sites that are very close to each other
#' Determined on basis of x y coordinates
#' Write list to file
duplicate_xy <- data[duplicated(data[,2:3]),]
data_unique <- data[!duplicated(data[,2:3]),]
print(paste("Number of removed superimposed points: ",(nrow(data)-nrow(data_unique))))
duplicate_xy
file_name <- paste(path_results,"/","duplicate_entries.txt",sep="")
write.table(duplicate_xy, file_name, row.names=FALSE)

x_vec <- data_unique$Easting
y_vec <- data_unique$Northing
ID <- data_unique$PlcIndex

#' Write file with list of PlcIndex used
PlcIndex_list <- matrix(ID,ncol=1)
file_name <- paste(path_results,"/","PlcIndex.csv",sep="")
write.table(data_unique$PlcIndex, file_name, row.names=FALSE, col.names="PlcIndex")

#' Number of points/nodes in file with no duplicates
n <- length(ID)
print(paste('number of points: ',n))

#' Create matrix of internodal distances
#' Columns: NodeId1, NodeId2, distance 1-2
col_list <- cbind('ID1','ID2','d12')
ni <- n-1
nj <- n
n_rows <- n*(n-1)/2
nodes_list <- matrix(, nrow = n_rows, ncol = 3)
colnames(nodes_list) <- col_list

row <- 0
i <- 1
for (i in 1:ni)
{
	j1 <- i+1
	for (j in j1:nj)
	{
		#' Compute distance between nodes - pythagoras
		#'  for full grid references this is in units of 1m
		d <- sqrt((abs(x_vec[i]-x_vec[j]))^2+(abs(y_vec[i]-y_vec[j]))^2)
		#' to give distance in m*unit, rounded to 2 decimal places; this also reduces file size
		d <- d/radius_unit #' this factors the values by the unit. 1 gives metres, 1000 gives km
		d <- round(d,2)
		#' Include only if less than limit of distances to be included
		if(d < limit)
			{row <- row+1
			nodes_list[row,] <- cbind(ID[i],ID[j],as.numeric(d))
		}
	}
}
t1 <- proc.time() - ptm
print('loop computed')
print(t1/60)

#' Output file name and location writes to a text file
file_name <- paste(path_results,"/","nodes_list_d.txt",sep="")
#' Remove the unused rows in the matrix
nodes_list <- nodes_list[-(row+1:n_rows),]
m_nodes <- as.matrix(nodes_list)
#' need to write this WITHOUT the row number
write.table(m_nodes, file_name, row.names=FALSE)

t2 <- proc.time() - ptm
print('matrix copied')
print(t2)



## ab hier clustering_script



mem_clust_by_r <- as.data.frame(data_unique$PlcIndex)
names(mem_clust_by_r)[1] <- "PlcIndex"



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

data_file <- paste(path_results,"/","nodes_list_d.txt",sep="")

# The data table of nodes and internode distances is a Text file, with headers
matrix_IDs_distance <- read.table(data_file,header=TRUE)
# Columns are: node Id 1, node Id 2, distance between them. Note that this is generated
#  with a limit to the maximum distance to reduce overall matrix size, and hence
#  creates a partial matrix

t <- as.vector(proc.time() - ptm)[3]
print(paste('time to get matrix in mins',t/60))
ptm2 <- proc.time()[3]

# Define range of percolation radius and step value to progressively reduce it
# Radius in defined unit value, as per radius_unit in radius_values.txt
radius_values <- seq(upper_radius,lower_radius,by=-step_value)

# Changes to accomodate non-integer radius values
radii_count <- length(radius_values)
loop_count <- seq(radii_count,1,by=-1)

# Repeat for each radius value

for (i in loop_count)
{
  print(i)
  radius <- radius_values[i]
  print(radius)
  # Create sub-matrix such that all internode distances d12<=radius
  matrix_radius <- matrix_IDs_distance[matrix_IDs_distance[,3]<=radius,]
  # Create graph (note that characters and numerics are treated differently)
  matrix2 <- matrix_radius[,1:2]
  matrix2[,1] <- as.character(matrix2[,1])
  matrix2[,2] <- as.character(matrix2[,2])
  # This creates a graph from the sub-matrix matrix2
  # Directed means that the edges will only be 'counted' once (tbc)
  # In order to interpret matrix2 as a matrix added 'as.matrix'
  g <- graph.edgelist(as.matrix(matrix2), directed=TRUE)
  
  #take subcomponents - description of how this works
  #http://stackoverflow.com/questions/20725411/finding-strong-and-weak-clusters-and-their-membership-in-r
  
  # Identifies clusters in the graph; creates list of nodes and associated cluster id
  # weak refers to the mechanism used to generate the clusters and relates to compuational efficiency
  # Note that this does not include clusters of 1 node, so the counts are not really meaningful at the lower limit
  member_of_cluster_id <- clusters(g, mode="weak")$membership
  # V processes vertices of graph
  m <- cbind(V(g)$name,member_of_cluster_id)
  new_col_name <- paste("ClstRad",radius,sep="")
  colnames(m) <- c("PlcIndex",new_col_name)
  mem_clust_by_r <-  merge(mem_clust_by_r, m, by="PlcIndex", all.x=TRUE)
  
  # following converts NA to zero, but not used
  #	mem_clust_by_r[,new_col_name] <- as.numeric(as.character(mem_clust_by_r[,new_col_name]))
  #	mem_clust_by_r[is.na(mem_clust_by_r)] <- 0
  
  ptm1 <- proc.time()[3]
  t_print <-  as.vector(ptm1-ptm2)/60
  print(paste("for radius=",radius,"it took to compute",sprintf("%4.2f",t_print),"mins"))
  
}

file_name <- paste(path_results,"/","member_cluster_by_radius.csv",sep="")
# Write file without row names
write.csv(mem_clust_by_r,file_name,quote=FALSE,row.names=FALSE)


}
