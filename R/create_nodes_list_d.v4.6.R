# CREATE NODES LIST SCRIPT
# Original code from Elsa Arcaute, CASA, UCL 28.2.15; extensively modified by Simon Maddison.
# Cite as:	arXiv:1504.08318 [physics.soc-ph]
# This is now version 3, incorporating corrections, performance and style improvements 10/5/15
# This code creates the file nodes_list_d.txt containing the list for all the distances between pairs, 
# that is to say a partial matrix of distances between nodes.
# Distances are computed for each node with each successive node, working through the list.
# For ca. 3500 nodes this is no longer practicable, as computation become too long and files too big.
# This uses an approach that creates file with only distances below a limit, e.g. 100km. 
# This will then not be a fully interconnected network; 
# 	need to ensure that this does not adversely affect clustering computation.
# Note original used the bigmemory library but this is not available for windows, so all refs removed.
# Version 3.1: 29/3/16 - add check for null values, display and strip out, also update on duplicate coordinates.
# 	This data written out to text files.
# Version 4.1: 4/4/2016 -
# NOTE: problems with a matrix with #REF entries in .csv - 
#	No check for this, and caused problems with matrix manipulation
# Removed conversions to and from matrix data format for data, from original source. Not necessary
# Version 4.2: 10/4/16 - Minor error in duplicate computation corrected. Old ref to matrix_xy; note also 
#     that input data now has additional Lat Long fields added, which are not used here
# Version 4.3: 6/5/16 - Change index name in source file to bring within character limit in ArcGIS -> PlcIndex
# Version 4.4: 8/5/16 - Change of source text file format to include map data
# Version 4.5: 18/5/16 - Change to include output PlcIndex as csv, cleaned up list with only those processed 
# Version 4.6: 07/11/17 - Change working directory for updated machine

setwd("D:/Iron_Age_Hillforts/Percolation")
# paths
path_source <- paste(getwd(),"/source_data",sep="")
path_results <- paste(getwd(),"/working_data",sep="")

# Read in source file name defined in file source_file.txt
file_name <- paste(path_results,"/","source_file.txt",sep="")
source_files <- read.csv(file_name,header=TRUE,stringsAsFactors=FALSE)
source_file_name <- source_files$source_file
# Edit this file for different input file name

# Original file has been generated from Hogg Index spreadsheet
# Irish file was derived from material from the Atlas, ca. 14 May 2015
# Later analyses run in the same way on Atlas data March 2016 and on published data November 2017

source_file <- paste(path_source,"/",source_file_name,sep="")

# Read in limit value from source file
# Edit the source file radius_values.txt to change this.
file_name <- paste(path_results,"/","radius_values.txt",sep="")
radius_values <- read.csv(file_name,header=TRUE)
limit <- radius_values$limit

ptm <- proc.time()
# For computation time

data <- read.csv(source_file)
# This file needs three columns: PlcIndex, Easting, Northing in this order.
# Additional fields Lat Long not used for this program.
# Note that the Easting and Northing need to be UK grid references with the Alphabetic 
#  grid square converted to a numeric prefix for each
# List any entries with null values and write to file
data_NA <- data[rowSums(is.na(data)) > 0,]
print(paste("Entries in source file with one or more null values: "))
data_NA
file_name <- paste(path_results,"/","null_entries.txt",sep="")
write.table(data_NA, file_name, row.names=FALSE)
# Remove rows with null values from data
data <- na.omit(data)

# Remove points that are superimposed and keep only the first ID 
# - This removes one of two sites that are very close to each other
# Determined on basis of x y coordinates
# Write list to file
duplicate_xy <- data[duplicated(data[,2:3]),]
data_unique <- data[!duplicated(data[,2:3]),]
print(paste("Number of removed superimposed points: ",(nrow(data)-nrow(data_unique))))
duplicate_xy
file_name <- paste(path_results,"/","duplicate_entries.txt",sep="")
write.table(duplicate_xy, file_name, row.names=FALSE)

x_vec <- data_unique$Easting
y_vec <- data_unique$Northing
ID <- data_unique$PlcIndex

# Write file with list of PlcIndex used
PlcIndex_list <- matrix(ID,ncol=1)
file_name <- paste(path_results,"/","PlcIndex.csv",sep="")
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

row <- 0
i <- 1
for (i in 1:ni)
{
	j1 <- i+1
	for (j in j1:nj)
	{
		# Compute distance between nodes - pythagoras
		#  for full grid references this is in units of 1m
		d <- sqrt((abs(x_vec[i]-x_vec[j]))^2+(abs(y_vec[i]-y_vec[j]))^2)
		# to give distance in km, rounded to 2 decimal places; this also reduces file size
		d <- d/1000
		d <- round(d,2)
		# Include only if less than limit of distances to be included
		if(d < limit)
			{row <- row+1
			nodes_list[row,] <- cbind(ID[i],ID[j],as.numeric(d))
		}
	}
}
t1 <- proc.time() - ptm
print('loop computed')
print(t1/60)

# Output file name and location writes to a text file
file_name <- paste(path_results,"/","nodes_list_d.txt",sep="")
# Remove the unused rows in the matrix
nodes_list <- nodes_list[-(row+1:n_rows),]
m_nodes <- as.matrix(nodes_list)
# need to write this WITHOUT the row number
write.table(m_nodes, file_name, row.names=FALSE)

t2 <- proc.time() - ptm
print('matrix copied')
print(t2)