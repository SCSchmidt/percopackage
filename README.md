
<!-- README.md is generated from README.Rmd. Please edit that file -->

# percopackage

Percolation Analysis as implemented here is a 2D point pattern analysis
technique for identifying clusters of any size and form. It recognises
noise and produces repeatable result.

Recently the algorithm has been used in geography to identify
metropolitan areas, based on population density. The City Clustering
Algorithm (CCA) has been developed out of percolation theory by
Rozenfeld et al. (2008) based on distance within a cellular lattice; it
has been further developed by Rozenfeld et al. (2011) to use the
Euclidean distance between points. The algorithm implemented here is
based on Arcaute et al. (2016), who adopted the technique for defining
urban areas, using the density of street interconnections rather than
population.

A cluster consists of at least two points. Around each point in the
given set a defined distance threshold is drawn as a radius and all
points falling within this threshold are connected to the cluster. The
test is then re-applied for each of these neighbours in turn, and any
further points meeting this criterion are also part of the cluster. This
technique can be applied at any scale, from the molecular to the
geographical and beyond. The implementation here can therefore use
radius values of decimal meters and kilometers.

The algorithm has been implemented as an exploratory tool, meaning that
in the percolate-function a set of radii is defined (values from, values
to, and step values), for which the percolation is run. The results can
be compared via maps (mapClusters-function) and graphs
(plotClustFreq-function).

Two archaeological case studies using this algorithm and software
package are to be published by Maddison & Schmidt 2020.

## Percolation Analysis via the R package “percopackage” – Workflow and Descriptions

Percopackage is an R-package, which takes scripts developed by Simon
Maddison and Elsa Arcaute to create three functions, which perform

1)  the actual percolation analysis (function called “percolate”),
2)  maps the clusters created by percolation (function “mapClusters”)
    and
3)  creates three different plots of the analysis results (function
    “plotClustFreq”)

The percolate function does the actual computing of the percolation
analysis and exports different tables on which the following functions
operate. Two different data inputs are possible: A set of points in
Euclidian space, for which a distance matrix will be calculated OR a
distance matrix, which may for example have been calculated using least
cost paths. If a distance matrix is provided, the point set is still
needed for the mapping function.

### Directory Structure

The functions create subdirectories of the working directory the user is
working in (using getwd() ), therefore we describe here the relative
structure, which is being set up:

  - /analysis\_results – output plots from the analysis, now principally
    the cluster size vs radius plots and the csv of the analyis results
    (described further down)

  - /maps – outputs from the mapping, showing clusters and nodes
    overlaid on given outline shape files, as well as output shape files
    with the point data (shape file called “map\_name.shp”)

  - /working\_data – working data generated by the programs and used to
    pass intermediate data between them.

### Percolate - function

The input needed for the percolate-function is as follows:

*percolate(data, distance\_table = NULL,upper\_radius, lower\_radius, step\_value, limit, radius\_unit)*

**data** = A data.frame, with the format below:

PlcIndex,Easting,Northing

1,350350,233050

2,354700,260200

3,358700,238900

*PlcIndex* is an assigned index number to distinguish the points. This
does not need to be sequential and can for example be drawn from a
source database.

*Easting* and *Northing* are self-explanatory, and are to the metre.
(Note that additional columns, if any, will be ignored).

**distance\_table** is set to be NULL as a default, which means, from
the data input (coordinates) a distance matrix is calculated using
Euclidean distance. If not NULL, the data frame given here needs to
adhere to this format:

ID1, ID2, d12

1,2,12.1

1,3,14.2

2,3,2.9

*ID1* and *ID2* are point IDs (i.e. the PlcIndex values) and

*d12* denotes the weighted distance between them, if a non-euclidian
distance is whished to be used.

**limit** is an integer, the value above which distances will not be
calculated between sites.

**radius\_unit** is either 1 for meter or 1000 for km for all input
radii.

**upper\_radius** is an integer, the upper value of the radius range to
be used.

**lower\_radius** is an integer, the lower value of the radius range to
be used.

**step\_value** is numeric or integer and the step value to be used
between lower and upper radius.

**Exported** are in the folder *working\_data*: tables (csv) with the
input-data, the list of PlcIndex, null\_entries and duplicate\_entries
to be used by following functions.

In the folder *analysis\_results* the tables (csv) analysis\_by\_radius
and member\_cluster\_by\_radius are exported. They will be used by the
following functions as well, but may be useful for other applications as
well, therefore they are considered analysis results.

#### Note

The function creates the Inter-site distance matrix – given the dataset
of sites and XY coordinates, this computes a partial matrix of
inter-site distances. The limit is set by the **limit** parameter in the
function call. The inter-site distance is computed using Pythagoras’
theorem. The output file is located in */working\_data* and is called:
nodes\_list\_d.txt. See below for examples and explanation of the limit
parameter.

It also creates a file listing the nodes that have been used in creating
this list, which omits duplicates and null values etc. This is called
PlcIndex.csv.

Cluster extraction is the second step: This program identifies all
clusters for a specific radius. It steps through a range of radii, the
range being set by input-parameters. The data is output as a text csv
file called: member\_cluster\_by\_radius.csv. This program is the heart
of the analysis; the method of creating a graph where nodes are
connected when the distance between them is the given radius or less,
then extracting each cluster, remains as written by Elsa Arcaute.

Note that for large numbers of points, and large number of selected
radius values, the computation time can be quite long. We recommend
experimenting initially with only a small number of steps, e.g. with a
small radius range, or with a large step value, in order to see which
range of values is the most interesting.

For Hillforts in Britain, Simon Maddison ran the analysis with 1 km
steps between 2 and 40km, but also ran an analysis for 3-6km with 0.1km
steps to investigate clusters in high density areas of Scotland and
Wales. To run with 0.1km steps over the full range would not only take
excessive time but also produce large numbers of very similar plots.

Duplicate and null entries are also listed in files with these names,
but should have no values if the data is clean. (Historically there was
a lot of cleaning to do, which is why this runs as a separate program).

#### Example:

*percolate(data = database,,upper\_radius = 40, lower\_radius = 2, step\_value = 1, limit = 50, radius\_unit = 1000)*


A coordinates dataset is the database input. In the example above,
percolation analysis will run for radius values between 40 and 2 km,
with a step value of 1km (determined by the radius\_unit parameter, set
to 1000). When computing the inter-site distances, values greater than
50km will not be stored. This restricts the intermediate file size and
improves performance. (There is no point in computing the distance
between sites in Scotland and Cornwall for example, when all sites are
within a single cluster above a radius of say 40km).

### mapClusters - function

This function will map the clusters created by percolation. It needs the
following input:

*mapClusters(shape, map\_name, source\_file\_name, dpi=300)*

**shape** is a SpatialPolygonDataFrame (sp), which will be used for the
plotting as a background map, and determining the coordinate system.

**map\_name** is a string, which will be used to title the maps.

**source\_file\_name** is a string, which will be used to label the maps
and describes the source data used.

**dpi** is the output file resolution. The default is set to 300.

**Exported** maps are in the folder /maps: As many png map files as
there are different radii values used in the percolation analysis. Also
a shp-file with the point data.

#### Note

Mapping clusters – this function maps clusters for each radius on an
outline map overlay. The same data can for example be mapped on a
coastal outline, an outline with modern counties, or an outline with
historic counties. The outline is defined by an ESRI shape file, which
needs to be read, and passed to the function as an appropriate data
frame. The outputs generated were originally a multi-page pdf file, but
this proved unwieldy, and it was changed to generate png files, which
are much easier to import into Microsoft Word and other applications, as
well as to display. Point data is also exported as an ESRI shape file
for import directly into ArcGIS etc. Outputs are placed in the /maps
sub-directory.

Clusters are ranked in order of size for each different threshold
radius. The top 15 clusters are assigned a colour code with Red the
largest, Blue the next and so on. A legend is plotted with the assigned
cluster number against the colour. The plots also include the map name
and the name of the source file used to generate the data, as well as
the radius value.

The map projection is extracted from the given shape file and used for
configuring the point plots and output shape files. This has not been
extensively tested for different projection/ coordinate systems, and may
need further development/modification for different input shape files.
If necessary experiment by omitting the shape file plot by commenting
out, and simply plot the points, in the first instance. Points and shp
need to have the same coordinate system and projection.

#### Example

*mapClusters(data = database, shape = outline\_GB, map\_name = “Britain
Iron Age Hillforts”, source\_file\_name = “Hillfort Sites”)*

The output-maps will be called e.g. “percolation\_plots\_Britain Iron
Age Hillforts\_rad\_40.png” (number changing regarding to which
percolation radius is depicted). On top of the map, the title will be:
“Britain Iron Age Hillforts” and it will have a subtitle with a string
like: “Source File: Hillfort Sites, percolation distance: 40 km”.

## plotClustFreq - function

Percolation cluster size and rank analyses – this program processes the
cluster data generated by the cluster extraction program to give various
different frequency and ranking plots.

*plotClutsFreq(source\_file\_name)*

**source\_file\_name** is a string, which will be used to label the
plots and describes the source data set. If already input in mapClusters
this doesn’t need to be done again.

Three png-files are generated:

1)  radius to maximum cluster size,
2)  radius to mean cluster size and
3)  radius to normalized max. cluster size

and they will include the information of the source file as a subtitle.

At the moment they are stored in the /analysis\_results directory as
“radius\_to\_max\_cluster\_size”, “radius\_to\_mean\_cluster\_size”
and “radius\_to\_norm\_max\_cluster\_size”.

#### Note

This program evolved significantly as various different statistical
analyses were experimented with, under the guidance of Mark Lake. The
final version reflects the latest work by Elsa Arcaute. The outputs are
placed in the /analysis\_results sub-directory. The most useful plot is
the cluster transition plot, showing the maximum cluster size
(normalized) vs. radius, which shows the transitions.

### CONTACT DETAILS:

Simon Maddison <simon.maddison@btinternet.com> +44 797 089 2171

Sophie C. Schmidt <s.c.schmidt@uni-koeln.de>

# Bibliography

Arcaute, E., Molinero, C., Hatna, E., Murcio, R., Vargas-Ruiz, C.,
Masucci, A P & Batty, M 2016 Cities and regions in Britain through
hierarchical percolation. Royal Society Open Science, 3 (4): 150691.DOI:
10.1098/rsos.150691

Rozenfeld, H. D., Rybski, D., Andrade, J. S., Batty, M., Stanley, H. E.
& Makse, H. A. 2008, Laws of population growth. Proceedings of the
National Academy of Sciences of the United States of America, 105 (48):
18702.DOI: 10.1073/pnas.0807435105

Rozenfeld, H. D., Rybski, D., Gabaix, X. & Makse, H. A. 2011, The Area
and Population of Cities: New Insights from a Different Perspective on
Cities. American Economic Review 101 (5): 2205-25.

### Citation

Please cite this compendium as:


> Authors, (2020). *Title of compendium*. Accessed 21 Apr 2020. Online
> at <https://doi.org/xxx/xxx>

### Installation

You can install percopackage from github with:

``` r
# install.packages("devtools")
devtools::install_github("SCSchmidt/percopackage")
```

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

### Contributions

We welcome contributions from everyone. Before you get started, please
see our [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.
