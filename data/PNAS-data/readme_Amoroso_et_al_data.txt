Description of supplemental data

Bottom trawl-fishing footprint on the world’s continental shelves: data format and usage

Background:

This project collected and analyzed high-resolution trawling activity data on the continental shelves and slopes (0-1000m depth). High-resolution grids and a series of assumptions about the fine-scale distribution of trawling activity were used to quantify the distribution and intensity of trawling in 24 regions where =70% of total fishing activity was recorded and 8 regions where <70% of activity was recorded. 

Data files:

Two data files are provided. First, an S4 (R) object containing high resolution trawling activity data and information on regional and depth boundaries used in the analyses (see main paper). Second, a csv file containing predicted trawling footprint and associated uncertainty (based on the uniform assumption, see main paper) as a function of the regional Swept Area ratio. 


Attributes of data files:

"TBPdata" is an S4 object. It can be loaded directly into an R workspace for viewing and further processing of summary statistics per region, shapefiles for the regions covered in the study, raster files of the depth strata and data on trawling intensity in each cell where trawling activities were recorded. Full detail of the analytical methods used to process these data are provided in the main paper and supporting information. 


Accessing data in "TBPdata"

"TBPdata" can be loaded into an R workspace with the command:

load(file.path(fileDirectory, "TBPdata.Rdata")

To fully access and manipulate the data contained in TBPdata the following R packages need to be installed "rgeos", "raster", "rgdal", "proj4". 


The TBPdata object contains the following information:

TBPdata@RegionName : vector with the names of the regions used in the study

# e.g. to obtain names of regions
TBPdata@RegionName
 [1] "Adriatic Sea (GFCM 2.1)"          "Aegean Sea (GFCM 3.1)"           
 [3] "Aleutian Islands"                 "Argentina"      

TBPdata @RegionShapefile: a named list containing a SpatialPolygonsDataFrame of each region. Please note that all the elements of the list do not have the same geographic projection. Individual elements can be extracted with the region name. 

# e.g. to view projections and associated information by region
TBPdata@RegionShapefile[["Aegean Sea (GFCM 3.1)"]]
class       : SpatialPolygonsDataFrame 
features    : 1 
extent      : 22.5243, 28.99953, 34, 41.00941  (xmin, xmax, ymin, ymax)
coord. ref. : +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 
variables   : 2
names       : ID, AreaName 
min values  : 218, Aegean Sea (GFCM 3.1) 
max values  : 218, Aegean Sea (GFCM 3.1)

The projections can be obtained with projection () 

# e.g. to define projection by region
projection(TBPdata@RegionShapefile[["Aegean Sea (GFCM 3.1)"]])

The polygons can be viewed with plot ()

# e.g. to view a regional polygon (regional boundaries)
plot(TBPdata@RegionShapefile[["Aegean Sea (GFCM 3.1)"]])

TBPdata @RegionDephtStrata: a named list containing a RasterLayer of the depth strata in each study area. The grid corresponds to the highest resolution grid used in each region (note that some grids were built using geographic coordinates and others using metres). The 'value' of the cells is 1 or 2, corresponding to depth strata of 0-200 m or 200-1000m respectively. The projection of the RasterLayer is the same than the SpatialPolygonsDataFrame. 

The depth strata raster for a defined region can be accessed and plotted, and viewed within the regional polygon.

# e.g. to access information on a RasterLayer for a given region
TBPdata@RegionDepthStrata[["Aegean Sea (GFCM 3.1)"]]
class       : RasterLayer 
dimensions  : 421, 389, 163769  (nrow, ncol, ncell)
resolution  : 0.01666667, 0.01666667  (x, y)
extent      : 22.51667, 29, 34, 41.01667  (xmin, xmax, ymin, ymax)
coord. ref. : +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 
data source : in memory
names       : ETOPO1_Ice_g_geotiff 
values      : 1, 2  (min, max)

# e.g. to view depth strata in a given region
plot(TBPdata@RegionDepthStrata[["Aegean Sea (GFCM 3.1)"]])

# e.g. to view depth strata in a given region alongside the regional polygon (regional boundaries)
plot(TBPdata@RegionDepthStrata[["Aegean Sea (GFCM 3.1)"]])
plot(TBPdata@RegionShapefile[["Aegean Sea (GFCM 3.1)"]],add=T)


TBPdata @ SummaryByRegion: a data frame containing an aggregated summary of regional and trawling activity statistics. Refer to the main paper and supporting information for definitions of the alternate assumptions used to estimate footprint. Variables can be listed with matrix()

# e.g. to list regional and trawling activity statistics
matrix(names(TBPdata @SummaryByRegion),ncol=1)
      [,1]                                       
 [1,] "RegionName"                               
 [2,] "AreaKm2_0_200"                            
 [3,] "AreaKm2_200_1000"                         
 [4,] "SweptAreaKm2_0_200"                       
 [5,] "SweptAreaKm2_200_1000"                    
 [6,] "AreaTrawledKm2_CellAssumption_0_200"      
 [7,] "AreaTrawledKm2_CellAssumption_200_1000"   
 [8,] "AreaTrawledKm2_PoissonAssumption_0_200"   
 [9,] "AreaTrawledKm2_PoissonAssumption_200_1000"
[10,] "AreaTrawledKm2_UniformAssumption_0_200"   
[11,] "AreaTrawledKm2_UniformAssumption_200_1000"
[12,] "AreaKm2_Sand_0_200"                       
[13,] "AreaKm2_Sand_200_1000"                    
[14,] "AreaKm2_Mud_0_200"                        
[15,] "AreaKm2_Mud_200_1000"                     
[16,] "AreaKm2_Gravel_0_200"                     
[17,] "AreaKm2_Gravel_200_1000"                  
[18,] "AreaKm2_0_1000"                           
[19,] "SAR"                             

TBPdata@CellData: a named list containing a dataframe with data on trawling activity and habitat type for each cell in the region for which trawling activity was recorded. Owing to confidentially restrictions requested by management and enforcement agencies we cannot provide the precise geographic position of each cell. The values are averaged for the duration of the study (typically three years, see supporting information for study years by region). The dataframe is ordered swept area. 

# e.g. to view cell by cell activity data for each region
head(TBPdata@CellData[["Aegean Sea (GFCM 3.1)"]], 5, add=T)
        CellAreaKM2 DepthStrata SweptAreaKm2 AreaTrawledKm2_PoissonAssumption AreaTrawledKm2_UniformAssumption  %Gravel    %Sand      %Mud
2025238    2.610402           2     194.7651                         2.610402                         2.610402 65.00000 35.00000   0.00000
2111380    2.767344           2     153.9223                         2.767344                         2.767344 43.33333 30.83333  25.83333
2112123    2.769093           2     135.9922                         2.769093                         2.769093 48.75000 33.25000  18.00000
2027175    2.614890           1     124.6086                         2.614890                         2.614890  0.00000  0.00000 100.00000
2111381    2.767344           2     123.2242                         2.767344                         2.767344 43.33333 30.83333  25.83333

# e.g. to access all cell by cell activity data for each region
TBPdata@CellData[["Aegean Sea (GFCM 3.1)"]]


Accessing data in "project_footprint_Amoroso_et_al"

project_footprint_Amoroso_et_al is a csv file providing estimates of footprint (assuming uniform spread of trawling activity in grid cells, see main paper) as a function of the regional swept area ratio (0.1 increments from 1 to 12). Prediction intervals were estimated from 10^6 simulations for each value of SAR. The columns in this file are:
SAR: swept area ratio 
lower: lower (5%) prediction interval for footprint
median: median estimate of footprint
upper: upper (95%) prediction interval for footprint

