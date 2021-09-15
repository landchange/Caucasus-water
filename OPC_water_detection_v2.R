library(reticulate)
library(raster)
library(sf)
library(anchors)
library(dplyr)
library(rpart)
library(gdata)
library(pROC)
library(pscl)
library(randomForest)
library(InformationValue)
library(lubridate)
library(googledrive)
library(geojsonio)
library(rgee)

# Initilize Earth Enginge session
ee_Initialize()
ee_Initialize(email = 'owsleybc@gmail.com')
# ee_check()


###
# Create list of years to include
# datesList <- c('2000-01-01','2001-01-01','2002-01-01','2003-01-01','2004-01-01',
#                '2005-01-01','2006-01-01','2007-01-01','2008-01-01','2009-01-01',
#                '2010-01-01','2011-01-01','2012-01-01','2013-01-01','2014-01-01',
#                '2015-01-01','2016-01-01','2017-01-01','2018-01-01','2019-01-01',
#                '2020-01-01')

# Beginning and Ending Date
dateBegin <- '2000-01-01'
dateEnd <- '2001-01-01'

monthBegin <- 5
monthEnd <- 10

clouds <- 30
# bands <- c('B3', 'B6') # Landsat 8
bands <- c('B2', 'B5') # Landsat 5, 7

slope <- 10

# Set paths and rows
# paths <- c(166, 167, 168, 169, 170, 171, 172, 173)
paths <- c(166:173)
rowsList <- list(c(032:033), c(031:033), c(031:033), c(031:033), 
                 c(030:032), c(030:031), c(030:031), c(029:030))

# Study region shapefile
caucasusAdmin <- ee$FeatureCollection('users/owsleybc/caucasus/caucasusAdmin1')

# Landsat Images
# landsatCollection <- 'LANDSAT/LE05/C01/T1_SR'
landsatCollection <- 'LANDSAT/LE07/C01/T1_SR'
# landsatCollection <- 'LANDSAT/LC08/C01/T1_SR'

# JRC Images
jrcCollection <- 'JRC/GSW1_2/MonthlyHistory'

# Validation data
validationData <- 'users/yajunzhang2/caucasus/validations/total_caucasus_valid'
# bring in validation data and remove NA values
validation_ee <- ee$FeatureCollection(validationData)
validationPointsSf <- ee_as_sf(validation_ee, "validationPointsCSV.csv", 
                               maxFeatures = 10000)
# transform data into binary format
validationPointsSf <- replace.value(validationPointsSf, 'waterClassb', 
                                    from = 1, to = as.integer(0))
validationPointsSf <- replace.value(validationPointsSf, 'waterClassb', 
                                    from = 3, to = as.integer(1))
validationPointsSf <- replace.value(validationPointsSf, 'validated', 
                                    from = 1, to = as.integer(0))
validationPointsSf <- replace.value(validationPointsSf, 'validated', 
                                    from = 3, to = as.integer(1))
validationPointsSf <- replace.value(validationPointsSf, 'validated', 
                                    from = 2, to = NA)
validationPointsSf <- na.omit(validationPointsSf)
# change from sf object to an EE object
validation_ee <- sf_as_ee(validationPointsSf)

# create column names for OPC and JRC error tables and statistics
colNames <- c('date', 'roc_auc', 'opt_cutoff', 'prob_0_0', 'prob_0_1', 'prob_1_0', 'prob_1_1', 
              'specificity', 'sensitivity', 'intercept', 'slope', 'roc_auc_jrc', 'val_0_0', 
              'val_0_1', 'val_1_0', 'val_1_1', 'specificityjrc', 'sensitivityjrc')

# Urban and Mountain Masks
dem <- ee$Image("USGS/SRTMGL1_003")
mountainMask <- ee$Terrain$slope(dem)$lte(slope)

# bring in GMIS dataset???????????????????????????????????????????????????????????????????????????
# GMIS <- 'users/yajunzhang2/caucasus/Urban_mask/GMIS_Urban_mask_clip'

# Create JRC Training Layer: Land and water
jrcImages <- ee$ImageCollection(jrcCollection)$
  filterDate(dateBegin, dateEnd)$
  filter(ee$Filter$calendarRange(monthBegin, monthEnd, 'month'))

# reclass list for JRC dataset validation
var1 <- ee$List(c(0, 1, 2))
var2 <- ee$List(c(-1, 0, 1))


###
# Cloud masking function
maskSR <- function(image) {
  # Get the pixel QA band.
  qa = image$select('pixel_qa')
  # Both flags should be set to zero, indicating clear conditions.
  mask = qa$bitwiseAnd(2)$Or(qa$bitwiseAnd(4))
  return(image$updateMask(mask))
}


###
# j <- 1
for(j in 1:length(paths)) {
  path <- paths[j]
  rows <- rowsList[[j]]
  
  # File path for saving stats
  statsFileName = paste('F:/Caucasus/Water/rgee_output/stats_', path, '.csv', sep = '')
  
  # read in images from Landsat  based on date, conditions and Landsat path/row tiles
  imageCollection <- ee$ImageCollection(landsatCollection)$
    filterDate(dateBegin, dateEnd)$
    filter(ee$Filter$calendarRange(monthBegin, monthEnd, 'month'))$
    filter(ee$Filter$And(ee$Filter$eq('WRS_PATH', path),ee$Filter$inList('WRS_ROW', rows)))$
    filter(ee$Filter$lt('CLOUD_COVER', clouds))
  
  # Dates of images in image collection
  dateList <- ee_get_date_ic(imageCollection)
  # Sort the dates
  dateList <- dateList[order(dateList$time_start), ]
  dates <- as.Date(dateList$time_start)
  dateList <- cbind(dateList, dates)
  uniqueDates <- unique(dates)
  
  
  stats <- data.frame()
  frequencyMaps <- ee$Image(0)
  countMaps <- ee$Image(0)
  # for loop to generate OPC water index maps from images in path/row set
  #i <- 1
  for(i in 1:length(uniqueDates)) {
    month <- month(uniqueDates[i])
    date <- as.character(uniqueDates[i])
    datePlus1 <- as.character(uniqueDates[i] + 1)
    # create path/row set and mosaic
    imageCollection <- ee$ImageCollection(landsatCollection)$
      filterDate(date, datePlus1)$
      filter(ee$Filter$And(ee$Filter$eq('WRS_PATH', path),ee$Filter$inList('WRS_ROW', rows)))$
      map(maskSR)$
      select(bands)
    
    imageInfo <- ee_print(imageCollection)
    if(imageInfo$ic_nimages < 2) next
    
    imageMosaic <- imageCollection$mosaic()
    # Map$addLayer(imageMosaic)
    
    extent <- imageCollection$geometry()$dissolve()
    
    # Select JRC image with corresponding date and clip to Landsat Tile extent
    jrc <- jrcImages$filter(ee$Filter$eq('month', month))$
      first()$
      clip(extent)
    
    #reclassification of JRC image for use in validation of JRC to validation points
    jrc2 <- jrc$remap(var1, var2)
    
    # Calculate MNDWI
    MNDWI <- imageMosaic$normalizedDifference(bands)
    # Map$addLayer(MNDWI)
    
    # Generate stratified random points from JRC for training
    trainingSample <- jrc$stratifiedSample(
      classBand = 'water', 
      classValues = c(0, 1, 2), 
      classPoints = c(0, 750, 750), 
      numPoints = 1500, 
      seed = 9, 
      region = extent, 
      geometries = TRUE)
    
    # extract mndwi values
    trainingPoints <- ee$Image$sampleRegions(MNDWI, 
                                             collection = trainingSample, 
                                             scale = 30)
    
    # convert EE object into SF object
    trainingPointsSf <- ee_as_sf(trainingPoints, "trainingPointsCSV.csv")
    
    # remove NA values
    training <- na.omit(trainingPointsSf)
    
    # add training and rename columns  
    training$id <- seq(1:length(training$nd))
    # training$training <- 1 #??????????????????????????????????????????????????????????????????????????
    # colnames(training) <- c('mndwi', 'water', 'geometry', 'id', 'training')
    colnames(training) <- c('mndwi', 'water', 'geometry', 'id')
    
    # sample landsat tile for MNDWI values
    validationPoints <- MNDWI$sampleRegions(
      collection = validation_ee, 
      scale = 30)
    
    # convert EE object to SF 
    # remove NA values
    validationPointsSf <- ee_as_sf(validationPoints,"validationPointsCSV.csv")
    validationPointsSf$id <- seq(1:length(validationPointsSf$geometry))
    
    # clean data to match rows
    # validation <- validationPointsSf[ ,c('nd', 'validated', 'geometry', 
    #                                      'id', 'Training')]
    validation <- validationPointsSf[ ,c('nd', 'validated', 'geometry', 
                                         'id')]
    # rename columns
    # colnames(validation) <- c('mndwi', 'water', 'geometry', 'id', 'training')
    colnames(validation) <- c('mndwi', 'water', 'geometry', 'id')
    if(max(validation$water) == 0) next
    
    validationDf <- as.data.frame(validation) #?????????????????????????????????????????????????
    
    # # Create Single Training and Validation dataset
    # t_v_points <- bind_rows(training, validation)
    # # separate points into training and validation datasets
    # validation <- t_v_points[which(t_v_points$training == 0), ]
    # training <- t_v_points[which(t_v_points$training == 1), ] ???????????????????
    
    
    # generate logistic regression model using selected index
    mylogit <- glm(as.factor(water) ~ mndwi, 
                   family = binomial(link = "logit"), 
                   maxit = 100, 
                   data = training)
    # summary(mylogit)
    
    # generate predicted value
    predicted <- predict(mylogit, validation, type = "response")
    # plotROC(validation$water, predicted)
    
    roc <- roc(validation$water, predicted)
    # use predicted value and validation dataset for optimal probability cut-off value
    optCutoff_L <- optimalCutoff(validation$water, predicted)[1]
    
    # concordance value
    Concordance(validation$water, predicted)
    
    # generate confusion matrix
    confMatrix <- InformationValue::confusionMatrix(validation$water, predicted, 
                                                    threshold=optCutoff_L)
    # sensitivity
    sensitive <- InformationValue::sensitivity(validation$water, predicted, 
                                               threshold = optCutoff_L)
    # specificity
    specific<- InformationValue::specificity(validation$water, predicted, 
                                             threshold = optCutoff_L)
    
    # create probability map using logistic regression probability equation and water index
    probabilityMap <- MNDWI$expression('1 / (1 + exp(- (s * nd + c)))', 
                                       list('nd' = MNDWI$select('nd'), 
                                            's' = mylogit$coefficients[2], 
                                            'c' = mylogit$coefficients[1]))
    # Map$addLayer(probabilityMap)
    
    # mountain mask
    probabilityMap <- probabilityMap$updateMask(mountainMask)
    # urban mask
    # probabilityMap <- probabilityMap$updateMask(GMIS)
    
    
    ###
    # # validate JRC water data
    # validationPointsJRC <- jrc2$sampleRegions(collection = validation_ee, 
    #                                           scale = 30)
    # 
    # # change class from EE object to sf object
    # validationPointsJRCsf <- ee_as_sf(validationPointsJRC, "validationPointsJRC.csv")
    # # remove '-1' for binary dataset of water/non-water
    # validationPointsJRCsf <- validationPointsJRCsf[validationPointsJRCsf$remapped != -1, ]
    # 
    # # change class to dataframe
    # validationPointsJRCsf <- as.data.frame(validationPointsJRCsf) #????????????????????????????????????????????
    # 
    # validationPointsJRCsf$id <- seq(1:length(validationPointsJRCsf$geometry))
    # # match number of validation points to OPC validation points for comparision
    # validationPointsJRCsf <- merge(x = validationDf, y = validationPointsJRCsf, by = 'id')
    # 
    # # JRC error matrix
    # valMatrix <- InformationValue::confusionMatrix(validationPointsJRCsf$validated, 
    #                                                validationPointsJRCsf$remapped)
    # # JRC ROC
    # rocJRC <- roc(validationPointsJRCsf$validated, 
    #               validationPointsJRCsf$remapped)
    # # JRC specificity
    # sensitiveJRC <- InformationValue::sensitivity(validationPointsJRCsf$validated, 
    #                                               validationPointsJRCsf$remapped)
    # # JRC sensitivity
    # specificJRC <- InformationValue::specificity(validationPointsJRCsf$validated, 
    #                                              validationPointsJRCsf$remapped)
    
    
    ###
    # use optimal cutoff value and cast raster as float for frequency map
    frequencyMap <- probabilityMap$gte(optCutoff_L)$
      unmask()$
      toFloat()#$
    # rename('Frequency_Map')
    
    # use optimal cutoff value and cast raster as float for extent map
    countMap <- probabilityMap$gte(0)$ #????????????????????????????????????????????????????????????????
      unmask()$
      toFloat()#$
      # rename('Count_Map')
    
    # # dataframe of OPC and JRC statistics
    # stat <- data.frame(c(date, roc$auc, optCutoff_L, confMatrix[1, ], confMatrix[2, ], 
    #                      specific, sensitive, mylogit$coefficients[1], mylogit$coefficients[2], 
    #                      rocJRC$auc, valMatrix[1, ], valMatrix[2, ], specificJRC, sensitiveJRC))
    # # combine column names to stat dataframe
    # names(stat) <- colNames
    # stats <- rbind(stats, stat)
    
    # Combine extent and frequency maps with all observations
    frequencyMaps <- frequencyMaps$add(frequencyMap)
    countMaps <- countMaps$add(countMap)
    
  }
  # create csv with error and accuracy metrics
  # write.csv(stats, statsFile) #???????????????????????????????????????????????????????????????????????????????????
  
  
  # ????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  # create frequency map for the path
  waterMap <- frequencyMaps$divide(countMaps)
  # add count map and frequency map to single raster
  # waterMaps <- waterMaps$addBands(countMaps)
  # Rename count map to extent map and cast as float
  # extentMap <- waterMaps$select('Count_Map')$rename('Extent_Map')$toFloat() #??????????????????????????????????????
  # add extent map to water maps raster
  # waterMaps <- ee$Image(waterMaps)$addBands(extentMap)
  # Map$addLayer(waterMaps$select('Count_Map'))
  
  ###
  varName <- paste('water', path, sep = '')
  # assign(varName, waterMaps$updateMask(waterMaps$gt(0)))
  assign(varName, waterMap$updateMask(waterMap$gt(0)))
  
}

caucasusMosaic <- ee$ImageCollection(c(water166, water167, water168, water169, 
                                       water170, water171, water172, water173))$
  # select('Extent_Map')$
  mosaic()$
  unmask()$
  clip(caucasusAdmin)

Map$addLayer(caucasusMosaic)

rectangle <-ee$Geometry$Rectangle(40, 38, 50.5, 43.7)

## drive - Method 01
# img_02 <- ee_as_raster(
#   image = img,
#   region = rectangle,
#   scale = 30,
#   crs= "EPSG:32638",
#   via = "drive"
# )

# Map$addLayer(waterMaps166)

taskImage <- ee_image_to_drive(
  image = caucasusMosaic,
  description = "waterMap",
  # folder = 'rgee_raster',
  region = rectangle,
  scale = 30,
  crs = 'EPSG:32638',
  maxPixels = 10000000000000
)
taskImage$start()
ee_monitoring(taskImage)

# taskImage <- ee_image_to_asset(
#   image = CaucasusMosaic$select('Extent_Map'),
#   assetId = assetPath,
#   crs = 'EPSG:32638',
#   scale = 30,
#   maxPixels = 10000000000000
# )
# taskImage$start()
# ee_monitoring(taskImage)


# 
# ###############################################################################################################
# ###############################################################################################################
# 
# path167 <- 167
# rows167 <- c(031, 032, 033)
# 
# ###############################################################################################################
# ###############################################################################################################
# 
# path168 <- 168
# rows168 <- c(031, 032, 033)
# 
# ###############################################################################################################
# ###############################################################################################################
# 
# path169 <- 169
# rows169 <- c(031, 032, 033)
# dateEnd <- '2019-09-05'
# 
# ###############################################################################################################
# ###############################################################################################################
# 
# path170 <- 170
# rows170 <- c(030, 031, 032)
# dateEnd <- '2019-10-01'
# 
# ###############################################################################################################
# ###############################################################################################################
# 
# path171 <- 171
# rows171 <- c(030, 031)
# 
# ###############################################################################################################
# ###############################################################################################################
# 
# path172 <- 172
# rows172 <- c(030, 031)
# 
# ###############################################################################################################
# ###############################################################################################################
# 
# path173 <- 173
# rows173 <- c(029,030)
# dateEnd <- '2019-10-03'
