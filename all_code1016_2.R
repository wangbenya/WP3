## load the library 
rm(list=ls())


library(raster)
library(sp)
library(randomForest)
library(magrittr)
library(rgdal)
library(gstat)
library(ggplot2)
library(mlr)
library(SemiPar)
library(Hmisc)
library(foreign)
library(maptools)
library(prettymapr)
library(mlrMBO)
library(parallelMap)
library(caret)
library(automap)
library(reshape2)

## start the parallel 
parallelStartSocket(16)

WGS84 <- CRS("+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

study_area <- shapefile("~/WP2/data/study_area.shp")
water <- shapefile("~/WP2/data/water.shp")

## load the veg, soil, land use, groundwater subarea,
## surface water subare, catchment
soil <- raster("~/WP2/data/soil1.ovr")
veg <- raster("~/WP2/data/vegetation.ovr")
land_use <- raster("~/WP2/data/landuse1.ovr")
ss <- raster("~/WP2/data/Ssubarea.ovr")
gs <- raster("~/WP2/data/Gsubarea.ovr")
cat <- raster("~/WP2/data/catch_name2.tif.ovr")

## define the function 
## preprocess 
study_area <- spTransform(study_area, WGS84)
extent <- c(study_area@bbox[1, 1:2], study_area@bbox[2, 1:2])

water <- spTransform(water, WGS84)
#plot(water)

study_area_withW <- symdif(study_area, water)
#plot(study_area_withW)

pre <- function(x) {
    projection(x) <- WGS84
    extent(x) <- extent
    x <- raster::mask(x, study_area)
    return(x)
}

read_points <- function(read_data) {
    SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
    SP <- spTransform(SP, WGS84)
    SP@bbox <- study_area@bbox
    if (length(zerodist(SP)) >= 1) {
        SP <- SP[-(zerodist(SP)[, 1]),]
    }
 #   plot(study_area_withW)
  #  points(SP@coords)
    return(SP)
}

read_pointDataframes <- function(read_data) {
    SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
    SPD <- SpatialPointsDataFrame(SP, read_data)
    SPD <- spTransform(SPD, WGS84)
    SPD@bbox <- study_area@bbox
    if (length(zerodist(SPD)) >= 1) {
        SPD <- SPD[-(zerodist(SPD)[, 1]),]
    }
   # plot(study_area_withW)
    #points(SPD@coords)
    return(SPD)
}


reclass <- function(df, i, j) {
    df[, "DON"][df[, "DON"] <= i] <- "Low"
    df[, "DON"][df[, "DON"] < j] <- "Medium"
    df[, "DON"][(df[, "DON"] != "Low") & (df[, "DON"] != "Medium")] <- "High"
    df[, "DON"] <- factor(df[, "DON"], levels = c("Low", "Medium", "High"))
    return(df)
}


reclass2 <- function(df) {
    df[, "DON"][df[, "DON"] == 1] <- "Low"
    df[, "DON"][df[, "DON"] == 2] <- "Medium"
    df[, "DON"][(df[, "DON"] != "Low") & (df[, "DON"] != "Medium")] <- "High"
    df[, "DON"] <- factor(df[, "DON"], levels = c("Low", "Medium", "High"))
    return(df)
}

reclass3 <- function(df, i, j) {
  df[, "DON_k"][df[, "DON_k"] <= i] <- "Low"
  df[, "DON_k"][df[, "DON_k"] < j] <- "Medium"
  df[, "DON_k"][(df[, "DON_k"] != "Low") & (df[, "DON_k"] != "Medium")] <- "High"
  df[, "DON_k"] <- factor(df[, "DON_k"], levels = c("Low", "Medium", "High"))
  return(df)
}

# Add X and Y to training 
add_S1S2 <- function(dataset) {
    dataset$s1 <- coordinates(dataset)[, 1]
    dataset$s2 <- coordinates(dataset)[, 2]
    return(dataset)
}

## preprocess the landscape raster
soil <- pre(soil)
veg <- pre(veg)
land_use <- pre(land_use)
ss <- pre(ss)
gs <- pre(gs)
cat <- pre(cat)

v_veg<-values(veg)
v_veg[v_veg %in% c(2,3,4)]=1
v_veg[v_veg %in% c(8,9)]=8
v_veg[v_veg %in% c(12,13)]=12
v_veg[v_veg %in% c(18,19,20)]=18
values(veg)<-v_veg

v_land<-values(land_use)
v_land[v_land %in% c(1,2,13)]=1
v_land[v_land %in% c(5,7,12,6,11)]=5
v_land[v_land %in% c(8,10)]=8
values(land_use)<-v_land

# Create an empty grid where n is the total number of cells
r <- raster(study_area)
res(r) <- res(soil) # 10 km if your CRS's units are in km
base_grid <- as(r, 'SpatialGrid')
#plot(base_grid)

## M2, using RF to predict the DON
depth <- read.csv("~/WP2/data/sampling_depth.csv",header=T) %>% read_pointDataframes(.)

# Define the 1st order polynomial equation
f_depth <- as.formula(sampling_d ~ s1 + s2)
# Add X and Y to training 
depth<-add_S1S2(depth)
# Run the regression model
lm.depth <- lm(f_depth, data = depth)
# variogram on the de-trended data.
var.depth <- variogram(f_depth, depth)
#plot(var.depth)
dat.fit_depth <- fit.variogram(var.depth,vgm(c("Sph","Exp","Gau","Lin","Spl")))
#plot(var.depth, dat.fit_depth, xlim = c(0, 40000))

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
depth_k <- krige(f_depth, depth, base_grid, dat.fit_depth) %>% raster(.) %>% raster::mask(., study_area)
#plot(depth_k)
depth_k@data@names<-"GW_depth"

#convert the raster to points for plotting
map_depth <- rasterToPoints(depth_k)
#Make the points a dataframe for ggplot
df <- data.frame(map_depth)
#Make appropriate column headings
colnames(df) <- c("Longitude", "Latitude", "GW_depth")

#Now make the map
### distance
water <- raster::rasterize(water, depth_k)
water_distance <- raster::mask(distance(water),study_area)
water_distance@data@names<-"Distance_to_water"
#plot(water_distance)
## load the data 
landscapes<-stack(soil,veg,land_use,ss,gs,cat,depth_k,water_distance)
names(landscapes) <- c("Soil", "Veg", "Landuse", "SS", "GS", "Catchment", "GW_depth", "Distance")

## set the parameters for mlr
seed=35
set.seed(seed)
reg_rf = makeLearner("regr.randomForest")
#reg_rf$par.vals<-list(importance=T)

class_rf = makeLearner("classif.randomForest")
#class_rf$par.vals<-list(importance=T)

rdesc = makeResampleDesc("CV", iters = 3)

## define the parameters search mehtods
ctrl = makeTuneControlIrace(maxExperiments = 400L)

## define the parameter spaces for RF
para_rf = makeParamSet(
  makeIntegerParam("ntree", lower = 300, upper = 600),
  makeIntegerParam("nodesize", lower = 3, upper = 12),
  makeIntegerParam("mtry", lower = 2, upper = 8)
)

model_build <- function(dataset, n_target, method) {
  set.seed(35)
  if (method == "reg") {
    ## define the regression task for DON 
    WP3_target = makeRegrTask(id = "WP3_target", data = dataset, target = n_target)
    ## cross validation
    ## 10-fold cross-validation
    rin = makeResampleInstance(rdesc, task = WP3_target)
    ## tune the parameters for rf and xgboost
    res_rf = mlr::tuneParams(reg_rf, WP3_target, resampling = rdesc, par.set = para_rf, control = ctrl,
                             show.info = FALSE, measures = rsq)
    
    ## set the hyperparameter for rf and xgboost 
    lrn_rf = setHyperPars(reg_rf, par.vals = res_rf$x)
    
  } else {
    ## define the regression task for DON 
    WP3_target = makeClassifTask(id = "WP3_target", data = dataset, target = n_target)
    ## cross validation
    ## 10-fold cross-validation
    rin = makeResampleInstance(rdesc, task = WP3_target)
    res_rf = mlr::tuneParams(class_rf, WP3_target, resampling = rdesc, par.set = para_rf, control = ctrl,
                             show.info = FALSE, measures = acc)
    lrn_rf = setHyperPars(class_rf, par.vals = res_rf$x)
  }
  
  ## train the final model 
  set.seed(35)
  rf <- mlr::train(lrn_rf, WP3_target)
  return(rf)
}

## load the data 
all_results<-data.frame()
set.seed(66)
seed.list<-sample(1:1000,50,replace =F)

all_points<-read.csv("~/WP2/data/all_data1210.csv",header = T)
extra_n<-read.csv("~/WP2/data/extra_n.csv",header = T)
extra_n<-subset(extra_n,!(extra_n$WIN_Site_ID %in% all_points$WIN_Site_ID))

for (tt in seq(1,5)){
  
print(tt)
seeds<-seed.list[tt]
set.seed(seeds)
trainIndex <- createDataPartition(all_points$DON, p = .8, list = FALSE)

training <- all_points[trainIndex,]
testing <- all_points[-trainIndex,]

## load the point data 
training_df <- training[,c(1,2,3,5)] %>% rbind(.,extra_n[,c(1,2,3,5)]) %>%
  subset(.,.[,"DON"]!="NA") %>% read_pointDataframes(.) 

testing_df <- testing[,c(1,2,3,5)] %>% read_pointDataframes(.) 

training_points<-training[,c(1,2,3,5)] %>% rbind(.,extra_n[,c(1,2,3,5)]) %>%
  subset(.,.[,"DON"]!="NA") %>% read_points(.) 
testing_points <- read_points(testing)

## map1, using kringing for DON interpolation
f.1 <- as.formula(log10(DON) ~ s1 + s2)

# Add X and Y to training 
training_df<-add_S1S2(training_df)
testing_df<-add_S1S2(testing_df)

# Compute the sample variogram; note that the f.1 trend model is one of the
var.smpl1 <- variogram(f.1, training_df)
# Compute the variogram model by passing the nugget, sill and range value
dat.fit1 <- fit.variogram(var.smpl1, vgm(c("Exp","Sph","Gau","Lin","Spl")))
# Perform the krige interpolation (note the use of the variogram model
map1 <- krige(f.1, training_df, base_grid, dat.fit1) %>% raster(.) %>% raster::mask(., study_area)
values(map1) <- 10 ^ (values(map1))
dat.krg_DON<-map1


map1_predict <- data.frame(observed_DON=testing_df@data$DON,predicted_DON=raster::extract(map1, testing_points))

for (t in c(1,2)){
  map1_predict[, t][map1_predict[, t] <=0.6] <- "Low"
  map1_predict[, t][map1_predict[, t] < 1.2] <- "Medium"
  map1_predict[, t][(map1_predict[, t] != "Low") & (map1_predict[, t] != "Medium")] <- "High"
  map1_predict[, t] <- factor(map1_predict[, t], levels = c("Low", "Medium", "High"))
  
}

confusionMatrix(map1_predict[,2],map1_predict[,1])$overall


## M2, using RF to predict the DON
landscape_train <- raster::extract(landscapes, training_points)
landscape_test <- raster::extract(landscapes, testing_points)

M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON","Longitude","Latitude")])
M2_test <- cbind(as.data.frame(landscape_test), testing_df@data[c("DON","Longitude","Latitude")])

names(M2_train) <- colnames(M2_test)

## preprocess for catogoried variables 
transCat<-function(i){
  var_name<-names(M2_train)[i]
  a<-as.list(unique(M2_train[,var_name]))
  b<- data.frame()
  
  for (x in seq(1, length(a))) {
    value = a[[x]]
    data1 <- subset(M2_train, M2_train[,var_name] == value)
    c <- data.frame(value, mean(data1$DON), length(data1$DON))
    b <- rbind(b,c)
  }
  
  colnames(b) <- c(var_name, paste0("mean_DON_",var_name),"NO.")
  return(b)
}

b1<-transCat(1)
b2<-transCat(2)
b3<-transCat(3)
b4<-transCat(4)
b5<-transCat(5)
b6<-transCat(6)

soil_max = b1[which(b1$NO. == max(b1$NO.)),][, 2]
veg_max = b2[which(b2$NO. == max(b2$NO.)),][, 2]
landuse_max = b3[which(b3$NO. == max(b3$NO.)),][, 2]
ss_max = b4[which(b4$NO. == max(b4$NO.)),][, 2]
GS_max = b5[which(b5$NO. == max(b5$NO.)),][, 2]
cat_max = b6[which(b6$NO. == max(b6$NO.)),][, 2]

base1 <- (merge(b1, M2_train, by = 'Soil', all.y = T))
base2 <- (merge(b2, base1, by = 'Veg', all.y = T))
base2 <- base2[, - c(3, 6)]

test1 <- (merge(b1, M2_test, by = 'Soil', all.y = T))
test1[is.na(test1)] <- soil_max
test1 <- test1[, - c(3)]

test2 <- (merge(b2, test1, by = 'Veg', all.y = T))
test2[is.na(test2)] <- veg_max
test2 <- test2[, - c(3)]

base3 <- (merge(b3, base2, by = 'Landuse', all.y = T))
base4 <- (merge(b4, base3, by = 'SS', all.y = T))
base4 <- base4[, - c(3, 6)]

test3 <- (merge(b3, test2, by = 'Landuse', all.y = T))
test3[is.na(test3)] <- landuse_max
test3 <- test3[, - c(3)]

test4 <- (merge(b4, test3, by = 'SS', all.y = T))
test4[is.na(test4)] <- ss_max
test4 <- test4[, - c(3)]

base5 <- (merge(b5, base4, by = 'GS', all.y = T))
base6 <- (merge(b6, base5, by = 'Catchment', all.y = T))
base6 <- base6[, - c(3, 6)]

test5 <- (merge(b5, test4, by = 'GS', all.y = T))
test5[is.na(test5)] <- GS_max
test5 <- test5[, - c(3)]

test6 <- (merge(b6, test5, by = 'Catchment', all.y = T))
test6[is.na(test6)] <- cat_max
test6 <- test6[, - c(3)]

## create the training and testing sets 
WP2Train <- base6[, c(12, 10, 8, 6, 4,2, 13:17)]
WP2Test <- test6[, c(12, 10, 8, 6, 4,2, 13:17)]

WP2Train<-reclass(WP2Train,0.6,1.2)
WP2Test<-reclass(WP2Test,0.6,1.2)

## build the model for map2
names(WP2Train)<-c("Soil", "Veg", "Landuse", "SS", "GS", "Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
names(WP2Test)<-c("Soil", "Veg", "Landuse", "SS", "GS", "Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")

set.seed(seeds)
rf_DON_m2 <- model_build(WP2Train, "DON","clf")

map2_predict <- predict(rf_DON_m2, newdata = WP2Test)
confusionMatrix(map2_predict$data$response, map2_predict$data$truth)$overall

## map4, kriging first and then rf
# kriging for DOC
f.DOC <- as.formula(log10(DOC) ~ s1 + s2)
training_DOC <- training[,c(1,2,3,4)] %>% rbind(.,extra_n[,c(1,2,3,4)]) %>%
  subset(.,.[,"DOC"]!="NA") %>% read_pointDataframes(.) 
training_DOC<-add_S1S2(training_DOC)
var.smpl_DOC <- variogram(f.DOC, training_DOC)
#plot(var.smpl_DOC)

dat.fit_DOC <- fit.variogram(var.smpl_DOC,vgm(c("Exp","Sph","Gau","Lin","Spl")))
#plot(var.smpl_DOC,dat.fit_DOC)
# Perform the krige interpolation (note the use of the variogram model
dat.krg_DOC <- krige(f.DOC, training_DOC, base_grid, dat.fit_DOC) %>% raster(.) %>% raster::mask(., study_area)
values(dat.krg_DOC) <- 10 ^ (values(dat.krg_DOC))

# kriging for NH4
f.NH4 <- as.formula(log10(NH4) ~ s1 + s2)

training_NH4 <- training[,c(1,2,3,8)] %>% rbind(.,extra_n[,c(1,2,3,7)]) %>%
  subset(.,.[,"NH4"]!="NA") %>% read_pointDataframes(.) 
training_NH4<-add_S1S2(training_NH4)

var.smpl_NH4 <- variogram(f.NH4, training_NH4)

#plot(var.smpl_NH4)
dat.fit_NH4 <- fit.variogram(var.smpl_NH4, vgm(c("Exp","Sph","Gau","Lin","Spl")))

# Perform the krige interpolation (note the use of the variogram model
dat.krg_NH4 <- krige(f.NH4, training_NH4, base_grid, dat.fit_NH4) %>% raster(.) %>% raster::mask(., study_area)
values(dat.krg_NH4) <- 10 ^ (values(dat.krg_NH4))

# kriging for NOx
f.NOx <- as.formula(log10(NOx) ~ s1 + s2)

training_NOx <- training[,c(1,2,3,7)] %>% rbind(.,extra_n[,c(1,2,3,6)]) %>%
  subset(.,.[,"NOx"]!="NA") %>% read_pointDataframes(.) 
training_NOx<-add_S1S2(training_NOx)

var.smpl_NOx <- variogram(f.NOx, training_NOx)
#plot(var.smpl_NOx)

dat.fit_NOx <- fit.variogram(var.smpl_NOx,vgm(c("Sph","Exp","Gau","Lin","Spl")))
# Perform the krige interpolation 
dat.krg_NOx <- krige(f.NOx, training_NOx, base_grid, dat.fit_NOx) %>% raster(.) %>% raster::mask(., study_area)
values(dat.krg_NOx) <- 10 ^ (values(dat.krg_NOx))

## create rasterstack with kriging data
kriging_nutrietn<-stack(dat.krg_DON,dat.krg_DOC, dat.krg_NH4, dat.krg_NOx)
names(kriging_nutrietn) <- c("DON_k", "DOC_k", "NH4_k", "NOx_k")

## extract the data from landscapes_withN
landscape_train_withKN <- raster::extract(kriging_nutrietn, read_points(base6[,15:17]))
landscape_test_withKN <- raster::extract(kriging_nutrietn, read_points(test6[,15:17]))

m22<-predict(rf_DON_m2, newdata = WP2Train)

M4_train_withKN <- cbind(base6[, c(12, 10, 8, 6, 4,2, 13:17)],as.data.frame(landscape_train_withKN))
M4_test_withKN <- cbind(test6[, c(12, 10, 8, 6, 4,2, 13:17)],as.data.frame(landscape_test_withKN))


names(M4_test_withKN) <- names(M4_train_withKN)

M4_train_withKN <- reclass(M4_train_withKN,0.6,1.2)
M4_test_withKN <- reclass(M4_test_withKN,0.6,1.2)

#M4_train_withKN <- reclass3(M4_train_withKN,0.6,1.2)
#M4_test_withKN <- reclass3(M4_test_withKN,0.6,1.2)

## 
set.seed(seeds)
DON_rf_m4<-model_build(M4_train_withKN,"DON","clf")

## map3 predict accuracy
map4_predict<-predict(DON_rf_m4,newdata=M4_test_withKN)

confusionMatrix(map4_predict$data$response,map4_predict$data$truth)$overall

predict_results<-data.frame(seeds,map1_predict,map2_predict$data,map4_predict$data)

all_results<-rbind(all_results,predict_results)

}

write.csv(all_results,file="~/WP2/results/all_results1013.csv",row.names = F)

seedss<-unique(all_results$seeds)

all_acc<-data.frame()

for (qq in seedss){
  sub_data<-subset(all_results,all_results$seeds==qq)
  print(dim(sub_data))
  acc_1<-postResample(sub_data[,3],sub_data[,2])[1]
  acc_2<-postResample(sub_data[,5],sub_data[,4])[1]
  acc_3<-postResample(sub_data[,7],sub_data[,6])[1]
  single_acc<-data.frame(qq,acc_1,acc_2,acc_3)
  all_acc<-rbind(all_acc,single_acc)
  
}

print(all_acc)
#all_acc<-melt(all_acc,id="qq")

#ggplot(data=all_acc,aes(x=all_acc$variable,y=all_acc$value))+geom_boxplot()







