

stat_27<-read.table("Station27_Monthly_Means.odf", header=TRUE)

yr<-data.frame(seq(from=1977, to=2018, by=1))
colnames(yr)<-"yr"
stat_27$month<-substr(stat_27$INFO,4,6)
stat_27$yr<-substr(stat_27$INFO, 8,11)

stat_27<-stat_27[!stat_27$TEMP == -99, ]

depths<-unique(stat_27$DEPTH)
mean_temp<-matrix(nrow=42, ncol=length(depths))
jun_temp<-matrix(nrow=42, ncol=length(depths))
may_temp<-matrix(nrow=42, ncol=length(depths))
apr_temp<-matrix(nrow=42, ncol=length(depths))
for(i in 1:length(depths)){
  shallow<-subset(stat_27, DEPTH==depths[i] & yr>=1977)
  shallow.apr<-subset(shallow, month=="APR")
  shallow.apr2<-merge(yr, shallow.apr, by="yr", all=TRUE)
  shallow.may<-subset(shallow, month=="MAY")
  shallow.may2<-merge(yr, shallow.may, by="yr", all=TRUE)
  shallow.jun<-subset(shallow, month=="JUN")
  shallow.jun2<-merge(yr, shallow.jun, by="yr", all=TRUE)
  #mean_temp[,i]<-(shallow.may2$TEMP+shallow.jun2$TEMP)/2
  apr_temp[,i]<-shallow.may2$TEMP-shallow.apr2$TEMP
  may_temp[,i]<-shallow.may2$TEMP
  jun_temp[,i]<-shallow.may2$TEMP-shallow.jun2$TEMP
}

library(zoo)

#use rolling means from the surrounding 6 observations

apr_temp2<-apr_temp
for(i in 1:10){
  for(j in 1:42){
    if(is.na(apr_temp2[j,i]&j<=39)){apr_temp2[j,i]<-mean(c(apr_temp2[j-3,i], apr_temp2[j-2,i], apr_temp2[j-1,i],apr_temp2[j+1,i], apr_temp2[j+2,i], apr_temp2[j+3,i]), na.rm=TRUE)}
    if(is.na(apr_temp2[j,i]&j==40)){apr_temp2[j,i]<-mean(c(apr_temp2[j-3,i], apr_temp2[j-2,i], apr_temp2[j-1,i],apr_temp2[j+1,i], apr_temp2[j+2,i]), na.rm=TRUE)}
    if(is.na(apr_temp2[j,i]&j==41)){apr_temp2[j,i]<-mean(c(apr_temp2[j-3,i], apr_temp2[j-2,i], apr_temp2[j-1,i],apr_temp2[j+1,i]), na.rm=TRUE)}
    if(is.na(apr_temp2[j,i]&j==42)){apr_temp2[j,i]<-mean(c(apr_temp2[j-3,i], apr_temp2[j-2,i], apr_temp2[j-1,i]), na.rm=TRUE)}
  }
}  

jun_temp2<-jun_temp
for(i in 1:10){
  for(j in 1:42){
    if(is.na(jun_temp2[j,i]&j<=39)){jun_temp2[j,i]<-mean(c(jun_temp2[j-3,i], jun_temp2[j-2,i], jun_temp2[j-1,i],jun_temp2[j+1,i], jun_temp2[j+2,i], jun_temp2[j+3,i]), na.rm=TRUE)}
    if(is.na(jun_temp2[j,i]&j==40)){jun_temp2[j,i]<-mean(c(jun_temp2[j-3,i], jun_temp2[j-2,i], jun_temp2[j-1,i],jun_temp2[j+1,i], jun_temp2[j+2,i]), na.rm=TRUE)}
    if(is.na(jun_temp2[j,i]&j==41)){jun_temp2[j,i]<-mean(c(jun_temp2[j-3,i], jun_temp2[j-2,i], jun_temp2[j-1,i],jun_temp2[j+1,i]), na.rm=TRUE)}
    if(is.na(jun_temp2[j,i]&j==42)){jun_temp2[j,i]<-mean(c(jun_temp2[j-3,i], jun_temp2[j-2,i], jun_temp2[j-1,i]), na.rm=TRUE)}
  }
}  


for(i in 1:10){
  for(j in 1:42){
    apr_temp[j,i]<-replace(apr_temp[j,i],is.na(apr_temp[j,i]),apr_temp2[j,i])
    jun_temp[j,i]<-replace(jun_temp[j,i],is.na(jun_temp[j,i]),jun_temp2[j,i])
  }
}

subset.dat<-read.csv("subset_survey_dat3.csv")
catch.dat<-read.csv("common_spp_dat3.csv")

# shallow$spec.map<-paste(shallow$spec,shallow$map, sep=".")
# spp.freq2<-data.frame(table(shallow$spec.map))
# spp.freq3<-data.frame(table(spp.freq2))
# ll<-subset(spp.freq2, Freq>1)


#install.packages("Rcpp")
library(RGeostats)
#library(devtools)
library(rgeos)
library(rgdal)
library(raster)
library(plyr)



#some tows caught none of the common spp in any year
#sanity check to see that they are still there, consider removing later?
length(unique(subset.dat$map))
length(unique(catch.dat$map))

#get a vector of the most common spp
most.common<-as.character(unique(catch.dat$name))

#only 5s
only.5s<-subset(subset.dat, rec==5)

april_sub<-subset(only.5s, month==4)
april_sub$timing<-april_sub$month+april_sub$day/30

may_sub<-subset(only.5s, month==5)
may_sub$timing<-may_sub$month+may_sub$day/31

june_sub<-subset(only.5s, month==6)
june_sub$timing<-june_sub$month+june_sub$day/30

new_days<-rbind(april_sub, may_sub, june_sub)
avg_month<-aggregate(new_days$timing~new_days$year, FUN=mean)
colnames(avg_month)<-c("yr","timing")

useful_data<-data.frame(only.5s$lat.start, only.5s$long.start, only.5s$year, only.5s$month, only.5s$set.depth.mean, only.5s$bot.temp)
colnames(useful_data)<-c("lat","lon","yr","month","depth","temp")
depths_code<-c(depths, Inf)
useful_data$depthmap<-cut(useful_data$depth, depths_code, labels=c(1,2,3,4,5,6,7,8,9,10))
useful_data<-merge(useful_data, avg_month, by="yr")

unique_yrs<-unique(useful_data$yr)

apr_temp<-data.frame(apr_temp)
apr_temp$yr<-yr$yr
apr_temp<-merge(apr_temp, avg_month, by="yr")
jun_temp<-data.frame(jun_temp)
jun_temp$yr<-yr$yr
jun_temp<-merge(jun_temp, avg_month, by="yr")

par(mfrow=c(2,1))
plot(apr_temp$X1~apr_temp$yr, type="l", lwd=2, ylim=c(-0.5,4), xlab="Year", ylab="Temp Adjustment",cex.lab=1.5, cex.axis=1.25, main="April")
lines(apr_temp$X10~apr_temp$yr, lwd=2, col="blue")
#legend("topright", legend=c("Surface", "Bottom"), lwd=2, col=c("black","blue"), bty="n")

plot(jun_temp$X1~jun_temp$yr, type="l", lwd=2, ylim=c(-6,1), xlab="Year", ylab="Temp Adjustment",cex.lab=1.5, cex.axis=1.25, main="June")
lines(jun_temp$X10~jun_temp$yr, lwd=2, col="blue")
legend(1998,-3.75, legend=c("Surface", "Bottom"), lwd=2, col=c("black","blue"), bty="n")




depths<-c(1,2,3,4,5,6,7,8,9,10)

depth_sub<-matrix(nrow=1,ncol=9)
colnames(depth_sub)<-c(colnames(useful_data),"temp_adj")
for(d in 4:10){
  depth_calc<-subset(useful_data, depthmap==depths[d])
  if (nrow(depth_calc) == 0) next
  for(y in 1:length(unique_yrs)){
    year_calc<-subset(depth_calc, yr==unique_yrs[y])
    if (nrow(year_calc) == 0) next
    if(year_calc$timing<=5.25){year_calc$temp_adj<-apr_temp[y,d+1]}
    if(year_calc$timing>=5.75){year_calc$temp_adj<-jun_temp[y,d+1]}
    if(year_calc$timing<5.75 & year_calc$timing>5.25){year_calc$temp_adj<-0}
    depth_sub<-rbind(depth_sub, year_calc)
  }
}

depth_sub<-depth_sub[2:10899,]
depth_sub$temp_adj[is.na(depth_sub$temp_adj)] <- 0
depth_sub$mod_temp<-depth_sub$temp+depth_sub$temp_adj



library(tmap)
library(grid)

library(RGeostats)
#library(marmap)
#library(maptools)
#library(oce)
#library(ocedata)
library(rgeos)
require(rgdal)
library(raster)
library(RColorBrewer)

#bathymetry data
deep.w <- readOGR(dsn = "ne_10m_bathymetry_J_1000", layer = "ne_10m_bathymetry_J_1000")
land <- readOGR(dsn = "ne_10m_bathymetry_L_0", layer = "ne_10m_bathymetry_L_0")


x<-c(-59,-43)
y<-c(40,52)
xy<-cbind(x,y)
S<-SpatialPoints(xy)
bounds<-bbox(S)
crop.land<-crop(land, S)
#plot(crop.land, col="red")

crop.deep<-crop(deep.w, S)
#plot(crop.deep, col="red")

union_sp<-gDifference(crop.land, crop.deep)


#open spatial list
load("spatial_list.Rdata")
poi<-Spatial_List$loc_g

sp_poi<-SpatialPoints(poi, proj4string = CRS("+proj=utm +zone=22 +datum=WGS84 +units=km +no_defs")) 

sp_poidf<-SpatialPointsDataFrame(poi, data.frame(seq(from=1, to=50, by=1)))

#NOT SURE IF I NEED THESE SPECIFIC VALS
#r <- raster(ncols=100, nrows=100, xmn=229.8978,xmx=807.9266,ymn= 4745.421, ymx=5437.54, crs= CRS("+proj=utm +zone=22 +units=km"))
r <- raster(ncols=100, nrows=100, xmn=229.8978,xmx=810,ymn= 4500, ymx=5437.54, crs= CRS("+proj=utm +zone=22 +units=km"))

newthing <- as(r,"SpatialPoints")

union_sp_new<-spTransform(union_sp, CRS("+proj=utm +zone=22 +units=km"))

cropped_sp<-crop(newthing, union_sp_new)

spix<-SpatialPixelsDataFrame(cropped_sp, data.frame(temp=rep(NA,length(cropped_sp))))

proj4string(spix) = CRS(proj4string(new_spdf))


depth_sub$yr<-as.numeric(depth_sub$yr)

yr<-unique(depth_sub$yr)


new_yr<-yr[i]
one_yr<-subset(depth_sub, yr==new_yr)
s<-cbind(-one_yr$lon, one_yr$lat)

one_yr<-na.omit(one_yr, cols="temp")

sp_coords<-SpatialPoints(data.frame(-one_yr$lon, one_yr$lat), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
spdf<-SpatialPointsDataFrame(sp_coords, data.frame(temp=one_yr$temp))


new_spdf<-spTransform(spdf, CRS("+proj=utm +zone=22 +units=km"))


# PbKrig <- autoKrige(temp~1,new_spdf, spix)
# #plot(PbKrig)
# 
# prediction_spdf = PbKrig$krige_output
# 
# jpeg(paste(new_yr,"temp.jpeg", sep="_"),  width=1000, height=1000)
# tm_shape(prediction_spdf)+
#   tm_raster("var1.pred", breaks=c(-3,-2,-1,0,1,2,3,4,5,6,7,8), midpoint=0, title="Temperature (C)",palette="-RdYlBu")+
#   tm_layout( title=new_yr, title.position = c("RIGHT", "BOTTOM"), title.size=2)+
#   tm_legend(legend.outside=TRUE)
# dev.off()

library(automap)

yr<-unique(depth_sub$yr)

temp_mat<-matrix(nrow=50, ncol=35)
for( i in 1:35){
new_yr<-yr[i]
one_yr<-subset(depth_sub, yr==new_yr)
s<-cbind(-one_yr$lon, one_yr$lat)

one_yr<-na.omit(one_yr, cols="temp")

sp_coords<-SpatialPoints(data.frame(-one_yr$lon, one_yr$lat), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
spdf<-SpatialPointsDataFrame(sp_coords, data.frame(temp=one_yr$temp))


new_spdf<-spTransform(spdf, CRS("+proj=utm +zone=22 +units=km"))


PbKrig <- autoKrige(temp~1,new_spdf, spix)
#plot(PbKrig)

prediction_spdf = PbKrig$krige_output


new_rast <- raster(prediction_spdf)

temp_mat[,i]<-extract(new_rast, sp_poi,df=TRUE)[,2] 
}


new_mat2<-matrix(nrow=35, ncol=35)
for(k in 1:35){
  for(i in 1:35){
    new_mat2[i,k]<-cor(temp_mat[,i], temp_mat[,k], method="pearson", use="complete.obs")
    #new_mat2[i,k]<-MSE(temp_mat[,i], temp_mat[,k])
  }
}

new_mat3<-new_mat2
for(i in 1:35){
  new_mat3[i,i:35]<-NA
}

axismath<-function(num){
  ans<-0+(((1/34))*num)
  return(ans)}
jpeg("heatmaps_spearman_temps_no_buff.jpeg", height=1000, width=1000)
par(mar=c(6,3,2,6))
fields::image.plot(new_mat3, col = heat.colors(12), axes=F, horizontal = TRUE)#, breaks=c(-1,-.75,-.5,-.25,0,.25,.5,.75,1))
axis(1, at=c(axismath(1), axismath(10), axismath(18), axismath(25), axismath(32)), labels=c(1978, 1990, 1998,
                                                                                            2005, 2014), cex.axis=2)
axis(4, at=c(axismath(1), axismath(10), axismath(18), axismath(25), axismath(32)), labels=c(1978, 1990, 1998,
                                                                                            2005, 2014), las=1, cex.axis=2)
dev.off()





temp_df<-data.frame(temp=apply(temp_mat, 2, mean, na.rm=TRUE), sd=apply(temp_mat, 2, sd, na.rm=TRUE))


temp_df$yr<-yr

all_yrs<-data.frame("yr"=seq(from=1977, to=2018, by=1))

merged_temp<-merge(temp_df, all_yrs, by="yr", all=TRUE)

png("temp_trend.png", width=1000, height=600)
plot(merged_temp$temp~merged_temp$yr, type="l", lwd=5, bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
dev.off()

lines(merged_temp$temp+merged_temp$sd~merged_temp$yr)

 #convert into more useable objects




prediction_spdf2 <- as(prediction_spdf,"SpatialPointsDataFrame")

tm_shape(prediction_spdf2)+
  tm_dots("var1.pred")
