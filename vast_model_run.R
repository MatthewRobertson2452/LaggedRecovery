
# Install INLA using currently recommended method
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(devtools)
# Install package
#install_github("james-thorson/VAST")
# Load package
library(VAST)
library(TMB)
library(INLA)


#install.packages("sendmailR")

subset.dat<-read.csv("conv_subset_survey_dat4.csv")
catch.dat<-read.csv("conv_spp_dat2.csv")
#get a vector of the most common spp
most.common<-as.character(unique(catch.dat$name))

#only 5s
only.5s<-subset(subset.dat, rec==5 & year>=1985)

####For some reason this kept data for 1978, remove that next time

catch.dat<-subset(catch.dat, year>=1985 )
#map_yr<-data.frame(only.5s$map, only.5s$lat.start, only.5s$long.start, only.5s$year)
#colnames(map_yr)<-c("map", "lat.start", "long.start", "year")

#dataframe for years of interest
uni.yrs<-unique(catch.dat$year)

#choose the species to be examined
species<-c(most.common[1],most.common[2])


new_list<-list()
#subset catch data to only be for that species
for(i in 1:2){
  plaice.1975<-subset(catch.dat, name==species[i] & year>=1985 )
  
  #take out the map and the number caught
  plaice.1975<-plaice.1975[,c(35,36,49,50)]
  colnames(plaice.1975)<-c("num","wgt", "map","name")
  
  #merge catch data with set data to ensure sets with no catch are present
  spp1.1975<-plyr::join(plaice.1975, only.5s, by='map', type="full", match="first")
  
  #give any na's from the join a zero, actual tows with no catch of this species 
  spp1.1975$num[is.na(spp1.1975$num)] <- 0
  spp1.1975$wgt[is.na(spp1.1975$wgt)] <- 0
  
  spp1.1975$name[is.na(spp1.1975$name)]<-most.common[i]
  
  new_list[[i]]<-spp1.1975
}


full.1975<-rbind(new_list[[1]],new_list[[2]])

full.1975<-full.1975[!is.na(full.1975$lat.start),]
sum(is.na(full.1975$long.start))

full.1975$name<-factor(full.1975$name, levels=sort(unique(full.1975$name)))


#####PREP VAST

#this is really just for the session info file
Version = get_latest_version( package="VAST" )

Method = c("Grid", "Mesh", "Spherical_mesh")[2]#specify type of mesh
grid_size_km = 25 #specify grid size
n_x = 50  # Specify number of stations (a.k.a. "knots")

#should have had the same number of factors as spp
FieldConfig = c("Omega1"=2, "Epsilon1"=2, "Omega2"=2, "Epsilon2"=2) #omega is spatial var, epsilon is spatio-temporal var
#Omega1 is the spatial var in the encounter prob eq, Omega2 is the spatial var in the pos catch rate eq
#equal to a 1 indicates include, 0 indicates remove
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=4, "Epsilon2"=4) #specify correlation structure across time, see ?VAST::make_data
OverdispersionConfig = c("Eta1"=0, "Eta2"=0) #this is the "vessel effect", default is turned off (0), again see ?VAST::make_data
#can use AR1 instead of 0, eta1 is encounter prob, eta2 is pos catch rates
ObsModel = c(2,1)  #specify the link function and distribution for data, same help as noted aboved
#lots of options for this, first is for pos catch rates, second is for encounter prob

Options =  c("SD_site_density"=0, "SD_site_logdensity"=0, "Calculate_Range"=1, "Calculate_evenness"=1,
             "Calculate_effective_area"=1, "Calculate_Cov_SE"=1, 'Calculate_Synchrony'=1, 'Calculate_Coherence'=1)

strata.limits <- data.frame('STRATA'="All_areas")


#setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Dissertation Plans\\VAST\\New_converted_models\\Fixed Data\\100 knots")
DateFile = paste0(getwd())
#dir.create(DateFile)

Record = ThorsonUtilities::bundlelist( c("Version","Method","grid_size_km","n_x","FieldConfig","RhoConfig","OverdispersionConfig","ObsModel","Options") )
save( Record, file=file.path("converted_twospp_twofactors_50knots.RData"))
capture.output( Record, file="Record.txt")


full.1975$name<-factor(full.1975$name)
#data( EBS_pollock_data, package="FishStatsUtils" )
Data_Geostat = data.frame(  "spp"=full.1975$name,
                            "Catch_KG"=full.1975$wgt, "Year"=full.1975$year, 
                            "Vessel"=full.1975$data.series, "AreaSwept_km2"=full.1975$dist.towed,
                            "Lat"=full.1975$lat.start, "Lon"=-full.1975$long.start, "Pass"=0)

#Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits )


lat_lon_dat<-data.frame(Lat=full.1975$lat.start, Lon=-full.1975$long.start)

Extrapolation_List<-FishStatsUtils::make_extrapolation_info(observations_LL=lat_lon_dat,Region="Other")

Spatial_List = make_spatial_info( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon_i=Data_Geostat[,'Lon'], Lat_i=Data_Geostat[,'Lat'], 
                                  Extrapolation_List=Extrapolation_List, Save_Results=FALSE, fine_scale=FALSE)

save(Spatial_List, file=file.path("spatial_list.RData"))
save(Extrapolation_List, file=file.path("extrap_list.RData"))
# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

#double check the spatial domain first, if this is wrong then the whole model will be wrong
#plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat )


#include a "Q_ik" argument to add covars that affect catchability, see https://github.com/nwfsc-assess/geostatistical_delta-GLMM/wiki/Adding-new-covariates
#if the covar affects density it is put in 'X_xtp', however it needs to match the location of knots for every year, needs to be an array where
#X_xtp[x,t,p] is the value of the covariate for knot x at time t and covariate p (if you have only one covariate, then X_xtp is an array with 3rd dimension equal to one)
#to do this, Thorson says to use a nearest-neighbor algorithm, again see the link just above
TmbData = make_data("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, "RhoConfig"=RhoConfig,
                    "ObsModel"=ObsModel, "c_i"=as.numeric(Data_Geostat[,'spp'])-1, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'],
                    "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "spatial_list"=Spatial_List,
                    "Options"=Options )

#to speed up runs, I can use steps found here https://github.com/nwfsc-assess/geostatistical_delta-GLMM/wiki/Speed-up-simulations


#this is the actual model construction step...long step #1
TmbList = make_model("TmbData"=TmbData, "Version"=Version, "RhoConfig"=RhoConfig,
                     "loc_x"=Spatial_List$loc_x, "Method"=Spatial_List$Method)
Obj = TmbList[["Obj"]]

start.time <- Sys.time()
#the optimizer step, can also take a while, i.e. long step #2
Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=FALSE,
                           newtonsteps=1, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#check AIC
Opt$AIC
Opt$max_gradient


Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file="converted_50knots_2factors.RData")


Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]


#make sure no NAs
#if SigmaE, SigmaO, or SigmaVT sd's hit lower bound, these effects are basically being turned off
sdrep<-sdreport(Obj) #takes a while to run so only use when needed
save(Save, file="sdrep_new.RData")

#DateFile = paste0(getwd(),'/VAST_output_yellowtail/log_norm_vessel_AR1')

Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat)

Q = plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP="Posterior_Predictive",
                              FileName_Phist="Posterior_Predictive-Histogram", 
                              FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist")




# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="Other", "spatial_list"=Spatial_List, "Extrapolation_List"=Extrapolation_List )
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

#not sure why these don't save
TmbData$n_x<-n_x
TmbData$s_i<-Data_Geostat$knot_i-1

plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=TmbData, Report=Report, Q=Q,
               savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]],
               MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]],
               FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], 
               Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), 
               oma=c(3.5,3.5,0,0), cex=1.8, mfrow=c(7,6), spatial_list = Spatial_List, extrapolation_list = Extrapolation_List)

plot_anisotropy( FileName="Aniso.png", Report=Report, TmbData=TmbData )

Cov_List = FishStatsUtils::summarize_covariance(Report=Report, ParHat=Obj$env$parList(), Data=TmbData, SD=Opt$SD, plot_cor=TRUE,
                                                category_names=levels(Data_Geostat[,'spp']),plotTF=FieldConfig, mgp=c(2,0.5,0), tck=-0.02,
                                                oma=c(0,5,2,2) )


Cov_List = summarize_covar( Report=Report, ParHat=Obj$env$parList(), Data=TmbData, SD=Opt$SD, plot_cor=TRUE,
                            category_names=levels(Data_Geostat[,'spp']),plotTF=FieldConfig, mgp=c(2,0.5,0), tck=-0.02,
                            oma=c(0,5,2,2) )



#log predicted density
Dens_xt = FishStatsUtils::plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD,
                                    PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], 
                                    Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include,
                                    Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]],
                                    zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8,
                                    category_names = levels(full.1975$name))

#spatio-temporal var in encounter probability
Dens_xt_encounter = plot_maps(plot_set=6, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD,
                              PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], 
                              Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include,
                              Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]],
                              zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=TRUE,
                              category_names = levels(full.1975$name))


#spatio-temporal var in log-positive catch rates
Dens_xt_encounter = plot_maps(plot_set=7, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD,
                              PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], 
                              Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include,
                              Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]],
                              zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=TRUE,
                              category_names = levels(full.1975$name))

#
Dens_xt_encounter = plot_maps(plot_set=13, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD,
                              PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], 
                              Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include,
                              Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]],
                              zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=TRUE,
                              category_names = levels(full.1975$name))

View(Report$Epsilon2_gct)


Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)],
                 "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

Dens_DF[1:6,]

Index = plot_biomass_index(TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, 
                           category_names = levels(full.1975$name))

Index$Table[,c("Year","Fleet","Estimate_metric_tons","SD_log","SD_mt")]

plot_range_index(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), Year_Set=Year_Set, 
                 category_names = levels(full.1975$name))


MapDetails_List$fine_scale=FALSE

fac_save<-FishStatsUtils::plot_factors( Report=Report, ParHat=Obj$env$parList(), Data=TmbData, SD=Opt$SD, mapdetails_list=MapDetails_List, 
                                        Year_Set=Year_Set, category_names=levels(full.1975$name))


jpeg("all_loadings.jpeg", width=1000, height=1000)
par(mfrow=c(2,2))
plot(fac_save$Rotated_loadings$Omega1[,1], fac_save$Rotated_loadings$Omega1[,2], col="white", xlab="Factor 1", ylab="Factor 2", 
     cex.lab=1.5, cex.axis=1.5, main="Spatial Encounter", xlim=c(-2,2), ylim=c(-2,2))
abline(h=0, lwd=2, lty=3)
abline(v=0, lwd=2, lty=3)
text(fac_save$Rotated_loadings$Omega1[,1], fac_save$Rotated_loadings$Omega1[,2], labels=seq(from=1, to=4, by=1), cex=2)

plot(fac_save$Rotated_loadings$Omega2[,1], fac_save$Rotated_loadings$Omega2[,2], col="white", xlab="Factor 1", ylab="Factor 2", 
     cex.lab=1.5, cex.axis=1.5, main="Spatial Positive Density",xlim=c(-2,2), ylim=c(-2,2))
abline(h=0, lwd=2, lty=3)
abline(v=0, lwd=2, lty=3)
text(fac_save$Rotated_loadings$Omega2[,1], fac_save$Rotated_loadings$Omega2[,2], labels=seq(from=1, to=4, by=1), cex=2)

plot(fac_save$Rotated_loadings$Epsilon1[,1], fac_save$Rotated_loadings$Epsilon1[,2], col="white", xlab="Factor 1", ylab="Factor 2", 
     cex.lab=1.5, cex.axis=1.5, main="Spatio-temporal Encounter",xlim=c(-2,2), ylim=c(-2,2))
abline(h=0, lwd=2, lty=3)
abline(v=0, lwd=2, lty=3)
text(fac_save$Rotated_loadings$Epsilon1[,1], fac_save$Rotated_loadings$Epsilon1[,2], labels=seq(from=1, to=4, by=1), cex=2)

plot(fac_save$Rotated_loadings$Epsilon2[,1], fac_save$Rotated_loadings$Epsilon2[,2], col="white", xlab="Factor 1", ylab="Factor 2", 
     cex.lab=1.5, cex.axis=1.5, main="Spatio-temporal Positive Density",xlim=c(-2,2), ylim=c(-2,2))
abline(h=0, lwd=2, lty=3)
abline(v=0, lwd=2, lty=3)
text(fac_save$Rotated_loadings$Epsilon2[,1], fac_save$Rotated_loadings$Epsilon2[,2], labels=seq(from=1, to=4, by=1), cex=2)
dev.off()



beta2<-subset(Opt$diagnostics, Param=="beta2_ft")

plot(exp(beta2$MLE[1:35])~uni.yrs, type="l")

str(Report$mean_D_cyl)
jpeg("den_plot1.jpg", width = 350, height = 350)
plot(log(Report$D_gcy)~seq(from=1, to=42, by=1), type="l")
dev.off()


library(corrplot)
library(RColorBrewer)
M <-cor(mtcars)
corrplot(Report$lowercov_uppercor_epsilon2, type="lower")

med_fac<-apply(fac_save$Rotated_factors$Epsilon2[,1,], 2, median)

jpeg("temporal_epsilon2_1.jpg", width = 1000, height = 600)
plot(med_fac~seq(from=1977, to=2018, by=1), type="l", lwd=2, ylim=c(-3,3), ylab="Pos Density Factor 1", xlab="Year", cex.lab=1.25, cex.axis=1.5)
for(i in 1:266){
  lines(fac_save$Rotated_factors$Epsilon2[i,1,]~seq(from=1977, to=2018, by=1), col="grey", lwd=0.5)
}
lines(med_fac~seq(from=1977, to=2018, by=1), lwd=2)
points(med_fac~seq(from=1977, to=2018, by=1), pch=19, cex=1.5)
dev.off()

#range of spatial variation, whole bank is ~700 km N/s and 500 km E/W. distance at which locations have a correlation of 10%
Report$Range_raw1

#range of spatio-temporal variation
Report$Range_raw2

#Synchrony metrics are under phi
#this is for each spp
1-Save$Report$phi_cz

#this is for each knot
1-Save$Report$phi_gz

#likely the average across spp
1-Save$Report$phi_cbar_z

#likely the average across space
1-Save$Report$phi_gbar_z

#average across space and spp
1-Save$Report$phi_z

#spatial/species buffering index
Save$Report$maxsdB_cz

#coherence, ranges from 0 to 1, where 0 indicates all factors are equal, and 1 indicates that the first factor explains all variance
Report$psi


Report$Index_cyl
all.yrs<-seq(from=1977, to=2018, by=1)

jpeg("engel.jpeg")
par(mfrow=c(2,2))
plot(Report$Index_cyl[,,1][1,1:20]~all.yrs[1:20], type="l")
plot(Report$Index_cyl[,,1][2,1:20]~all.yrs[1:20], type="l")
plot(Report$Index_cyl[,,1][3,1:20]~all.yrs[1:20], type="l")
plot(Report$Index_cyl[,,1][4,1:20]~all.yrs[1:20], type="l")
dev.off()



knot_var<-matrix(ncol=50, nrow=34)
for(i in 1:50){
  knot_var[,i]<-rollapply(Report$D_gcy[i,1,], width = 2, FUN = sd, fill = NA)  
}

dev.off()
plot(knot_var[,33], type="l", ylim=c(0,50))
for(i in 2:50){
  lines(knot_var[,i])
}


plot(loc_sp, col="white")
text(loc_sp, labels=seq(from=1, to=50, by=1), col="blue")
