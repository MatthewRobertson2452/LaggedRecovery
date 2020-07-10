
load("converted_50knots_2factors.RData")
load("spatial_list.RData")
load("temp_mat.RData")
load("temp_yr.RData")



local_d<-Save$Report$D_gcy
global_d<-Save$Report$Index_cyl
locs<-Spatial_List$loc_x

#local_d_switch<-1 #turn to 1 if using local d and 0 if using global d
yt_density<-local_d[,2,1:33]
yt_global_density<-global_d[,,1][2,1:33]
ampl_density<-local_d[,1,1:33]
ampl_global_density<-global_d[,,1][1,1:33]


model <- 
  paste("
// To demonstrate AD. A function that would be difficult to differentiate
// analytically.

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(a);
  PARAMETER(b);
  //PARAMETER(c);
  PARAMETER(log_sigma);

  DATA_VECTOR(global_d);
  DATA_VECTOR(local_d);

  int n = local_d.size();

  //Type a = exp(log_a);
  //Type b = exp(log_b);
  Type sigma = exp(log_sigma);
  vector<Type> pred_d(n);
  vector<Type> resid(n);

  Type nll = 0.0;

  pred_d = a*pow(global_d,b);
  resid = local_d-pred_d;


  nll = -dnorm(resid, 0.0, sigma, true).sum();

  REPORT(a);
  REPORT(b);
  //REPORT(c);
  REPORT(log_sigma);

  return(nll);
}
")
writeLines(model,"exponential_fxn.cpp")


dyn.unload("exponential_fxn")
compile("exponential_fxn.cpp")
dyn.load("exponential_fxn")
obj <- MakeADFun(data=list(local_d=log(as.numeric(ampl_density)), global_d=rep(log(ampl_global_density), each=50)), 
                 parameters=list(a=1, b=1, log_sigma=log(5)), DLL='exponential_fxn')

fit<-nlminb(obj$par,obj$fn,,obj$gr)

fit$par

new_seq<-seq(from=0, to=10, length=10)

plot(data$local_d~data$global_d, pch=19, xlim=c(0,10))
lines(new_seq,fit$par[1]*new_seq^fit$par[2], col="red")

a<-2
b<-1.1
plot(a*new_seq^b, pch=19)


yt_local_d<-log(yt_density)+10
grid_cell<-seq(from=0, to=49, by=1)
yt_global_d<-log(yt_global_density)+10
nyr=33

model <- 
  paste("
        #include <TMB.hpp>
        template<class Type>
        Type objective_function<Type>::operator() ()
        {
        PARAMETER_VECTOR(log_a);
        PARAMETER_VECTOR(log_b);
        PARAMETER_VECTOR(log_sigma);
        
        DATA_MATRIX(local_d);
        DATA_VECTOR(global_d);
        DATA_INTEGER(n_grid);
        DATA_VECTOR(grid);
        DATA_INTEGER(nyr);

        //int n = local_d.size();
        
        vector<Type> a = exp(log_a);
        vector<Type> b = exp(log_b);
        vector<Type> sigma = exp(log_sigma);
        vector<Type> pred_d(n_grid);
        //vector<Type> resid(n_grid);
        vector<Type> resid_quick(n_grid);
        matrix<Type> resid(n_grid,nyr);
        
        Type nll = 0.0;
        Type zero = 0.0;
        using namespace density; //package for various multivariate Gaussian distribs
        
      for (int j=0; j<n_grid; j++){
        vector<Type> quick_local_d = local_d.row(j);

        pred_d = a(j)*pow(global_d, b(j));
        resid_quick = quick_local_d-pred_d;
        resid.row(j) = resid_quick;
        
        nll -= dnorm(resid_quick, zero, sigma(j), true).sum();
      }

        REPORT(log_a);
        REPORT(log_b);
        REPORT(log_sigma);
        REPORT(resid);
        ADREPORT(log_b);
        return(nll);
        }
        ")
writeLines(model,"exponential_fxn_RE.cpp")


dyn.unload("exponential_fxn_RE")
compile("exponential_fxn_RE.cpp")
dyn.load("exponential_fxn_RE")
obj <- MakeADFun(data=list(local_d=yt_local_d, global_d=yt_global_d, n_grid=50, grid=grid_cell, nyr=nyr), 
                 parameters=list(log_a=rep(log(1),50), log_b=rep(log(1),50), log_sigma=rep(log(3),50)), 
                 random=c("log_a","log_b"), DLL='exponential_fxn_RE')

fit<-nlminb(obj$par,obj$fn,,obj$gr)

fit$message

rep = obj$report()
sd.rep<-sdreport(obj)

save_list<-list(
  rep_obj=rep,
  sd_rep_obj=sd.rep
)
save(save_list, file="yt_re_output.RData")


b_mat<-matrix(nrow=50, ncol=4)
b_mat[,1] = exp(sd.rep$value)
b_mat[,2] =  sd.rep$sd
b_mat[,3]  = exp(sd.rep$value - qnorm(0.975)*b_mat[,2])
b_mat[,4]  = exp(sd.rep$value + qnorm(0.975)*b_mat[,2])

par(mfrow=c(1,1))
plot(seq(from=1, to=50, by=1), b_mat[,1], pch=19, ylim=c(-5,5))
abline(h=1, col="red")
segments(x0=seq(from=1, to=50, by=1), y0=b_mat[,1]-1.96*b_mat[,2] , x1=seq(from=1, to=50, by=1), y1=b_mat[,1]+1.96*b_mat[,2])


#if this is declining then it is 
plot(rowMeans(yt_local_d), b_mat[,1], pch=19, ylim=c(-2,10))


yt_resids<-

obj$gr(obj$par)

fit$par

new_seq<-seq(from=10, to=20, length=100)

par(mfrow=c(5,10),mar=c(0,0,1,0),oma=c(4,4,4,4))
for(i in 1:50){
plot(yt_local_d[i,]~yt_global_d, pch=19, main=i, xaxt="n", yaxt="n")
lines(new_seq,exp(rep$log_a[i])*new_seq^exp(rep$log_b[i]), col="red")
}


yt_slopes_df<-data.frame(long=locs[,1], lat=locs[,2], slopes_yt=exp(rep$log_b))
yt_slopes_df$id<-seq(from=1, to=50, by=1)

ggplot(yt_slopes_df, aes(x=long, y=lat, fill=slopes_yt))+
  geom_point(size=3, pch=21, colour="black")+
  scale_fill_gradient2(midpoint=1, low="blue", mid="white",
                       high="red")+
  labs(fill='Slope') +
  theme_classic()

yt_df<-data.frame(long=locs[,1], lat=locs[,2], rep$resid)
colnames(yt_df)[3:35]<-seq(from=1985, to=2017, by=1)
yt_df$id<-seq(from=1, to=50, by=1)

resid_df<-reshape(yt_df,
                  direction = "long",
                  varying = list(names(yt_df)[3:35]),
                  v.names = "Resid",
                  idvar = "id",
                  timevar = "Year",
                  times = 1985:2017)

ggplot(resid_df, aes(x=long, y=lat, fill=Resid))+
  geom_point(size=3, pch=21, colour="black")+
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                       high="red")+
  labs(fill='Residuals') +
  theme_classic()+
  facet_wrap(~Year, ncol = 6)



YR<-YR
yrs<-seq(from=1985, to=2017, by=1)

true_temp<-data.frame(yrs=YR, temp=t(temp_mat))
yrs<-as.data.frame(yrs)

merged_temps<-merge(true_temp, yrs, by="yrs", all.y=TRUE)

new_mat2<-matrix(nrow=35, ncol=33)
new_mat_p<-matrix(nrow=35, ncol=33)
for(k in 1:33){
  for(i in 1:35){
    new_mat2[i,k]<-cor.test(temp_mat[,i],as.numeric(yt_df[,3:35][,k]), method="spearman")$estimate
    new_mat_p[i,k]<-cor.test(temp_mat[,i],as.numeric(yt_df[,3:35][,k]), method="spearman")$p.value
    #if(new_mat_p[i,k]>0.05){new_mat2[i,k]<-NA}
  }
}

new_mat3<-new_mat2
for(i in 1:33){
  new_mat3[i,i:33]<-NA
}

cor_test1<-cor.test(temp_mat[,i],as.numeric(yt_df[,3:35][,k]), method="spearman")


par(mfrow=c(1,1))
layout.matrix <- matrix(c(1, 1,
                          2, 2), nrow = 2, ncol = 2, byrow = TRUE)

layout(mat = layout.matrix,
       heights = c(2, 1), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns
#layout.show(2)
par(mar=c(0,2,0,1))
image(flip(flip(raster(new_mat3),1),1), axes=F,col = heat.colors(12))
par(mar=c(2,2,0,1))
plot(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), bty="n",pch=19)
abline(v=2011, lwd=2, col="grey", lty=2)
abline(v=1991, lwd=2, col="grey", lty=2)
points(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), pch=19)
lines(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), lwd=2)




new_mat2<-matrix(nrow=35, ncol=33)
new_mat_p<-matrix(nrow=35, ncol=33)
for(k in 1:33){
  for(i in 1:35){
    new_mat2[i,k]<-cor.test(temp_mat[,i],ampl_density[,k], method="spearman")$estimate
    new_mat_p[i,k]<-cor.test(temp_mat[,i],ampl_density[,k], method="spearman")$p.value
    #if(new_mat_p[i,k]>0.05){new_mat2[i,k]<-NA}
  }
}

new_mat3<-new_mat2
for(i in 1:33){
  new_mat3[i,i:33]<-NA
}



par(mfrow=c(1,1))
layout.matrix <- matrix(c(1, 1,
                          2, 2), nrow = 2, ncol = 2, byrow = TRUE)

layout(mat = layout.matrix,
       heights = c(2, 1), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns
#layout.show(2)
par(mar=c(0,2,0,1))
image(flip(flip(raster(new_mat3),1),1), axes=F,col = heat.colors(12))
par(mar=c(2,2,0,1))
plot(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), bty="n",pch=19)
abline(v=2011, lwd=2, col="grey", lty=2)
abline(v=1991, lwd=2, col="grey", lty=2)
points(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), pch=19)
lines(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), lwd=2)



