Post\_VAST\_Analyses
================
Matt Robertson
10/07/2020

# Introduction

Fan and I have decided to divide the methods of this paper into four
sections:

1.  Run a multi-species VAST for yellowtail flounder and American plaice
2.  Test for the existence of density-dependent habitat selection (DDHS)
    by examining residuals from a single linear model
3.  If the initial test indicates DDHS, create knot specific models to
    account for it
4.  Correlate temperature distribution with either a) raw local density
    (i.e. distribution) data if no evidence of DDHS, or b) the residuals
    of knot-specific relationships if there is evidence of DDHS.

I will not detail the methods for the VAST because they have not
changed. Furthermore, details about new methods will be given
throughout.

# Testing for density-dependent habitat selection

Myers & Stokes (1989) as well as Shackell et al. (2005) tested for the
existence of density-dependent habitat selection by examining the
relationship between local and global biomass/abundance. They considered
the null hypothesis to be that variation in local density responses is
independent of habitat suitability (i.e. a single positive linear or
exponential relationship). This relationship would indicate that as
global biomass increases, so does density at all locations.

Alternatively, a density dependent habitat selection process (i.e. the
alternate hypothesis) would involve local population increasing at a
faster rate in marginal habitats than in prime/core habitats. This would
indicate that marginal habitats only become occupied when density
increases, while core areas are more stable to global biomass changes.
Therefore, a single global model of local density vs. global biomass
would violate assumptions of normality.

To apply this hypothesis test to the output of a VAST I will first
examine the relationship between all local densities (50 knots) and
global biomass each year. If there is a positive relationship, and there
is no evidence of model misfit (e.g. non-normal residuals would indicate
locally varying relationships) this will serve as evidence for the null
hypothesis (density-dependent habitat selection does not play a large
role in population distribution). If there is no relationship or
evidence of model misfit, this will serve as evidence for the
alternative hypothesis and we will test whether there are location
specific relationships between local density and global biomass.

## Yellowtail flounder

``` r
local_d<-Save$Report$D_gcy
global_d<-Save$Report$Index_cyl
locs<-Spatial_List$loc_x

yt_density<-local_d[,2,1:33]
yt_global_density<-global_d[,,1][2,1:33]
ampl_density<-local_d[,1,1:33]
ampl_global_density<-global_d[,,1][1,1:33]
```

A linear model initially appears to fit the yellowtail flounder data,
although the \(R^{2}\) value is quite low (\~0.07).

``` r
par(mar=c(4,4,0,0))
plot(log(as.numeric(yt_density))~log(rep(yt_global_density, each=50)), pch=19, xlab="log(global biomass)", ylab="log(local density)")
lm1<-lm(log(as.numeric(yt_density))~log(rep(yt_global_density, each=50)))
abline(a=lm1$coefficients[1], b=lm1$coefficients[2], col="red", lwd=2)
```

![Relationship between yellowtail flounder log(local density) and
log(global biomass). Solid red line is the predicted relationship from a
linear model. There are 50 local density estimates for every one global
biomass
estimate.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-3-1.png)

``` r
#setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Dissertation Plans\\Chapter 1\\VAST Manuscript\\Plots")
#png("yt_linmod.png", width = 4, height = 4, units = 'in', res=300)
par(mar=c(4,4,0,0))
plot(log(as.numeric(yt_density))~log(rep(yt_global_density, each=50)), pch=19, xlab="log(global biomass)", ylab="log(local density)")
lm1<-lm(log(as.numeric(yt_density))~log(rep(yt_global_density, each=50)))
abline(a=lm1$coefficients[1], b=lm1$coefficients[2], col="red", lwd=2)
```

![Relationship between yellowtail flounder log(local density) and
log(global biomass). Solid red line is the predicted relationship from a
linear model. There are 50 local density estimates for every one global
biomass
estimate.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-3-2.png)

``` r
#dev.off()
```

    ## 
    ## Call:
    ## lm(formula = log(as.numeric(yt_density)) ~ log(rep(yt_global_density, 
    ##     each = 50)))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -7.7084 -3.9295  0.7498  3.5695  7.9159 
    ## 
    ## Coefficients:
    ##                                        Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                            -11.8273     0.9538  -12.40   <2e-16 ***
    ## log(rep(yt_global_density, each = 50))   1.5833     0.1347   11.75   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.954 on 1648 degrees of freedom
    ## Multiple R-squared:  0.07732,    Adjusted R-squared:  0.07676 
    ## F-statistic: 138.1 on 1 and 1648 DF,  p-value: < 2.2e-16

When examining the residual histogram we can see a multimodal
distribution. This not only means that the residuals aren’t normally
distributed but also that there appear to be two unique distributions
underlying the poor fit. We see further evidence of model
misspecification in the qq plots as well.

![Residuals from the yellowtail flounder linear
model.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-5-1.png)

![Model diagnostic plots from the yellowtail flounder linear
model.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-6-1.png)

## American plaice

American plaice initally appears to have a more apparent linear trend
than we had observed for yellowtail flounder. This is further indicated
by the higher \(R^2\) value of \~0.38.

``` r
par(mar=c(4,4,0,0))
plot(log(as.numeric(ampl_density))~log(rep(ampl_global_density, each=50)), pch=19, xlab="log(global biomass)", ylab="log(local density)")
lm1<-lm(log(as.numeric(ampl_density))~log(rep(ampl_global_density, each=50)))
abline(a=lm1$coefficients[1], b=lm1$coefficients[2], col="red", lwd=2)
```

![Relationship between yellowtail flounder log(local density) and
log(global biomass). Solid red line is the predicted relationship from a
linear model. There are 50 local density estimates for every one global
biomass
estimate.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-7-1.png)

``` r
#setwd("C:\\Users\\mroberts\\OneDrive - Memorial University of Newfoundland\\Marine Institute\\Dissertation Plans\\Chapter 1\\VAST Manuscript\\Plots")
#png("ampl_linmod.png", width = 4, height = 4, units = 'in', res=300)
par(mar=c(4,4,0,0))
par(mar=c(4,4,0,0))
plot(log(as.numeric(ampl_density))~log(rep(ampl_global_density, each=50)), pch=19, xlab="log(global biomass)", ylab="log(local density)")
lm1<-lm(log(as.numeric(ampl_density))~log(rep(ampl_global_density, each=50)))
abline(a=lm1$coefficients[1], b=lm1$coefficients[2], col="red", lwd=2)
```

![Relationship between yellowtail flounder log(local density) and
log(global biomass). Solid red line is the predicted relationship from a
linear model. There are 50 local density estimates for every one global
biomass
estimate.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-7-2.png)

``` r
#dev.off()
```

    ## 
    ## Call:
    ## lm(formula = log(as.numeric(ampl_density)) ~ log(rep(ampl_global_density, 
    ##     each = 50)))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.84435 -0.58829  0.00699  0.55324  3.06245 
    ## 
    ## Coefficients:
    ##                                          Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                              -5.48268    0.25643  -21.38   <2e-16
    ## log(rep(ampl_global_density, each = 50))  1.15299    0.03606   31.98   <2e-16
    ##                                             
    ## (Intercept)                              ***
    ## log(rep(ampl_global_density, each = 50)) ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9172 on 1648 degrees of freedom
    ## Multiple R-squared:  0.3829, Adjusted R-squared:  0.3825 
    ## F-statistic:  1022 on 1 and 1648 DF,  p-value: < 2.2e-16

Unlike yellowtail flounder, the residuals for American plaice appear
normally distributed with no indication of model misspecification.

![Residuals from the American plaice linear
model.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-9-1.png)

![Model diagnostic plots from the American plaice linear
model.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-10-1.png)

## Conclusions

This simple test appears to indicate that yellowtail flounder does not
fit to the null hypothesis of proportional increase in local and global
biomass while American plaice does fit the null hypothesis. Therefore
indicating that American plaice does not appear to have
density-dependent habitat selection, while yellowtail flounder may.
Myers & Stokes (1989) had the same conclusion for American plaice and
Simpson & Walsh (2004) had a similar conclusion (using different
methods) for yellowtail flounder.

# Testing for locally varying relationships

To examine the potential for density-dependent habitat selection in
yellowtail flounder, I created a non-linear random effects model in TMB
using, \[\hat{y}_{k}=a_kx^{b_k}\] to examine local variability in the
density-dependent habitat selection relationship (Myers & Stokes, 1989).
Where \(\hat{y}_{k}\) is the log of knot density +10, with the +10
ensuring that all local densities are positive, and \(x\) representing
the log of total biomass +10. Both \(a_k\) and \(b_k\) are random
effects to provide unique estimates for all knots.

We would expect that the estimate of \(b_K\) for locations would be
shallowest/concave in prime/core habitats and highest/convex in marginal
habitats. Weaker responses (i.e. shallow \(b_k\)) would specifically
indicate that a location is stable with regional biomass changes.

The model converged well, where the gradient of the approximated
marginal log-likelihood for all fixed effects was \<10-6 and the Hessian
matrix was positive definite at the maximum-likelihood estimates. We can
observe the specific fits separately here,

![Yellowtail flounder local density based on global biomass across the
fifty knots (panels) in VAST. The red lines represent the predicted
local relationship from the random effects model and grey polygon
represents the root mean square error for those
relationships.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-12-1.png)![Yellowtail
flounder local density based on global biomass across the fifty knots
(panels) in VAST. The red lines represent the predicted local
relationship from the random effects model and grey polygon represents
the root mean square error for those
relationships.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-12-2.png)

We can then examine the spatial trends for the slopes. Again, we would
expect higher slopes (and/or concave relationships) for marginal
habitats and shallow slopes (and/or convex relationships) for optimal
habitats. Figure 3.2 shows that the shallowest slopes are close to the
southeast shoal with values increasing with increased distance from that
area. Therefore indicating that ideal habitat for yellowtail flounder is
on the southeast shoal. Again, this matches the findings of Simpson &
Walsh (2004).

![The estimate of \(b_k\) from the knot specific exponential models
comparing yellowtail flounder local density based on global biomass
across the fifty knots in the
VAST.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-13-1.png)![The
estimate of \(b_k\) from the knot specific exponential models comparing
yellowtail flounder local density based on global biomass across the
fifty knots in the
VAST.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-13-2.png)

We can also examine whether the slopes are truly significantly different
than 1 and/or 0.

![Model estimates of \(b_k\) at each knot. The segments represent +/-
1.96\*SD. The solid red line indicates 1, the dashed black line
indicates
0.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-14-1.png)![Model
estimates of \(b_k\) at each knot. The segments represent +/- 1.96\*SD.
The solid red line indicates 1, the dashed black line indicates
0.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-14-2.png)

## Conclusions

There were locally varying relationships between local density and
global biomass, indicating that yellowtail flounder is influenced by
density-dependent habitat selection. Furthermore, it appears that
optimal habitat for yellowtail flounder are centered on the southeast
shoal and the most marginal habitat is in division 3L.

# Distribution shifts

## Removing the effect of density-dependent habitat selection

Yellowtail flounder distribution appears to be influenced by
density-dependent habitat selection. To ensure that our results are not
affected by this mechanism we can examine the residuals from this
relationship. They should represent deviations from the
density-dependent relationship (i.e. the density-independent influence
on distribution).

Figure 4.1 shows a general trend of large positive residuals in 3L from
1985-1992, followed by large negative residuals in 3L until \~2004,
relatively low residuals until 2010, and then large positive residuals
in 3L onwards. Residuals outside of 3L are generally lower and do not
tend to show strong temporal correlation. Therefore indicating that the
largest variations from the predictions from the density-dependent
habitat selection model occur in the marginal habitats which may be most
influenced by environmental conditions.

![Yellowtail flounder local density residuals from the density dependent
habitat selection
model](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-15-1.png)![Yellowtail
flounder local density residuals from the density dependent habitat
selection
model](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-15-2.png)

## Comparing distributions through space and time

Since American plaice were not influenced by density-dependent habitat
selection, and we were able to remove the effect of density-dependent
habitat selection from yellowtail flounder, we can directly examine the
density-independent distribution of both species.

First, we will examine self-correlation through time for temperature,
and then both fish distributions to understand how the distributions
have changed in relative terms.

Temperature has been highly correlated (\>0.6) throughout the entire
time-series (Figure 4.2). There are some years where the correlation is
lower than average, however it appears that despite mean temperatures
changing, the distribution of temperature has remained relatively
constant.

``` r
YR<-YR
yrs<-seq(from=1985, to=2017, by=1)

true_temp<-data.frame(yrs=YR, temp=t(temp_mat))
yrs<-as.data.frame(yrs)

merged_temps<-merge(true_temp, yrs, by="yrs", all.y=TRUE)

new_mat2<-matrix(nrow=35, ncol=35)
new_mat_p<-matrix(nrow=35, ncol=35)
rev_seq<-seq(from=35, to=1)
for(k in 1:35){
  for(i in 1:35){
    new_mat2[rev_seq[i],k]<-cor.test(temp_mat[,i],temp_mat[,k], method="spearman")$estimate
    new_mat_p[rev_seq[i],k]<-cor.test(temp_mat[,i],temp_mat[,k], method="spearman")$p.value
    #if(new_mat_p[rev_seq[i],k]>0.05){new_mat2[rev_seq[i],k]<-NA}
  }
}

new_mat3<-new_mat2
for(i in 1:35){
  new_mat3[i,1:(35-i)]<-NA
}


par(mfrow=c(1,1), mar=c(0,0,0,4))

plot(raster(new_mat3), axes=F,col = heat.colors(12), bty="n", box=FALSE)


par(mfrow=c(1,1), mar=c(0,0,0,4))

plot(raster(new_mat3), axes=F,col = heat.colors(12), bty="n", box=FALSE)
```

![Spearman correlation of temperature distribution with temperature
distribution across
time](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-16-1.png)

Yellowtail flounder show three distinct periods, an initial state that
remained constant for about 7 years, followed by new state that was
negatively correlated with the inital state and lasted for \~20 years,
and the most recent state which again appears correlated with the first
state but not with the second (Figure 4.3).

![Spearman correlation of yellowtail flounder distribution with
yellowtail flounder distribution across
time](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-17-1.png)

American plaice appear to have two distributional states, an initial
state that was correlated with itself for 5-7 years and a second that
was negatively correlated with the first state and has persisted to this
day (Figure 4.4).

![Spearman correlation of American plaice distribution with American
plaice distribution across
time](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-18-1.png)

## Comparing distribution shifts with temperature distributions

Since the distribution of temperature has remained relatively contant
through time it will serve as a spatial reference. Where positive
correlations with fish distributions will indicate that most fish are
found in warm water, and negative correlations will indicate that most
fish are found in cold water. If fish distributions have non-lagged
patterns then each temperature condition (denoted by a mean temperature)
should have the same correlation with fish distribution. Therefore, we
would predict that extreme temperatures may generate a change in
correlation, but that soon after this change, the fish distribution
should return when temperatures return.

Yellowtail flounder show a southward distributional shift after 1991
that was maintained until temperatures hit their highest recorded value
in 2010 when they moved northwards again (Figure 4.5). This pattern
indicates that yellowtail flounder showed a persistent distributional
shift following a decline in temperature.

``` r
legend.col <- function(col, lev){
  
  opar <- par
  
  n <- length(col)
  
  bx <- par("usr")
  
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx+0.02, yy, col = col[i], border = "black")
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev,na.rm=TRUE), max(lev,na.rm=TRUE)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .75)
  par <- opar
}
```

``` r
new_mat2<-matrix(nrow=35, ncol=33)
new_mat_p<-matrix(nrow=35, ncol=33)
rev_seq<-seq(from=35, to=1)
for(k in 1:33){
  for(i in 1:35){
    new_mat2[rev_seq[i],k]<-cor.test(temp_mat[,i],as.numeric(yt_df[,3:35][,k]), method="spearman")$estimate
    new_mat_p[rev_seq[i],k]<-cor.test(temp_mat[,i],as.numeric(yt_df[,3:35][,k]), method="spearman")$p.value
    #if(new_mat_p[rev_seq[i],k]>0.05){new_mat2[rev_seq[i],k]<-NA}
  }
}

new_mat3<-new_mat2
for(i in 1:29){
  new_mat3[i,1:(29-i)]<-NA
}


par(mfrow=c(1,1))
layout.matrix <- matrix(c(1, 1,
                          2, 2), nrow = 2, ncol = 2, byrow = TRUE)

layout(mat = layout.matrix,
       heights = c(2, 1), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns

par(mar=c(0,2,0,4))
image(raster(new_mat3), axes=F,col = heat.colors(12))
legend.col(col = heat.colors(12), lev = new_mat3)
par(mar=c(2,2,0,4))
plot(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), bty="n",pch=19)
abline(v=2011, lwd=2, col="grey", lty=2)
abline(v=1991, lwd=2, col="grey", lty=2)
points(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), pch=19)
lines(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), lwd=2)
```

![Spearman correlation of yellowtail flounder distribution (x) with
temperature distribution (y) across time. Time-series on the bottom
panel represents mean temperature across the Grand Bank during the
survey. Years for temperature and fish distribution are not matched 1:1
due to earlier data for temperature and missing years for temperature.
The x-axis in the time-series matches the correlation
plot.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-20-1.png)

``` r
par(mfrow=c(1,1))
layout.matrix <- matrix(c(1, 1,
                          2, 2), nrow = 2, ncol = 2, byrow = TRUE)

layout(mat = layout.matrix,
       heights = c(2, 1), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns
par(mar=c(0,2,0,4))
image(raster(new_mat3), axes=F,col = heat.colors(12))
legend.col(col = heat.colors(12), lev = new_mat3)
par(mar=c(2,2,0,4))
plot(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), bty="n",pch=19)
segments(x0=2011, y0=0, x1=2011, y1=2.8, lwd=2, col="grey", lty=2)
segments(x0=1991, y0=0, x1=1991, y1=2.8, lwd=2, col="grey", lty=2)
points(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), pch=19)
lines(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), lwd=2)
```

![Spearman correlation of yellowtail flounder distribution (x) with
temperature distribution (y) across time. Time-series on the bottom
panel represents mean temperature across the Grand Bank during the
survey. Years for temperature and fish distribution are not matched 1:1
due to earlier data for temperature and missing years for temperature.
The x-axis in the time-series matches the correlation
plot.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-20-2.png)

American plaice also show a southward distributional shift after 1991,
however that shift has never been fully reversed (Figure 4.6). This
pattern indicates that American plaice also showed a persistent
distributional shift following the decline in temperature and that their
threshold value to return to the north is higher than the threshold for
yellowtail flounder.

``` r
new_mat2<-matrix(nrow=35, ncol=33)
new_mat_p<-matrix(nrow=35, ncol=33)
rev_seq<-seq(from=35, to=1)
for(k in 1:33){
  for(i in 1:35){
    new_mat2[rev_seq[i],k]<-cor.test(temp_mat[,i],ampl_density[,k], method="spearman")$estimate
    new_mat_p[rev_seq[i],k]<-cor.test(temp_mat[,i],ampl_density[,k], method="spearman")$p.value
    #if(new_mat_p[rev_seq[i],k]>0.05){new_mat2[rev_seq[i],k]<-NA}
  }
}

new_mat3<-new_mat2
for(i in 1:29){
new_mat3[i,1:(29-i)]<-NA
}

par(mfrow=c(1,1))
layout.matrix <- matrix(c(1, 1,
                          2, 2), nrow = 2, ncol = 2, byrow = TRUE)

layout(mat = layout.matrix,
       heights = c(2, 1), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns

par(mar=c(0,2,0,4))
image(raster(new_mat3), axes=F,col = heat.colors(12))
legend.col(col = heat.colors(12), lev = new_mat3)
par(mar=c(2,2,0,4))
plot(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), bty="n",pch=19)
abline(v=2011, lwd=2, col="grey", lty=2)
abline(v=1991, lwd=2, col="grey", lty=2)
points(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), pch=19)
lines(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), lwd=2)
```

![Spearman correlation of American plaice distribution (x) with
temperature distribution (y) across time. Time-series on the bottom
panel represents mean temperature across the Grand Bank during the
survey. Years for temperature and fish distribution are not matched 1:1
due to earlier data for temperature and missing years for
temperature.The x-axis in the time-series matches the correlation
plot.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-21-1.png)

``` r
par(mfrow=c(1,1))
layout.matrix <- matrix(c(1, 1,
                          2, 2), nrow = 2, ncol = 2, byrow = TRUE)

layout(mat = layout.matrix,
       heights = c(2, 1), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns

par(mar=c(0,2,0,4))
image(raster(new_mat3), axes=F,col = heat.colors(12))
legend.col(col = heat.colors(12), lev = new_mat3)
par(mar=c(2,2,0,4))
plot(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), bty="n",pch=19)
segments(x0=2011, y0=0, x1=2011, y1=2.8, lwd=2, col="grey", lty=2)
segments(x0=1991, y0=0, x1=1991, y1=2.8, lwd=2, col="grey", lty=2)
points(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), pch=19)
lines(yrs$yrs, rowMeans(merged_temps[,2:51], na.rm = TRUE), lwd=2)
```

![Spearman correlation of American plaice distribution (x) with
temperature distribution (y) across time. Time-series on the bottom
panel represents mean temperature across the Grand Bank during the
survey. Years for temperature and fish distribution are not matched 1:1
due to earlier data for temperature and missing years for
temperature.The x-axis in the time-series matches the correlation
plot.](post_vast_analyses_markdown_files/figure-gfm/unnamed-chunk-21-2.png)

# Overall conclusions

Yellowtail flounder appears to be influenced by density-dependent
habitat selection while American plaice does not. Ideal habitat for
yellowtail flounder is on the southeast shoal with marginal habitat in
the north (division 3L). Since American plaice are not affected by
density-dependent habitat selection, it does not need to be accounted
for in later analyses. Meanwhile, we can examine the spatial
distribution of yellowtail flounder without the effects of
density-dependent habitat selection by examining the residuals from the
density-dependent relationship. These analyses show that both American
plaice and yellowtail flounder shifted their distributions following the
coldest point in the time-series (1991). Neither species shifted back,
even when temperatures returned to their previous values. Yellowtail
only shifted back to a similar distribution once temperatures reached
their maximum in 2011, while American plaice has still never shifted
back.

# References

Myers, R.A. and Stokes, K. (1989). Density-dependent habitat utilization
of groundfish and the improvement of research survey. ICES Committee
Meeting D15.

Simpson, M.R. and Walsh, S.J. (2004). Changes in the spatial structure
of Grand Bank yellowtail flounder: Testing MacCall’s basin hypothesis.
Journal of Sea Research, 51(3-4), 199-210.

Shackell, N.L., Frank, K.T., and Brickman, D.W. (2005). Range
contraction may not always predict core areas: an example from marine
fish. Ecological applications, 15(4), 1440-1449.
