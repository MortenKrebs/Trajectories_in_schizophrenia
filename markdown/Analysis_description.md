This is a document providing the supplementary figures and tables for the paper "Biological Predictors of Psychiatric Trajectories in a Danish Population Cohort of Individuals with Schizophrenia" by Krebs et al.

R code is available at [<https://github.com/MortenKrebs/Trajectories_in_schizophrenia>](https://github.com/MortenKrebs/Trajectories_in_schizophrenia).

The R-package `diagtract` can downloaded from [<https://github.com/MortenKrebs/diagtraject>](https://github.com/MortenKrebs/diagtraject).

Sequence Analysis in Schizophrenia
==================================

Data Material
-------------

``` r
rm(list=ls())

library(rmarkdown)
library(knitr)
library(kableExtra)
library(devtools)

library(data.table)
library(dplyr)


library(mstate)
library(TraMineR)
library(TraMineRextras)
library(vegan)
library(WeightedCluster)
library(fpc)
library(WGCNA)
```

    ## ==========================================================================
    ## *
    ## *  Package WGCNA 1.51 loaded.
    ## *
    ## *    Important note: It appears that your system supports multi-threading,
    ## *    but it is not enabled within WGCNA in R. 
    ## *    To allow multi-threading within WGCNA with all available cores, use 
    ## *
    ## *          allowWGCNAThreads()
    ## *
    ## *    within R. Use disableWGCNAThreads() to disable threading if necessary.
    ## *    Alternatively, set the following environment variable on your system:
    ## *
    ## *          ALLOW_WGCNA_THREADS=<number_of_processors>
    ## *
    ## *    for example 
    ## *
    ## *          ALLOW_WGCNA_THREADS=28
    ## *
    ## *    To set the environment variable in linux bash shell, type 
    ## *
    ## *           export ALLOW_WGCNA_THREADS=28
    ## *
    ## *     before running R. Other operating systems or shells will
    ## *     have a similar command to achieve the same aim.
    ## *
    ## ==========================================================================

``` r
library(class)


library(ggplot2)
library(corrplot)
library(png)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(scales)
```

``` r
library(devtools)
install_github("MortenKrebs/diagtraject")
library(diagtraject)
```

The Study was performed using the iPSYCH cohort (See [Pedersen et al](https://www.nature.com/articles/mp2017196)).

``` r
# 3 data.frames with one row pr individual and collums containing: 
# id, birthday and ##date of first diagnosis for each category of diagnoses

load("diagdates.Rda") # individuals diagnosed with Schizophrenia before Dec 31, 2012  
load("diagdates_random.Rda") # random population sample
load("diagdates2016.Rda") # individuals diagnosed with Schizophrenia 2013-2016
```

The primary cohort included all indiviuals born in Denmark between May 1, 1981 and December 31, 2002, and diagnosed with schizophrenia before December 31, 2012. (N= 5432)

``` r
colnames(diagdates)
##  [1] "pid"       "birthdate" "F1"        "F3"        "F4"       
##  [6] "F50"       "F60"       "F70"       "F84"       "F9"       
## [11] "SCZ"       "gender"
nrow(diagdates)
## [1] 5432

nrow(diagdates2016)
## [1] 870


dt <- data.table(diagdates)


fdate <- function(x) as.Date(x,"%d/%m/%Y" )

dt[,"birthdate":= fdate(birthdate) ]
dt[,c("F1","F3", "F4","F50","F60", "F70", "F84", "F9","SCZ") :=
     lapply(list(F1,F3,F4,F50,F60,F70,F84,F9,SCZ),function(x){  
       (as.numeric(fdate(x))- as.numeric(birthdate))/365})]

#Calculate age 31/12/2016
dt[,censored:= (as.numeric(fdate("31/12/2016"))-as.numeric(birthdate))/365] 
```

An additional non-overlapping cohort of indiviuals born in Denmark between May 1, 1981 and December 31, 2005, and diagnosed with Schizophrenia before December 31, 2016 was obtained for replication. (N=870).

``` r
## Chisquare of comorbidity:
# 
# summary(xtabs(N~.,dt[,.N,.(F1=!is.na(F1),F3=!is.na(F3),
#                            F4=!is.na(F4),F50=!is.na(F50),
#                            F60=!is.na(F60),F70=!is.na(F70),
#                            F84=!is.na(F84),F9=!is.na(F9))]))
# 


# using the 'mstate' library
mat <-pm <-  matrix(NA, 8,8)
tmat <- transMat(x = list(2:3,3, c()), names = c("None",  "D1", "D2"))
covs <- c("gender","birthdate")

for(i in 1:8 ){ 
dt_s <- dt
dt_s$x <- dt_s[,c("F1","F3", "F4","F50","F60", "F70", "F84", "F9")[i],with=F]
dt_s[,D1stat :=as.numeric(!is.na(x))]
dt_s[D1stat==0,D1time := censored]
dt_s[D1stat==1,D1time := x]

for( j in 1:8 ) if(!i==j)  { 
dt_s$x <- dt_s[,c("F1","F3", "F4","F50","F60", "F70", "F84", "F9")[j],with=F]
dt_s[,D2stat :=as.numeric(!is.na(x))]
dt_s[D2stat==0,D2time := censored]
dt_s[D2stat==1,D2time := x]
dt_s<- dt_s[!(D1stat==1 & D2stat==1 & D1time==D2time)]

msbmt <- msprep(time = c(NA, "D1time", "D2time"), status = c(NA,  "D1stat", "D2stat"), data = dt_s, trans = tmat, keep = covs)
cox1 <- coxph(Surv(Tstart, Tstop, status) ~I(trans==3)+gender+birthdate, data = msbmt, subset=trans==2|trans==3, method = "breslow")
mat[i,j]<- exp(cox1$coefficients[1])
pm[i,j] <- summary(cox1)$coefficients[1,5]
}}


mat2 <-log(mat)
pm[pm<1e-6] <- 1e-6

col <- colorRampPalette(c("#67001F", "#B2182B", "#FFFFFF", "#2166AC", "#053061")[5:1])(20)
diag(mat2) <- 0
diag(pm) <- 1
rownames(mat2)<- colnames(mat2) <- c("F1","F3", "F4","F50","F60", "F70", "F84", "F9")


knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
Color is the log(HR) of disorder 2 (x-axis) given disorder 1 (y axis). Estimates are based on survival analysis with disorder 1 as time-dependent covariat and adjusting for age and gender. Individuals with ties (=diagnosis 1 and 2 given on the same date) are excluded. * indicates p< 0.05/56",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})


par(mfrow=c(1,1))
corrplot(mat2, is.corr=FALSE, p.mat=pm,
                      method="p.square", insig="pch", col=col,
                      sig.level=0.05/56, pch="*", pch.cex=1, mar= c(4,0,4,0),
                      cl.ratio=0.25, cl.align.text="l", cl.lim=c(-1.5,1.5), 
                      tl.col="black", bg="#fffffc",ssl.pos = "r",ssl.lev=6,full_col=F)
text(-1.5,5,labels = "First Diagnosis",srt=90)
text(9.5,9,labels = "HR")
mtext("Second Diagnosis",side=3,line=1)
```


Diagnosis State Sequences:
--------------------------

Using increaments of 12 month, comorbidity trajectories were transformed to sequences. In case of multiple comorbidities, alphabet extensions were used.

``` r
dt.m <-melt(dt, id.vars = 'pid', direction = "long", 
            measure.vars  = list(c(3:10,13)), value.name = "time",
            variable.name = "event")
dt.m<- dt.m[order(pid)]
dt.m[,time:=time]
dt.m <- dt.m[!is.na(time)]

dt.mm <- dt.m 

events <- levels(dt.m$event)
drop <-  matrix(FALSE, nrow = length(events), ncol = length(events),
                dimnames = list(events,events))
drop['censored',] <- T
drop[,'censored'] <- T
diag(drop) <- F

e2sm <-seqe2stm(events = events, dropMatrix = drop)
```

``` r
seq_mis <- seqdef(TSE_to_STS(dt.m,
                             id = "pid", timestamp = "time", event= "event",
                             tmin=1, tmax=ceiling(max(dt.m$time)),stm = e2sm)
                  , firstState ="None", missing="censored")
```

``` r
n <- nrow(seq_mis);n
## [1] 5432
ran <- range(seqlength(seq_mis));ran
## [1] 16 35
alp <-length(alphabet(seq_mis));alp
## [1] 193
```

Among the 5432 included subjects, we observe trajectories ranging from 16 to 35 years. The alphabet include 193 different states.

``` r
setkey(dt,pid)
seq25 <- seq_mis[seqlength(seq_mis)>=25 & dt[,SCZ] < 25,]
seq25 <- seq25[,seq(1,25,3)]

cols <- brewer.pal(9, "Set1")
tab <-seqtab(seq25,tlim = 1:20)
lab <- attributes(seq25)$labels[
  which(lapply(attributes(seq25)$labels, function(x) sum(
    grepl(x, rownames(attributes(tab)$weights))))>0)]
o<-c(1,4,6,7,2,3,5,8)
lab[o][1:4] <- paste(lab[o][1:4], 
                     c("Substance Abuse", "Mood Disord.",
                       "Neurot. Disord.", "Personal. Disord."), 
                     sep= " - ")


attributes(seq25)$cpal <- rep("grey", length(attributes(seq25)$labels) )
attributes(seq25)$cpal[which(
  lapply(attributes(seq25)$labels, function(x) sum(
    grepl(x, rownames(attributes(tab)$weights))))>0)]<- c(
      cols[c(1,2,5,3,7,4,6)], "white")

knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
Top-20 most frequent sequences in first 25 years for individuals followed more than 25 years.",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})
par(mfrow=c(1,2))
seqfplot(seq25,tlim=20:1, pbarw=F, with.legend=F, yaxis="pct")
plot(0,type='n',axes=FALSE,ann=FALSE)
legend(x=0.7,y=1, lab[o],
      fill = c(cols[c(1,2,5,3,7,4,6)][o], "white"),cex=.8)
```


``` r
# For later use a version keeping `"censored"` as a seperate state was created:

seq <- seqdef(TSE_to_STS(dt.m,id = "pid", timestamp = "time", 
                         event="event",
                         tmin=1, tmax=ceiling(max(dt.m$time)),
                         stm = e2sm),
              firstState ="None")
```

``` r
#Identifying the state at time of censoring: 

drop['censored',] <- F
drop[,'censored'] <- T
e2sm <-seqe2stm(events = events, dropMatrix = drop)

seq_last <- seqdef(
  TSE_to_STS(dt.m, id = "pid", timestamp = "time", 
             event= "event", tmin=1,tmax=ceiling(max(dt.m$time))+1,
             stm = e2sm),
  firstState ="None", missing="censored")
seq_mis_last <- seq_mis
seq_mis_last$last <- seq_last[,ncol(seq_last)]
attributes(seq_mis_last)$alphabet <- c(alphabet(seq_mis_last),  
                                       unique(seq_mis_last$last)[which(
                                         !unique(seq_mis_last$last) %in%
                                           alphabet(seq_mis_last))])
attributes(seq_mis_last)$alphabet <- attributes(
  seq_mis_last)$alphabet[order(attributes(seq_mis_last)$alphabet)]
attributes(seq_mis_last)$labels <- attributes(seq_mis_last)$alphabet
```

Sequence Dissimilarities:
-------------------------

### Transition Rates:

To obtain population-wide estimates of transition probabilities, the random population cohort of 30000 individuals was included (See [Pedersen et al](https://www.nature.com/articles/mp2017196) for details) and estimates were computed using inverse sampling probability weighting (Bowman II).

``` r
load("diagdates_random.Rda")
load("diagdates2016.Rda")
```

``` r
dt_pop <- data.table(rbind(diagdates, diagdates2016, diagdates_random))

dt_pop[,year := format(birthdate,'%Y')]
dt_pop[,year := as.factor(as.numeric(year))]
sample_dist <- dt_pop[pid %in% diagdates_random$pid,.N,c("year","gender")]
```

To estimate sampling probability, the distribution of age and sex in the Danish population was obtained (publicly available from [Statistics Denmark](www.statistikbanken.dk)).

``` r
DK_pop <- data.table("year"=1981:2005,
                     "M"= c(27117*7/12,27063,26001,26572,
                            27465,28434,29079,30324,31475,
                            32620,33005,34812,34609,35639,
                            35886,34819,34749,34059,33879,
                            34432,33497,32964,33167,33077,32827),
                     "F"= c(25972*7/12,25595,24821,25228,
                            26284,26878,27142,28520,29876,
                            30813,31353,32914,32760,34027,
                            33885,32819,32899,32115,32341,
                            32652,31961,31111,31432,31532,31455))
DK_pop <-melt(DK_pop, id.vars = "year", value.name = "n")
DK_pop <-DK_pop[year<2002]
```

Inverse sampling probability weighting was performed taking age and sex into account.

``` r
sample_dist <- sample_dist[order(gender,year)]
sample_dist[,pop_weights := DK_pop[,n]/sample_dist[,N]]
dt_pop <- merge(dt_pop, sample_dist, by=c("year","gender"))
dt_pop[!is.na(SCZ), pop_weights := 1]
setkey(dt_pop,pid)
```

``` r

dt_pop[,"birthdate":= fdate(birthdate) ]
dt_pop[,c("F1","F3", "F4","F50","F60", "F70", "F84", "F9","SCZ") :=
     lapply(list(F1,F3,F4,F50,F60,F70,F84,F9,SCZ),function(x){  
       (as.numeric(fdate(x))- as.numeric(birthdate))/365})]
dt_pop[,censored:= (as.numeric(fdate("31/12/2016"))-
                      as.numeric(birthdate))/365] 
dt.m <-melt(dt_pop, id.vars = 'pid', direction = "long", 
            measure.vars  = list(c(5:12,16)), value.name = "time",
            variable.name = "event")
dt.m<- dt.m[order(pid)]
dt.m[,time:=time]
dt.m <- dt.m[!is.na(time)]

events <- levels(dt.m$event)
drop <-  matrix(FALSE, nrow = length(events), 
                ncol = length(events), dimnames = list(events,events))
drop['censored',] <- T
drop[,'censored'] <- T
diag(drop) <- F

e2sm <-seqe2stm(events = events, dropMatrix = drop)

seq_pop <- seqdef(TSE_to_STS(dt.m,id = "pid", 
                             timestamp = "time", event= "event",
                             tmin=1, tmax=ceiling(max(dt.m$time)),
                             stm = e2sm),
                  firstState ="None")


seq_mis_pop <- seqdef(TSE_to_STS(dt.m, id = "pid", 
                                 timestamp = "time", event= "event",
                                 tmin=1, tmax=ceiling(max(dt.m$time)),
                                 stm = e2sm), 
                      firstState ="None", missing="censored")

drop['censored',] <- F
drop[,'censored'] <- T
e2sm <-seqe2stm(events = events, dropMatrix = drop)

seq_last <- seqdef(
  TSE_to_STS(dt.m, id = "pid", timestamp = "time", 
             event= "event", tmin=1,tmax=ceiling(max(dt.m$time))+1,
             stm = e2sm),
  firstState ="None", missing="censored")
seq_mis_last_pop <- seq_mis_pop
seq_mis_last_pop$last <- seq_last[,ncol(seq_last)]
attributes(seq_mis_last)$alphabet <- c(alphabet(seq_mis_last),  
                                       unique(seq_mis_last$last)[which(
                                         !unique(seq_mis_last$last) %in% 
                                           alphabet(seq_mis_last))])
attributes(seq_mis_last)$alphabet <- attributes(
  seq_mis_last)$alphabet[order(attributes(seq_mis_last)$alphabet)]
attributes(seq_mis_last)$labels <- attributes(seq_mis_last)$alphabet
```

Time-varying Transition Rates were computed with `TraMineR::seqtrate()`.

``` r
seq_mis_last_uni_pop <- unique(seq_mis_last_pop)
attributes(seq_mis_last_uni_pop)$weights <- match(
  seqconc(seq_mis_last_pop),seqconc(seq_mis_last_uni_pop))

tr_t <- seqtrate(seq_mis_last_uni_pop, time.varying = T,weighted = T)
```

### Substitution Costs:

To accomodate the large sequence alphabet substitution cost were selcted to take into account overlap in states using Jaccard distance:

\[ d_J(A,B) = 1-{{|A \cap B|}\over{|A \cup B|}} \]

``` r
alp <- seqdef(as.data.frame(alphabet(seq_pop)),stsep="[.]")
rownames(alp) <- alphabet(seq_pop)
sub.cost_jacc <- matrix(0, nrow(alp),nrow(alp))
alp <- alp[order(seqlength(alp)),]
for(i in 1:nrow(alp))
  for(j in 1:nrow(alp))
    if(i>=j){
      sub.cost_jacc[i,j] <-  1-
        sum(unlist(alp[i,1:seqlength(alp[i,])]) %in% 
              unlist(alp[j,1:seqlength(alp[j,])]))/
        length(unique(c(unlist(alp[i,1:seqlength(alp[i,])]), 
                        unlist(alp[j,1:seqlength(alp[j,])]))))} else {
                          sub.cost_jacc[i,j] <- NA}
sub.cost_jacc <- as.matrix(as.dist(sub.cost_jacc ))
rownames(sub.cost_jacc)  <- colnames(sub.cost_jacc) <- rownames(alp)
sub.cost_jacc <- sub.cost_jacc[order(rownames(sub.cost_jacc)),
                               order(rownames(sub.cost_jacc))]
sub.cost <- seqsubm(seqdata = seq_pop, method = "CONSTANT")
rownames(sub.cost_jacc) <- rownames(sub.cost)
colnames(sub.cost_jacc) <- colnames(sub.cost)
sub.cost_jacc["censored->",] <- sub.cost_jacc[,"censored->"] <- 0

dim(sub.cost_jacc)
## [1] 202 202
```

``` r
kable(sub.cost_jacc[1:5,1:5],format="latex",booktabs=T, caption = "Substitution cost for the first five states in the alphabet")
```

### Right-censoring:

To calculate dissimilarities between sequences with right censoring, we used inferred states weighted by the probabilities of that state, given the last observed state in the sequence. When calculating the dissimilarity between two sequences \(i\) and \(j\) of unequal length, the dissimilarity, \(D(i,j)\) was defined as:
\[ D(i,j)  = d_{obs} + d_{inf} \] where \(d_{obs}(i,j)\) is the dissimilarity between the sequences before right censoring occurs and \(d_{inf}\) is the sum of the substitution cost matrix weighted by the inferred probabilities.

``` r
  sms <- which(rownames(sub.cost_jacc) %in% 
                 rownames(seqsubm(seq,"CONSTANT")))
  d_OM <- seqdist(seq, method = "OM", indel=.5, 
                  sm = sub.cost_jacc[sms,sms])
```

### Imputation:

For the right censored states, probabilities weighted substitution costs were used:
\[ d_{inf}(i,j) = \sum\limits_{t=1}^{t_{max}} (Pr(i)_t Pr(j)_t)^T \times  SC \] where \(Pr(i)_t\) and \(Pr(j)_t\) are the vectors of the probabilities of sequence \(i\) and \(j\) being in each state of the alphabet a time \(t\). \(^T\) indicate the outer product and \(SC\) is the substition cost matrix.

The time dependent probability vectors \(Pr(i)_t\) are estimated assuming a first-order markov property \(Pr(X_{t+1}=s) = Pr(X_{t+1}=s|X_t=s_t)\) for all states \(s_1, s_2, ... , s_t, s\).

Software to perform these computations is available in `diagtraject::mis.cost()`.

``` r
ms <- mis.cost(seq_mis, last, tr_t, sm = sub.cost_jacc, 
               cens.type = "right", imp.length = "max", 
               diag=F, sum_to_1=T, resol.comp = resol.comp, 
               resol.ratio = resol.reduc, mc.cores=27)
```

``` r
# Overall dissimilarities
dist_OM <- d_OM + as.matrix(ms$dist)
```

Multidimensional scaling
------------------------

Metric multidimensional scaling using `vegan::wcmdscale()`.

``` r
# finding unique sequences:
mcor <- match(seqconc(seq_mis_last),seqconc(unique(seq_mis_last))) 
uni <- (!duplicated(mcor))
dist <- list("dist"=dist_OM[uni,uni],"weight"= table(mcor),"mcor"= mcor)
```

``` r
# vegan::wcmdscale
wmd <- wcmdscale(dist$dist,w = table(dist$mcor), eig = T, k=13)
```

``` r
R2_list <- lapply(1:13,function(x) {
NewDists  <- dist(wmd$points[mcor,1:x], diag=TRUE, upper=TRUE)
r <- cor(c(as.dist(dist$dist[mcor, mcor], 
                   diag=TRUE, upper=TRUE)), c(NewDists))
r^2})

R2 <- data.frame(k=1:13,R2=do.call('c', R2_list) )       
ggplot(data=R2, aes(x=k,y=R2)) + 
  geom_line(color="blue") + 
  geom_point() +
  scale_y_continuous(limits =c(0,1))+
  labs(title="", x= "Number of Dimensions",y=
         expression(paste("Cumulative",R^{2})))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
```


``` r
wmd_all <-wcmdscale(dist_unique_OM$dist,w = table(mcor), eig = T)
par(mfrow=c(3,4))
for(i in 2:13)
stressplot(wmd_all, k=i, p.col="blue", l.col="red", lwd=2)
```

### Bootstrap Stability of MDS:

Bootstrap stability of MDS over \(n=100\) bootstrap samples of the dissimilarity matrix. Where stability is measured as the stability coefficient:

\[ ST = 1- \frac{\sum_{i=1}^{n}  \parallel X_{i}^{*}-\bar{X}^{*} \parallel^2}{\sum_{i=1}^{n}  \parallel X_{i}^{*} \parallel ^2} \]

where \(X_{i}^{*}\) is the bootstrap MDS solutions after a Procrustes transformation and \(\bar{X}^{*}\) is the average result across all bootstap solutions.

It can be interperpreted as the ratio of between and total variance (See [de Leeuw](https://link.springer.com/content/pdf/10.1007%2FBF01896814.pdf))

It was computed using `diagtraject::bootmds()`, highly inspired by the implementation in `smacof::bootmds`.

``` r
mdsboot <- diagtraject::bootmds(dist$dist,k = 3, nrep = 100, w = table(dist$mcor))
```

``` r
mdsboot
## [1] 0.9986071
```

``` r
dt_mds <- cbind(dt, wmd[[6]]$points[mcor,])
dt_mds[,n_dia := dt.m[,.N,pid]$N-1]
```

``` r
ggplot(dt_mds, aes(x=Dim1, y=Dim2, color=factor(n_dia)))+  
  geom_jitter(h=.05,w=.05)+ 
  scale_colour_brewer(palette = "Spectral")+
  labs(x="Dimension 1", y="Dimension 2",
       color= "Number of comorbid diagnoses")
```


``` r
first_dia <- melt(dt_mds, id.vars = 'pid', 
                  direction = "long", measure.vars =  list(c(3:12)),
                  value.name = "time")
setkey(first_dia,pid)
setkey(dt_mds,pid)
dt_mds[,age_min := first_dia[,min(time,na.rm=T),pid]$V1]
```

``` r
## Correlation with sex and year of birth 


dt_mds[,dia_year:=as.numeric(format(SCZ*365+birthdate,"%Y"))]

newdata=data.frame(birthdate=dt_mds[,birthdate], 
                   age_min=rep(15,nrow(dt_mds)),
                   gender=rep("M",nrow(dt_mds)))

lm <- lm(Dim1~birthdate, data=dt_mds)
conf_interval <- predict(lm, newdata=newdata, 
                         interval="confidence",level = 0.95)
dt_mds[,fit:=conf_interval[,1]]
dt_mds[,lwr:=conf_interval[,2]]
dt_mds[,upr:=conf_interval[,3]]

lm <- lm(Dim1~age_min+gender+birthdate, data=dt_mds)
conf_interval <- predict(lm, newdata=newdata, 
                         interval="confidence", level = 0.95)
dt_mds[,fit1:=conf_interval[,1]]
dt_mds[,lwr1:=conf_interval[,2]]
dt_mds[,upr1:=conf_interval[,3]]

p<- ggplot(dt_mds, aes(y=Dim1, x=birthdate))+
  geom_line(aes(x =birthdate, y=fit ),color="red")+
  geom_line(aes(x =birthdate, y=upr ),color="red", linetype="dashed")+
  geom_line(aes(x =birthdate, y=lwr ),color="red", linetype="dashed")+
  geom_line(aes(x =birthdate, y=fit1 ),color="blue")+
  geom_line(aes(x =birthdate, y=upr1 ),color="blue", linetype="dashed")+
  geom_line(aes(x =birthdate, y=lwr1 ),color="blue", linetype="dashed")+
  geom_jitter(h=.05,w=.5,size=.1)+labs(y="Dimension 1", x="Year of Birth")+
  annotate("text",x=as.Date("1998",format = "%Y"), 
           y=-1, parse=T,color= "blue",
           label = as.character(expression(
             Dim%~%birthdate+onset_age+sex)))+
  annotate("text",x=as.Date("1998",format = "%Y"), 
           y=-5, parse=T,color= "red",
           label = as.character(expression(Dim%~%birthdate)))


knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
To asses the association of resulting MDS coordinates with covariates, we tested for association with birthdate, and found a clear assoication with Dimension 1, however, after adjustment for age of diagnosis, the association disappears. We interpret this as a change in the pattern of diagnoses over time.",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})
p
```


To asses the association of resulting MDS coordinates with covariates, we tested for association with birthdate, and found a clear assoication with Dimension 1, however, after adjustment for age of diagnosis, the association disappears. We interpret this as a change in the pattern of diagnoses over time.
\end{centering}

\end{figure}
``` r
knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
Dimensions show associations with sex, which is included as a covariate in downstream analyses.",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})

ggplot(dt_mds, aes(x=Dim1, y=Dim2, 
                   color=factor(gender,labels = c("female","male"))))+
  geom_jitter(h=.05,w=.05)+ 
  scale_colour_manual(values = c("red","blue"))+
  labs(x="Dimension 1", y="Dimension 2",color= "Sex")
```

``` r
### Dimension Characteristics:
# To visualize how different variables affected the first three Dimensions. We created histograms where 
# all indivuduals are ordered by their coordinate of that dimension. We then create 50 bins and plot the 
# comorbidity pattern for each bin: 

dt_mds[,paste0("Dim",1:7,"r") := lapply(list(
  Dim1,Dim2,Dim3,Dim4,Dim5,Dim6,Dim7), rank,ties.method="random")]
dt_mds[,paste0("Dimension",1:7,"_f") := lapply(list(
  Dim1r,Dim2r,Dim3r,Dim4r,Dim5r,Dim6r,Dim7r), 
  function(x) cut(x, breaks=50))]

dim <- paste0("Dimension",1:7,"_f")

p_10 <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,':=' (count=.N), by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=cut(F10,breaks = seq(0,36,6))]
  t1 <- t1[,.(sum(!is.na(F10)),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
  setkey(t1,xv)
  setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3, fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_brewer(palette = "YlGn",drop=F)+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Substance Abuse ", x=sub("_f","",x),y="%")+ 
    theme(legend.position="none",
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))  
})

p_30 <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,':=' (count=.N), by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=cut(F30,breaks = seq(0,36,6))]
  t1 <- t1[,.(sum(!is.na(F30)),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
  setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3, fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_brewer(palette = "YlGn",drop=F)+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Mood Disorders", x= sub("_f","",x),y="%")+ 
    theme(legend.position="none",
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  
})

p_40 <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,':=' (count=.N), by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=cut(F40,breaks = seq(0,36,6))]
  t1 <- t1[,.(sum(!is.na(F40)),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
  setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3, fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_brewer(palette = "YlGn",drop=F)+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Neurotic Disorders", x= sub("_f","",x),y="%")+ 
    theme(legend.position="none",
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  
})

p_50 <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,':=' (count=.N), by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=cut(F50,breaks = seq(0,36,6))]
  t1 <- t1[,.(sum(!is.na(F50)),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
  setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3, fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_brewer(palette = "YlGn",drop=F)+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Eating Disorders", x= sub("_f","",x),y="%")+ 
    theme(legend.position="none",
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
})

p_60 <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,':=' (count=.N), by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=cut(F60,breaks = seq(0,36,6))]
  t1 <- t1[,.(sum(!is.na(F60)),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
  setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3, fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_brewer(palette = "YlGn",drop=F)+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Personality Disorders", x= sub("_f","",x),y="%")+ 
    theme(legend.position="none",
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  
})

p_70 <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,':=' (count=.N), by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=cut(F70,breaks = seq(0,36,6))]
  t1 <- t1[,.(sum(!is.na(F70)),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
  setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3, fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_brewer(palette = "YlGn",drop=F)+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Mental Retardation", x= sub("_f","",x),y="%")+
    theme(legend.position="none",
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  
})

p_84 <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,':=' (count=.N), by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=cut(F84,breaks = seq(0,36,6))]
  t1 <- t1[,.(sum(!is.na(F84)),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
  setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3, fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_brewer(palette = "YlGn",drop=F)+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Autism Spectrum Disorders", x= sub("_f","",x),y="%")+
    theme(legend.position="none",
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  
})

p_90 <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,':=' (count=.N), by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=cut(F90,breaks = seq(0,36,6))]
  t1 <- t1[,.(sum(!is.na(F90)),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
  setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3, fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_brewer(palette = "YlGn",drop=F)+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Childhood Disorders", x= sub("_f","",x),y="%", 
         fill= "Age at onset in years")+ 
    theme(#legend.position="none",
      legend.direction="vertical",
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))  +
    guides(fill=guide_legend(ncol=2))
})


p_n_dia  <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,count := .N, by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=factor(n_dia,levels=(c(1:7,0)))]
  t1 <- t1[,.(length(n_dia),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
#   t1 <- t1[, n_dia, x2]
#   t1[, V3:=xv]
 setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3,fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_manual(values=c(brewer.pal(7,"BuPu"),"white"))+
    scale_y_continuous(limits =c(0,101))+
    labs(title="Number of Diagnoses", x= sub("_f","",x),y="%", fill="Number of diagnoses")+ 
    theme(legend.direction="vertical",
          #legend.title=element_blank(),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
  guides(fill=guide_legend(ncol=2)) })


p_SCZ <- lapply(dim, function(x) {
  x2 <- eval(noquote(paste0(x)))
  t1 <- dt_mds
  t1[,':=' (count=.N), by=x2]
  t1[, "xv":=t1[,x2, with=F]]
  t1[, cut:=cut(SCZ,breaks = seq(0,36,6))]
  t1 <- t1[,.(sum(!is.na(SCZ)),mean(count)), by=.(xv,cut)]
  t1[, V3:=V1/V2*100]
  setkey(t1,xv)
  setkey(t1,cut)
  ggplot(data=t1, aes(x=xv, y=V3, fill=cut)) +
    geom_bar(stat="identity")+ 
    scale_fill_brewer(palette = "YlGn",drop=F)+
    scale_y_continuous(limits =c(0,101))+
    labs(title="Schizophrenia", x=sub("_f","",x),y="%")+ 
    theme(legend.position="none",
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))  
})

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

leg <- g_legend(p_n_dia[[1]])
leg2 <-g_legend(p_90[[1]]) 
p_n <- lapply(p_n_dia, function(x) x +theme(legend.position="none"))
p_90 <- lapply(p_90, function(x) x +theme(legend.position="none"))

emp <-textGrob("")

# All: 
plist_all<- c(p_n,list(emp),p_10,list(leg),p_30,list(emp),
              p_40,list(leg2),p_50,list(emp),p_60,list(emp),
              p_70,list(emp),p_84,list(emp),p_90,list(emp), p_SCZ, list(emp))

leg_comb <- arrangeGrob(leg,leg2, ncol=1)

#Dim1
plist1<- list(p_n[[1]], p_SCZ[[1]],leg_comb , p_30[[1]],p_60[[1]],p_90[[1]])
#Dim2
plist2<- list(p_n[[2]], p_SCZ[[1]],leg_comb,p_30[[2]],p_84[[2]],p_90[[2]])
#Dim3
plist3<- list(p_n[[3]], p_SCZ[[1]], leg_comb,p_10[[3]],p_30[[3]],p_84[[3]])





do.call("grid.arrange", c(plist1, ncol=3))
```


``` r
do.call("grid.arrange", c(plist2, ncol=3))
```
``` r
do.call("grid.arrange", c(plist3, ncol=3))
```

``` r
do.call("grid.arrange", c(plist_all, ncol=8))
```
``` r
# Barplot of n_dia
#ggplot(data = dt_mds[,mean((Dim1-min(Dim1))/diff(range(Dim1))),n_dia],aes(x= n_dia, y = V1))+geom_bar
#(stat = "identity")


# Regressions:

# n_dia ~ Dim1:
s1 <- summary(glm(formula = n_dia~I(scale(-Dim1)),family = "poisson",data=dt_mds))
c1 <- exp(confint(glm(formula = n_dia~I(scale(-Dim1)),family = "poisson",data=dt_mds)))
s2 <- summary(glm(formula = n_dia~I(scale(-Dim1)),family = "gaussian",data=dt_mds))
c2 <- confint(glm(formula = n_dia~I(scale(-Dim1)),family = "gaussian",data=dt_mds))

qplot(y=n_dia,x=scale(-Dim1),data=dt_mds[n_dia<7])+
   #theme_bw()+
  geom_smooth(aes(y=n_dia,x=scale(-Dim1),color="lm"),method="lm")+
  geom_smooth(aes(y=n_dia,x=scale(-Dim1),color="poisson"),method = "glm", formula = y~x,  method.args = list(family = "poisson"))+
  annotate("text",label=paste0("exp(beta)= ",  round(exp(s1$coefficients[2,1]),2),
                              " [",  round(c1[2,1],2), "-",round(c1[2,2],2), "]",
                              ", AIC=", round(s1$aic,0) ),x=1.5,y=8,col="#00BFC4",size=3)+
    annotate("text",label=paste("beta=",  round(s2$coefficients[2,1],2),
                              " [",  round(c2[2,1],2), "-",round(c2[2,2],2), "]",
                              ", AIC=", round(s2$aic,0) ),x=2.3,y=1,col="#F8766D",size=3)+
  labs(x="Dim1", y="Number of diagnoses",color="Model")
```

``` r
# Child_dia ~ Dim2:
  
s1 <- summary(glm(formula = !is.na(F90)~I(scale(Dim2)),family = "binomial",data=dt_mds))
c1 <- exp(confint(glm(formula = !is.na(F90)~I(scale(Dim2)),family = "binomial",data=dt_mds)))
s2 <- summary(glm(formula = !is.na(F30)~I(scale(Dim2)),family = "binomial",data=dt_mds))
c2 <- exp(confint(glm(formula = !is.na(F30)~I(scale(Dim2)),family = "binomial",data=dt_mds)))

ggplot(aes(fill=!is.na(F90),x=scale(Dim2)),data=dt_mds)+geom_density()+theme_bw()+
  labs(fill="Childhood disorder",x="Dim2") +
  annotate("text",label=paste0("exp(beta)= ",round(exp(s1$coefficients[2,1]),2),
                              " [",  round(c1[2,1],2), "-",round(c1[2,2],2), "]",
                              ", AIC=", round(s1$aic,0) ),x=2.5,y=.5,size=5)
```

``` r
ggplot(aes(fill=!is.na(F30),x=scale(Dim2)),data=dt_mds)+geom_density()+theme_bw()+labs(fill="Mood disorder",x="Dim2") +
  annotate("text",label=paste0("exp(beta)= ",round(exp(s2$coefficients[2,1]),3),
                              " [",  round(c2[2,1],3), "-",round(c2[2,2],3), "]",
                              ", AIC=", round(s2$aic,0) ),x=2.5,y=.5,size=5)
```

``` r
s1 <- summary(glm(formula = !is.na(F10)~I(scale(Dim3)),family = "binomial",data=dt_mds))
c1 <- exp(confint(glm(formula = !is.na(F10)~I(scale(Dim3)),family = "binomial",data=dt_mds)))



ggplot(aes(fill=!is.na(F10),x=scale(Dim3)),data=dt_mds)+geom_density()+theme_bw()+labs(fill="Substance Abuse",x="Dim3") +
  annotate("text",label=paste0("exp(beta)= ",round(exp(s1$coefficients[2,1]),5),
                              " [",  round(c1[2,1],5), "-",round(c1[2,2],5), "]",
                              ", AIC=", round(s2$aic,0) ),x=.5,y=1,size=3)
```

Clustering
----------

``` r
# Hierarchical clustering based on Sequence Analysis Results.
clust <- hclust(dist(dt_mds[,.(Dim1,Dim2,Dim3)]),method="ward.D2")
range <- as.clustrange(clust, diss = dist(dt_mds[,.(Dim1,Dim2,Dim3)]),ncluster = 8,stat="R2")
qplot(2:8,range$stats$R2)+geom_line()+labs(y="R2",x="Clusters")
```

``` r
dt3 <- dt[!is.na(case2012)]
dt3[,cl5:= factor(cutree(clust,5))]
dt3 <- merge(dt_mds, dt3[,.(pid,cl5)],by="pid")

cl5o <- order(tapply(dt3$Dim1,dt3$cl5,mean))
dt3[,cl5:= factor(cutree(clust,5),levels = cl5o, label = 1:length(cl5o))]  


myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(dt3[,cl5])
colScale <- scale_colour_manual(name = "Cluster",values = myColors)
fillScale <- scale_fill_manual(name = "Cluster",values = myColors)


p1 <-ggplot(dt3, aes(x=Dim1, y=Dim2, group=cl5,color=cl5,fill=cl5))+
 stat_density2d(aes(alpha=..level..),geom='polygon',contour = T,show.legend = F) +  
  geom_jitter(h=.05,w=.05)+colScale+fillScale+theme(legend.text=element_text(size=12))
p2 <- ggplot(dt3, aes(x=Dim3, y=Dim1, group=cl5,fill=cl5,color=cl5))+
 stat_density2d(aes(alpha=..level..),geom='polygon',contour = T) +  
  geom_jitter(h=.05,w=.05)+colScale+fillScale

t <- dt3[,.(mean(n_dia),mean(!is.na(F84)|!is.na(F90)),mean(!is.na(F10)),mean(!is.na(F30))),cl5]
t$cl <- factor(t$cl5,levels=c(1:5))
p3 <- ggplot(t, aes(x=cl, y=V1,fill=cl))+
  geom_bar(stat="identity")+fillScale+labs(y="N diagnoses",x="Cluster")

t$cl <- factor(t$cl5,levels=c(1:5))
p4 <- ggplot(t, aes(x=cl, y=V2,fill=cl))+geom_bar(stat="identity")+fillScale+
  labs(y="Childhood Disorders",x="Cluster")+lims(y=c(0,1))

t$cl <- factor(t$cl5,levels=c(1:5))
p5 <- ggplot(t, aes(x=cl, y=V3,fill=cl))+geom_bar(stat="identity")+fillScale+
  labs(y="Substance Abuse, prop.",x="Cluster")+lims(y=c(0,1))

t$cl <- factor(t$cl5,levels=c(1:5))
p6 <- ggplot(t, aes(x=cl, y=V4,fill=cl))+geom_bar(stat="identity")+fillScale+
  labs(y="Affective Disorder, prop.",x="Cluster")+lims(y=c(0,1))

leg <- g_legend(p1)

p_list<- list(p1,p2,p3,p4,p5,p6)

p_list <- lapply(p_list, function(x) x +theme(legend.position="none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")))

lay1 <- rbind(c(1,1:3,3))
lay2 <- rbind(1:4)
g1 <-arrangeGrob(p_list[[1]],leg,p_list[[2]],layout_matrix=lay1)
g2 <-arrangeGrob(p_list[[3]],p_list[[4]],p_list[[5]],p_list[[6]],layout_matrix=lay2)
grid.arrange(g1,g2, top="Results of  k=5 clustering")
```

### Cluster stability assessment

Clustering was repeated with 100 bootstap permutations of the MDS-dimensions. The overlap between the clustering of the full dataset and the bootstrap permutations was assessed using the mean jaccard coefficient.

``` r
boot<-clusterbootw(data = dist_OM,distances = T,
                   clustermethod=mdshclustCBI, k= 5, k.mds=3,
                   B=100, method="ward.D2", members=NULL, mc.cores = 13)
kable(data.frame(N=table(boot$partition),Jaccard=boot$bootmean,Recovered=boot$bootrecover, 
                 Dissolved= boot$bootbrd),format="latex",booktabs=T)
```

``` r
# Visualizing stability using `WGCNA`:
labels <- matrix(NA, ncol=length(boot$partition), nrow=101)

labels[1,] <-   boot$partition
for(i in 2:101) labels[i,] <- matchLabels(boot$bootpartition[,(i-1)],
                                          boot$partition)

colors <- labels2colors(labels,colorSeq= brewer.pal(7,"Set1"))


knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
Above: dendrogram of the hierarchical agglomerative clustering. Below: Results of clusterings repeated on 100 random subsets of the data. For each permutation k=3-MDS was computed, and subsequently k=5-clustering was perfomed.",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})

plotDendroAndColors(boot$result$result,t(colors), main="",
                    c("", unlist(lapply(1:10, function(x) c(rep("",9),paste("Resampl.",x*10))))),
                    dendroLabels =F, hang=.03,autoColorHeight = T)
```

``` r
# Load additional genetic and registry variables
load("additional.Rda")
```

``` r
source("~/trajectory/single/all/script/Nagelkerke.R")

dt_ld2<- merge(dt_ld, dt[,.(pid,case,birthdate,gender,random)], by="pid")
pval <- lapply(61:67, function(x) {
  tmp = dt_ld2 
  tmp$y = tmp[,x,with=F]
  summary(glm(!is.na(case)~y+gender+
                C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave, 
              data=tmp, family="binomial" ))[[13]][2,4]})

null <- glm(!is.na(case)~birthdate+gender+
              C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave, 
            data=dt_ld2, family="binomial" )

nag <-  lapply(61:67, function(x) {
  tmp = dt_ld2 
  tmp$y = tmp[,x,with=F]
  nagelkerke(fit = glm(!is.na(case)~y+birthdate+gender+
                         C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave, 
                       data=tmp, family="binomial" ),
             null = null)$Pseudo.R.squared.for.model.vs.null[3,]})
```

``` r
options(digits=3, scipen=T)
df_prs <-data.frame(Phenotype= name, p= unlist(pval), R2= unlist(nag))
names(df_prs)[2:3] <- paste0(names(df_prs)[2:3],footnote_marker_symbol(1)) 
df_prs %>%
kable(digits=c(1,50,4),format="latex",booktabs=T, caption= "Association with Schizophrenia",escape=F )%>%
row_spec(which(df_prs$p<0.0001), bold = T) %>% footnote(general=c(
"Polygenic Scores were computed for Chromosomes 1-22 using LDpred. Imputation was done using the",
" Hapmap3 (n=1457897) reference panel. All SNPs with an info Score >0.8 and MAF >0.05 were included",
"(No p-valusthresholding). LDradius was set to 200 (default in LDpred) PseudoR2 is calculated using",
" Nagelkerke's - comparing the full model to the model without the PGS.", 
paste0("P-values printed in bold indicates p","$<$","0.0001)."),
paste0("\\\\", footnote_marker_symbol(1), ": Based on a binomial logistic regression of 2861 cases and 18843"), "controls,adjusting for age, sex, 10 principal components of genetic similarity and genotype wave."),escape=F)
```

``` r
# Merge MDS with risk variables:

data <- merge(dt_mds[,1:33,with=F],additional,by="pid",all.x = F)
dim(data)
```

    ## [1] 5432   99

``` r
ind_birth <- c(34:41,86:87) 
ind_mat_infek <- 42:43
ind_prs <- (93:99)[pval<.0001]
ind_pc_wave <- 56:66
ind_sev <- 67:73
ind_infek <- 74:77
ind_rare <- 80:81
ind_parental <- 88:91
```

### Medical Birth Registry (MBR)

Online [documentation](https://sundhedsdatastyrelsen.dk/-/media/sds/filer/registre-og-services/nationale-sundhedsregistre/graviditet-foedsler-og-boern/foedselsregisteret/dokumentation-af-fdselsregisteret-1973-1996.xls?la=da).

``` r
data$apgar5=factor(data$apgar5,
                   levels=c("0","1","2","3","4","5",
                            "6","7","8","9","10","A"))
data$apgar5=as.character(data$apgar5)
data$apgar5[data$apgar5=="A"]=NA
data$apgar5=as.numeric(data$apgar5)

kable(t(data.frame("N"=table(data$apgar5)[table(data$apgar5)>4])),
      caption = "Distribution of APGAR Scores in cohort -omitting counts <5",
      format="latex",booktabs=T) %>% footnote(general=
"APGAR 5: Appearance, Pulse, Grimace, Activity, Respiration at 5 minutes")
```

``` r
#names(table(data$"Birth Length")[table(data$"Birth Length")>0])
leng=data$"Birth Length"
leng[leng %in%c("A","1","10")]=NA
leng=as.numeric(as.character(leng))
data$"Birth Length"=leng;rm(leng)
rare <- table(data$"Birth Length")[table(data$"Birth Length")<4]
ggplot(data[!`Birth Length` %in%as.numeric(names(rare))], aes(`Birth Length`))+
  geom_histogram(binwidth=2,fill="steelblue")+labs(title="Birth Length (N= 5315)")
```


``` r
regular_plot <- knit_hooks$get("plot")

# Custom knitr hook to add notes to the plot
knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
When calculating the value of birth weight score all singletons in the MBR born between 1981 and 2005 and 
having information of birth weight as well as gestational age were divided according to sex,
gestational age in weeks, and for persons with gestational age of 28 weeks or more also into the
following groups of calendar year at date of birth: 1981-1985, 1986-1990, 1991-1995, 1996-2000, 
2001-2005. The variable birth weight score contains the proportion of persons with same sex, gestational age,
and calendar period (only if gestational age is 28 weeks or larger) who have the same or a smaller
birth weight than the index person.",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})
data$`Birth Weight`=data$`Birth Weight`/100

ggplot(data, aes(`Birth Weight`))+geom_histogram(binwidth=.20,fill="steelblue")+labs(x="Birth Weight Score")+xlim(0,1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
```

``` r
data$C_RYGER[data$C_RYGER==99]=NA
data$"Maternal Smoking in Pregnancy" <- data$B_RYGER
data$"Maternal Smoking in Pregnancy"[!is.na(data$C_RYGER)] <-1
kable(t(data.frame("N"=data[,.N,`Maternal Smoking in Pregnancy`][,N],row.names = c("Unknown","No","Yes"))),
      caption ="Maternal Smoking during Pregnancy",format="latex",booktabs=T,align="l") %>%
  column_spec(1:2,width="20em") %>% footnote(general=c(
"Maternal smoking was not in the register before 1991.", 
"From 1991 to 1996 it was registered as a binary (smoker of non-smoker).", 
"In the period 1997 to 2005 the register had more refined data, but this was transformed in a binary",
"(non-smoker=non-smoker. smoker= smoker, stopped smoking during the first trimester, stopped smoking",
"after the first trimester, smokes at most 5 cigarettes daily, smokes 6-10 cigarettes daily, smokes",
"11-20 cigarettes daily, smokes more than 20 cigarettes daily or smoker quantity unspecified).",
"In the analysis smoking during pregnacy is treated as a binary."))
```

``` r

data$wave[is.na(data$wave)]="Other"
data <- data[,wave:=factor(wave)]

ind_birth <- c(ind_birth[-c(7:8)],100) 
```

``` r

data[,Sectio:= B_I11]
data[is.na(Sectio), Sectio := B_SECTIOU]

kable(t(data.frame("N"=data[,.N,`Sectio`][,N],row.names = c("No","Yes","Unknown"))),
      caption ="Sectio according to the Medical Birth Registry",format="latex",booktabs=T,align="l") %>%
  column_spec(1:2,width="20em") %>% footnote(general=c( 
"..."))
```

``` r

data$wave[is.na(data$wave)]="Other"
data <- data[,wave:=factor(wave)]

ind_birth <- c(ind_birth[-c(7:8)],101)
```

``` r
nms = grep("d2100",names(data))
nms <- names(data)[nms];nms
## [1] "d2100_ptype0_contacts" "d2100_ptype0_days"     "d2100_ptype1_contacts"
## [4] "d2100_ptype1_visits"   "d2100_ptype2_visits"   "d2100_ptype2_novisits"
## [7] "d2100_ptype3_contacts"
data <- data[,"Number of hospitalizations (SCZ) total":= 
               Reduce('+',.SD), .SDcols=nms[c(1,3,7)]]
data <- data[,"Number of hospitalizations (SCZ)":= `Number of hospitalizations (SCZ) total`/
                (as.numeric(as.Date("2016-12-31")-birthdate)/365.25 -  SCZ )]

rare <- table(data$"Number of hospitalizations (SCZ)")[table(data$"Number of hospitalizations (SCZ)")<4]

knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
All psychiatric hospital contacts with Schizophrenia as either main diagnosis
(aktionsdiagnose/hoveddiagnose) or basic diagnosis (grundmorbus) in the Danish Psychiatric Central
Research Register complete until December 31, 2016.",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})
grid.arrange(
ggplot(data[!`Number of hospitalizations (SCZ)` %in%as.numeric(names(rare))], aes(`Number of hospitalizations (SCZ)`))+
  geom_histogram(binwidth=1,fill="steelblue")+labs(x="Number of hospitalizations (SCZ) pr year")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")),
ggplot(data[!`Number of hospitalizations (SCZ)` %in%as.numeric(names(rare))]
       , aes(log10(`Number of hospitalizations (SCZ)`+1)))+
  geom_histogram(bins=20,fill="steelblue")+labs(x="log10(Number of hospitalizations (SCZ) pr year +1)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")))
```

``` r

data$"Number of hospitalizations (SCZ)" <- log10(data$`Number of hospitalizations (SCZ)`+1)
```

``` r
data <- data[,"Total time hospitalized (SCZ) total" := d2100_ptype0_days]


data <- data[,"Total time hospitalized (SCZ)":= `Total time hospitalized (SCZ) total`/
                (as.numeric(as.Date("2016-12-31")-birthdate)/365.25 -  SCZ )]


rare <-table(data$"Total time hospitalized (SCZ)")[table(
  data$"Total time hospitalized (SCZ)")<4]


knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
Contains the total number of days admitted as 24 hour inpatient with a Schizophrenia diagnosis.
For an admission date of the admission (unfinished admission)
the admission is presumed ongoing until December 31, 2016. An exception for this rule is when
the person has another contact in the Danish Psychiatric Central Research Register starting after
the first unfinished admission (ptype = 0). In this case all unfinished admissions with ptype = 0
for this person are counted as having a duration of 1 day only.
For persons with overlapping admissions with patient type 0 each day only counts once. ",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})
grid.arrange(
ggplot(data[!`Total time hospitalized (SCZ)` %in%as.numeric(names(rare))], 
       aes(`Total time hospitalized (SCZ)`))+
  geom_histogram(binwidth=10,fill="steelblue")+labs(x="Total time hospitalized (SCZ), days pr year")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")), 
ggplot(data[!`Total time hospitalized (SCZ)` %in%as.numeric(names(rare))], aes(log10(`Total time hospitalized (SCZ)`+1)))+
  geom_histogram(binwidth=.1,fill="steelblue")+labs(x="Total time hospitalized (SCZ) log10(days pr year+1)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")))
```


``` r
data$`Total time hospitalized (SCZ)` <- log10(data$`Total time hospitalized (SCZ)`+1)
ind_sev <- c(103,105); names(data)[ind_sev]
## [1] "Number of hospitalizations (SCZ)" "Total time hospitalized (SCZ)"
```

``` r

kable(
data.frame("Viral"= 
data[,.N,`Maternal Infection during Pregnancy, Viral`][,N], 
           "Bacterial"=
data[,.N,`Maternal Infection during Pregnancy, Bacterial`][,N], 
           row.names = c("no","yes")),caption = 
  "Maternal Infection during Pregnancy",format="latex",booktabs=T,align="l") %>%
  column_spec(1,width="20em") %>% footnote(general=c(
"From the National Patient Register, information an maternal infections during pregnancy was obtained.",
"It was defined as a maternal diagnosis of infection (bacterial and viral) within a period of time",
"prior to the date of birth that corresponded to the gestational age of the child."))
```


``` r
kable(
data.frame("Viral"= data[,.N,`Viral Infection`][,N], 
           "Bacterial"=data[,.N,`Bacterial Infection`][,N],
           "CNS"= data[,.N,`CNS Infection`][,N], 
           "Otitis"=data[,.N,`Otitis Infection`][,N],
           row.names = c("no","yes")),caption = "Infections",format="latex",booktabs=T,align="l") %>%
  column_spec(1:3,width="5em")  %>% 
  footnote(general=
             "National Patient Register was used to to identify patients diagnosed with infections")
```


``` r
kable( 
data.frame("Schizophrenia"= 
             data[,.N,`Paternal Diagnosis, Schizophrenia`][,N], 
           "Any Psychiatric"= 
             data[,.N,`Paternal Diagnosis, Any Psychiatric`][,N], 
           "Schizophrenia "= 
             data[,.N,`Maternal Diagnosis, Schizophrenia`][,N], 
           "Any Psychiatric "= 
             data[,.N,`Maternal Diagnosis, Any Psychiatric`][c(2,1,3),N], 
           row.names = c("no","yes","Unknown")),caption = 
  "Parental history of mental disorders",format="latex",booktabs=T) %>%
  add_header_above(c(" "=1,"Paternal"=2,"Maternal"=2))%>% 
  footnote(general=c( 
"National Patient Register was used to to identify patients",
"with parental history of mental disorders"))
```


``` r
### PGS distributions


prs = names(data)[ind_prs]  
data[!is.na(C1),c(prs):= lapply(.SD,scale),.SDcols=prs]

p <- lapply(ind_prs, function(x){
 tmp<- data 
 tmp$x1 <- tmp[,x,with =F]
 ggplot(tmp, aes(x1))+geom_histogram(fill="steelblue")+labs(x=names(data)[x])})
do.call('grid.arrange',p)
```


``` r
### Association Analysis


# Wrapper functions for `stats::manova` and `stats::anova`:


mancova.func <- function(tmp,ind_x,ind_y,mf){
tmp<- as.data.frame(tmp)
p.val <- Fstat <- N<-  rep(NA,length(ind_x))
for(p in (1:length(ind_x))){
  tmp$x=tmp[,ind_x[p]]
  tmp$y = as.matrix(tmp[,ind_y])
  manova.y=manova(mf, data = tmp[!is.na(tmp$x),]);
  N[p] = sum(!is.na(tmp[,ind_x[p]]))
  p.val[p] = summary(manova.y)$stats["x",6]   
  Fstat[p] = summary(manova.y)$stats["x",3]
  } 
names(N) <- names(p.val) <- names(Fstat) <-names(tmp)[ind_x]
list(N= N, Fstat=Fstat,p.val=p.val) }

ancova.func <- function(tmp,ind_x,ind_y,mf){
tmp<- as.data.frame(tmp)
p.val <- betas <- array(NA,dim=c(length(ind_y),length(ind_x)))
model <-  rep( list(list()), length(ind_x) )
for(p in (1:length(ind_x))){
  for(q in (1:length(ind_y))){
    tmp$x=tmp[,ind_x[p]]
    tmp$y = as.matrix(tmp[,ind_y[q]])
    lm = lm(mf, data = tmp[!is.na(tmp$x),])
    anova.y=anova(lm)
    p.val[q,p] = anova.y["x",5]
    betas[q,p] = coef(lm)["x"]
    model[[p]][[q]] =lm}} 
colnames(p.val)=names(tmp)[ind_x]
colnames(betas)=names(tmp)[ind_x]
list(p.val=p.val,beta=betas, model= model)}



clreg.func <- function(tmp,ind_x,pc=F,reg=c("lm","lr")){
p.val <- coef <- array(NA,dim=c(length(ind_x),length(levels(tmp$y))-1))
 for(q in (1:length(ind_x))){
    
    tmp$xx = tmp[,ind_x[q],with=F]
    
    if(pc)
    clreg <-  lm(xx ~ birthdate+gender 
             +C1 +C2 +C3 +C4 +C5 +C6 +C7 +C8 +C9 +C10+wave+y,data=tmp) else if(reg=="lm")
               clreg <-  lm(xx ~ birthdate+gender+y,data=tmp) else
                  clreg <-  glm(xx ~ birthdate+gender+y,data=tmp,family="binomial") 
    
    s<- summary(clreg)
    
    #print(p)
    p.val[q,]= s$coefficients[paste0("y",levels(tmp$y)[-1]),4]
    coef[q,]= s$coefficients[paste0("y",levels(tmp$y)[-1]),1]

    #if(reg!="lm") coef[q,] = exp(coef[q,])
    }
rownames(p.val)<-rownames(coef)<-names(tmp)[ind_x]
colnames(p.val)<-colnames(coef)<-paste0("cl",levels(tmp$y)[-1])
return(list(p.val=p.val,coef=coef)) #, model= model
     }
```

``` r
# Number of tests
n_test <- length(c(ind_prs, ind_rare[1],ind_sev,ind_birth,ind_infek,
                          ind_mat_infek, ind_parental))
```

``` r
# Polygenic Scores:

ind_y=grep("Dim",names(data))[1:3]
ind_x=ind_prs;names(data)[ind_x]
## [1] "PGS - Bipolar Affective Disorder" "PGS - Depressive Symptoms"       
## [3] "PGS - Education, Years"           "PGS - Neuroticism"               
## [5] "PGS - Schizophrenia"
mf = formula(y ~ birthdate+gender 
             +C1 +C2 +C3 +C4 +C5 +C6 +C7 +C8 +C9 +C10+wave+x)

manova.results_prs<-mancova.func(data[!is.na(data$C1),], ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_prs$p.val <= .01]
anova.results_prs<-ancova.func(data[!is.na(data$C1),], ind_x_sig,ind_y,mf)
#multinom.results_prs <- multinom.func(merge(data[!is.na(C1)],dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              #, ind_x_sig,pc=T)

clreg.results_prs <- clreg.func(merge(data[!is.na(C1)],dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              , ind_x_sig,pc=T,reg="lm")
```

``` r
# Association with rare mutations:  


ind_x=ind_rare;names(data)[ind_x]
```

    ## [1] "Disruptive or Damaging Mutations" "Synonymous Mutations"

``` r
mf = formula(y ~ birthdate*gender+C1+C2+C3+C4+C5+C6+C7+C8+C8+C9+C10+x)

manova.results_rare<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_rare$p.val <= .05/n_test]
#anova.results_rare<-ancova.func(data, ind_x_sig,ind_y,mf)
```

``` r
# Contacts and visits
ind_x=ind_sev;names(data)[ind_x]
## [1] "Number of hospitalizations (SCZ)" "Total time hospitalized (SCZ)"
mf = formula(y ~ birthdate*gender+x)

manova.results_sev<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_sev$p.val <= .05/n_test]
anova.results_sev<-ancova.func(data, ind_x_sig,ind_y,mf)
#multinom.results_sev <- multinom.func(merge(data,dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              #, ind_x_sig,pc=T)
clreg.results_sev <- clreg.func(merge(data[!is.na(C1)],dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              , ind_x_sig,pc=F,reg="lm")
```

``` r
# Association with birth & pregnancy variables:

ind_x=c(ind_birth);names(data)[ind_x]
## [1] "Maternal Age"                  "Paternal Age"                 
## [3] "Birth Length"                  "apgar5"                       
## [5] "Gestational Age"               "Birth Weight"                 
## [7] "Maternal Smoking in Pregnancy" "Sectio"
mf = formula(y ~ birthdate*gender+x)

manova.results_birth<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_birth$p.val <= .05/n_test]
anova.results_birth<-ancova.func(data, ind_x_sig,ind_y,mf)

#multinom.results_birth <- multinom.func(merge(data,dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              #, ind_x_sig,pc=T)

clreg.results_birth <- clreg.func(merge(data[!is.na(C1)],dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              , ind_x_sig,pc=F,reg="lm")
```

``` r
# Association with infection diagnoses during pregancy  

ind_x=ind_mat_infek;names(data)[ind_x]
## [1] "Maternal Infection during Pregnancy, Bacterial"
## [2] "Maternal Infection during Pregnancy, Viral"
mf = formula(y ~ birthdate*gender+x)

manova.results_mat_infek<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_mat_infek$p.val <= .05/n_test]
#anova.results_mat_infek<-ancova.func(data, ind_x_sig,ind_y,mf)
```

``` r

#Parental diagnosis of Schizophrenia and all Fxx diagnosis 

ind_x=ind_parental;names(data)[ind_x]
## [1] "Maternal Diagnosis, Any Psychiatric"
## [2] "Maternal Diagnosis, Schizophrenia"  
## [3] "Paternal Diagnosis, Any Psychiatric"
## [4] "Paternal Diagnosis, Schizophrenia"
mf = formula(y ~ birthdate*gender+x)

manova.results_fam<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_fam$p.val <= .05/n_test]
anova.results_fam<-ancova.func(data, ind_x_sig,ind_y,mf)
#multinom.results_fam <- multinom.func(merge(data,dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              #, ind_x_sig,pc=T)

clreg.results_fam <- clreg.func(merge(data[!is.na(C1)],dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              , ind_x_sig,pc=F,reg="lr")
```

``` r

# Association with infection diagnoses 
ind_x=ind_infek;names(data)[ind_x]
## [1] "Otitis Infection"    "Bacterial Infection" "CNS Infection"      
## [4] "Viral Infection"
mf = formula(y ~ birthdate*gender+x)

manova.results_infek<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_infek$p.val <= .05/n_test]
anova.results_infek<-ancova.func(data, ind_x_sig,ind_y,mf)
#multinom.results_infek <- multinom.func(merge(data,dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              #, ind_x_sig,pc=T)
clreg.results_infek <- clreg.func(merge(data[!is.na(C1)],dt3[,.(pid,y=relevel(cl5,ref = 5))],by="pid") 
              , ind_x_sig,pc=F,reg="lr")
```

``` r
#Heatmap of associations using modification of `corrplot::corrplot` included in `diagtraject::corrplot`:


anovas <- list(
  anova.results_prs,
  #anova.results_rare,
  anova.results_birth,
  anova.results_infek,
  #anova.results_mat_infek,
  anova.results_fam,
  anova.results_sev)
betas <- t(Reduce('cbind',lapply(anovas,function(x) x$beta)))
p.vals <- t(Reduce('cbind',lapply(anovas,function(x) x$p.val)))

colnames(betas) <- colnames(p.vals) <- paste("Dimension",1:ncol(betas))  

p.vals2 <- p.vals
p.vals2[p.vals<1e-6] <- 1e-6
betas2 <- betas
betas2[betas2< -.5] <- -.49
betas2[betas2>.5] <- .49
col <- colorRampPalette(c("#67001F", "#B2182B", "#FFFFFF", "#2166AC", "#053061"))(20)
corr.plot <- corrplot(betas2[,1:3], is.corr=FALSE, p.mat=p.vals2[,1:3],
                      method="p.square", insig="pch", col=col,
                      sig.level=0.05/nrow(betas), pch="*", pch.cex=1.2, 
                      cl.ratio=0.25, cl.align.text="l", cl.lim=c(-.5,.5), 
                      tl.col="black", bg="#fffffc",ssl.pos = "r",ssl.lev=6,full_col=F)
```

``` r
# Association with genetic and register variables is tested using MANCOVA with MDS dimensions as 
# dependent variables and risk variables as independent variables, adjusting for age, sex and for 
# genetic variable also genotype wave and 10 principal components of genetic ancestry. Bonferoni 
# correction for 27 tests.

Mancova_table <- rbind(
as.data.frame(manova.results_prs),
as.data.frame(manova.results_rare),
as.data.frame(manova.results_birth),
as.data.frame(manova.results_infek),
as.data.frame(manova.results_mat_infek),
as.data.frame(manova.results_fam),
as.data.frame(manova.results_sev)
)


options(digits=2, scipen=T)
# Full table:
kable(Mancova_table
,digits=c(100,1),format="latex",booktabs=T,,caption="MANCOVA Results") %>%  footnote(general =  c(
"MANCOVAs done with first three dimensions of MDS as dependent variable, adjusting for age and sex.",
"For genetic varibles adjustning additionally for 10 principal components and genotype wave."))
```

``` r
#Subset for main text:
tab2 <- cbind(Mancova_table[rownames(Mancova_table) %in% rownames(betas),],
            coef= betas[,1],p= p.vals[,1],coef =betas[,2],p= p.vals[,2],coef=betas[,3],p= p.vals[,3])
rownames(tab2) <- sub(" Diagnosis", "",rownames(tab2))

kable(tab2,
      digits=c(1,1,100,2,100,2,100,2,100),format="latex",booktabs=T,caption="Association tests") %>% 
  kableExtra::group_rows("Genetic variables",1,1) %>% 
  kableExtra::group_rows("National Birth Registry",2,6) %>% 
  kableExtra::group_rows("National Patient Registry",7,8) %>% 
  kableExtra::group_rows("Family History",9,11) %>% 
  kableExtra::group_rows("Outcome variables (PCRR)",12,13) %>%
  add_header_above(c(" "=4,"Dim1"=2,"Dim2"=2,"Dim3"=2)) %>%
  add_header_above(c(" "=2,"MANCOVA results"=2,"Post-hoc regressions"=6)) %>%
  footnote(general =  c( 
"MANCOVAs done with first three dimensions of MDS as dependent variable",
"Post-hoc regressions are linear regressions with the dimension score as dependent variable.",
"Only varibles with significant p-values after adjustment for 27 tests are included.",
"All analyses are adjusting for age and sex PGS-Education is adjusted additionally for 10 principal", "components and genotype wave.",
"PCRR: Psychiatric Cental Research Registry"))
```


Replication
-----------

### Replication dataset:

Sequence Analysis is performed in non-overlapping dataset consisting of patients diagnosed January 1 2013 - December 31 2016.

``` r
load("diagdates2016.Rda")  
diagdates_rep <- data.table(rbind(diagdates,diagdates2016))
```

``` r
# Sequence object is created as above for dataset 
# with replication samples and primary cohort:

dt_rep <- data.table(diagdates_rep)

dt_rep[,"birthdate":= fdate(birthdate) ]
dt_rep[,c("F1","F3", "F4","F50","F60", "F70", "F84", "F9","SCZ") :=
     lapply(list(F1,F3,F4,F50,F60,F70,F84,F9,SCZ),function(x){  
       (as.numeric(fdate(x))- as.numeric(birthdate))/365})]
dt_rep[,censored:= (as.numeric(fdate("31/12/2016"))-as.numeric(birthdate))/365] #Calculates age in 31/12/2016
```

``` r
dt.m <-melt(dt_rep, id.vars = 'pid', direction = "long", 
            measure.vars  = list(c(3:10,12)), value.name = "time",
            variable.name = "event")
dt.m<- dt.m[order(pid)]
dt.m[,time:=time]
dt.m <- dt.m[!is.na(time)]

events <- levels(dt.m$event)
drop <-  matrix(FALSE, nrow = length(events), ncol = length(events),
                dimnames = list(events,events))
drop['censored',] <- T
drop[,'censored'] <- T
diag(drop) <- F

e2sm <-seqe2stm(events = events, dropMatrix = drop)

seq_rep <- seqdef(TSE_to_STS(dt.m,id = "pid", timestamp = "time", 
                             event= "event",tmin=1, 
                             tmax=ceiling(max(dt.m$time)),stm = e2sm),
                  firstState ="None")

seq_mis_rep <- seqdef(TSE_to_STS(dt.m, id = "pid", timestamp = "time",
                                 event= "event",tmin=1, 
                                 tmax=ceiling(max(dt.m$time)),stm = e2sm), 
                      firstState ="None", missing="censored")

drop['censored',] <- F
drop[,'censored'] <- T
e2sm <-seqe2stm(events = events, dropMatrix = drop)

seq_last <- seqdef(TSE_to_STS(dt.m, id = "pid", timestamp = "time", 
                              event= "event", tmin=1,
                              tmax=ceiling(max(dt.m$time))+1,stm = e2sm),
                   firstState ="None", missing="censored")

seq_mis_last <- seq_mis
seq_mis_last$last <- seq_last[,ncol(seq_last)]
attributes(seq_mis_last)$alphabet <- c(alphabet(seq_mis_last),  
                                       unique(seq_mis_last$last)[which(
                                         !unique(seq_mis_last$last) %in% 
                                           alphabet(seq_mis_last))])
attributes(seq_mis_last)$alphabet <- attributes(seq_mis_last)$alphabet[
  order(attributes(seq_mis_last)$alphabet)]
attributes(seq_mis_last)$labels <- attributes(seq_mis_last)$alphabet

 sms <- which(rownames(sub.cost_jacc) %in% rownames(seqsubm(seq,"CONSTANT")))

d_OM <- seqdist(seq_rep, method = "OM", indel=.5, 
                sm = sub.cost_jacc[sms,sms])
ms <- mis.cost(seq_mis_rep, seq_mis_last$last, tr_t, sm = sub.cost_jacc, 
               cens.type = "right", imp.length = "max", 
               diag=F, sum_to_1=T, resol.comp = resol.comp, 
               resol.ratio = resol.reduc, mc.cores=27)

dist_OM_rep <- d_OM + as.matrix(ms$dist)

# finding unique sequences
mcor <- match(seqconc(seq_mis_last),seqconc(unique(seq_mis_last))) 
uni <- (!duplicated(mcor))
dist_unique_OM_rep <- list("dist"=dist_OM_rep[uni,uni],
                           "weight"= table(mcor),"mcor"= mcor)
```

``` r
# Weighted MDS is performed assigning weights of `10^-20` to replication sample.  
dt_rep[,mcor2:=dist_unique_OM_rep$mcor]
case_cor <- dt_rep[case==1,sum(case2012,na.rm = T)+1e-20,mcor2]
setkey(case_cor, mcor2)
w <- case_cor$V1
```

``` r

# Showing that MDS of Case2012 in unaltered 

wmd2012 <- wcmdscale(dist_unique_OM_rep$dist[w>=1,w>=1],w = w[w>=1], k=7) 
wmd_rep <- wcmdscale(dist_unique_OM_rep$dist,w = case_cor$V1, k=7,) 

#test
#wmd2012t<- wcmdscale(dist_OM, k=7) 
#head(wmd_rep[dt_rep$mcor2,][!is.na(dt_rep$case2012),])
#cor(wmd_rep[dt_rep$mcor2,][!is.na(dt_rep$case2012),],wmd2012t)

sum(abs(scale(wmd2012[,1])-scale(wmd_rep[w>=1,1])))
## [1] 6e-13

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

specify_decimal(diag(cor(scale(wmd2012)-scale(wmd_rep[w>=1,]))),15)
## [1] "1.000000000000000" "0.999999999999993" "0.999999999999998"
## [4] "1.000000000000000" "1.000000000000000" "1.000000000000000"
## [7] "1.000000000000000"
```

``` r
dt_rep[,paste0("Dim",1:7):=as.data.frame(wmd_rep[dt_rep$mcor2,])]

dt_rep[,Dim1 := Dim1 *diag(cor(
  dt_rep[!is.na(case2012),"Dim1",with=F],
  data[!is.na(case2012),"Dim1",with=F]))]

dt_rep[,Dim2 := Dim2 *diag(cor(
  dt_rep[!is.na(case2012),"Dim2",with=F],
  data[!is.na(case2012),"Dim2",with=F]))]

dt_rep[,Dim3 := Dim3 *diag(cor(
  dt_rep[!is.na(case2012),"Dim3",with=F],
  data[!is.na(case2012),"Dim3",with=F]))]


### Classifing sequences from replication data using K-nearest Neighboors 

# class::knn

knn <- knn(train = dt3[,.(Dim1,Dim2,Dim3)], test = dt_rep[,.(Dim1,Dim2,Dim3)], cl= dt3[,cl5],k=5)
dt_rep[,Cluster := knn]
dt_rep[,cl5_knn := knn]
dt_rep[,Cohort := factor(is.na(case2012),labels = c("Main","Replication"))]

knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
Projection of Replication Cohort into MDS space of primary cohort and subsequent classification using K-Nearest-Neighbors with k=5.",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})

ggplot(data = dt_rep, aes(x=Dim1,y=Dim2,col=Cluster,size=Cohort))+geom_jitter(h=.05,w=.05)
```

\begin{figure}

\begin{centering}
\caption{Replication Sample MDS}
\includegraphics[width=\maxwidth]{Analysis_description_files/figure-markdown_github/unnamed-chunk-56-1.pdf}

\

Projection of Replication Cohort into MDS space of primary cohort and subsequent classification using K-Nearest-Neighbors with k=5.
\end{centering}

\end{figure}
``` r
data_rep <- merge(dt_rep,additional,by="pid")
```

``` r
# Modified as above:

ind_birth <- c(34:41,86:87) 
ind_mat_infek <- 42:43
ind_prs <- (93:99)[pval<.0001]
ind_pc_wave <- 56:66
ind_sev <- 67:73
ind_infek <- 74:77
ind_rare <- 80:81
ind_parental <- 88:91

table(data_rep$apgar5)[table(data_rep$apgar5)>4]
data_rep$apgar5=factor(data_rep$apgar5,levels=c("0","1","2","3","4","5",
                                                "6","7","8","9","10","A"))
data_rep$apgar5=as.character(data_rep$apgar5)
data_rep$apgar5[data_rep$apgar5=="A"]=NA
data_rep$apgar5=as.numeric(data_rep$apgar5)
table(data_rep$apgar5)[table(data_rep$apgar5)>4]
#names(table(data_rep$"Birth Length")[table(data_rep$"Birth Length")>0])
leng=data_rep$"Birth Length"
leng[leng %in%c("A","1","10")]=NA
leng=as.numeric(as.character(leng))
data_rep$"Birth Length"=leng;rm(leng)
table(data_rep$"Birth Length")[table(data_rep$"Birth Length")>4]

#hist(data_rep$"Birth Weight",na.rm = T)
table(data_rep$B_RYGER)
names(table(data_rep$C_RYGER))
data_rep$C_RYGER[data_rep$C_RYGER==99]=NA
table(data_rep$C_RYGER)[table(data_rep$C_RYGER)>4]
data_rep$"Maternal Smoking in Pregnancy" <- data_rep$B_RYGER
data_rep$"Maternal Smoking in Pregnancy"[!is.na(data_rep$C_RYGER)] <-1
table(data_rep$"Maternal Smoking in Pregnancy")

data_rep$wave[is.na(data_rep$wave)]="Other"
data_rep <- data_rep[,wave:=factor(wave)]

ind_birth <- c(ind_birth[-c(7:8)],100); names(data_rep)[ind_birth]
nms = grep("d2100",names(data_rep))
nms <- names(data_rep)[nms];nms
data_rep <- data_rep[,"Number of hospitalizations (SCZ) total":= 
               Reduce('+',.SD), .SDcols=nms[c(1,3,7)]]
data_rep <- data_rep[,"Number of hospitalizations (SCZ)":= `Number of hospitalizations (SCZ) total`/
                (as.numeric(as.Date("2016-12-31")-birthdate)/365.25 -  SCZ )]
data_rep$"Number of hospitalizations (SCZ)" <- log10(data_rep$`Number of hospitalizations (SCZ)`+1)
table(data_rep$"Number of hospitalizations (SCZ)")[
  table(data_rep$"Number of hospitalizations (SCZ)")>4]
data_rep <- data_rep[,"Total time hospitalized (SCZ) total" := d2100_ptype0_days]
data_rep <- data_rep[,"Total time hospitalized (SCZ)":= `Total time hospitalized (SCZ) total`/
                (as.numeric(as.Date("2016-12-31")-birthdate)/365.25 -  SCZ )]
data_rep$"Total time hospitalized (SCZ)" <- log10(data_rep$"Total time hospitalized (SCZ)"+1)
table(data_rep$"Total time hospitalized (SCZ)")[
  table(data_rep$"Total time hospitalized (SCZ)")>4]
ind_sev <- c(102,104); names(data_rep)[ind_sev]
table(data_rep$"Maternal Infection during Pregnancy, Viral")
table(data_rep$"Maternal Infection during Pregnancy, Bacterial")
```

``` r
### Standardize PRS scores

prs = names(data_rep)[grep("PGS",names(data_rep))]  
data_rep[!is.na(C1),c(prs):= lapply(.SD,scale),.SDcols=prs]
```

``` r
#### Replication of Dim 1-3 Associations 

rep <- data.table(rownames(p.vals))

p_rep <- lapply(1:nrow(rep), function(i) {  
  tmp <- as.data.frame(data_rep[is.na(case2012)])
  tmp$x  <-tmp[,rep[i,V1]]
  tmp$y <- as.matrix(tmp[,c("Dim1","Dim2","Dim3")])
if(grepl("PGS",rep[i,V1]))
  mf <- y~gender*birthdate+wave+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+x else 
    mf <- y~gender*birthdate+x 

list(manova= summary(manova(lm(mf, tmp))), N= sum(!is.na(tmp$x)) )
})

rep[,N:=unlist(lapply(p_rep,function(x) x$N))]
#rep[,Fstat:=unlist(lapply(p_rep,function(x) x[[1]]["x",4]))]
rep[,p.val:=unlist(lapply(p_rep,function(x) x$manova$stats["x",6]))]
colnames(rep)[1] <- ""
```

``` r
kable(rep,digits=c(1,1,32),format="latex",booktabs=T,
      caption="Results of MANCOVAS in replication cohort") %>%  footnote(general =  c( 
"Mancova performed with MDS dimension 1-3 as dependent variable, adjusting for age and sex.",
"For genetic varibles adjustning additionally for 10 principal components and genotype wave."))
```

``` r
#### Replication of Single Dim Associations 

rep_coef <- data.table(which(p.vals[,1:3]==apply(p.vals[,1:3],1,min),
                               arr.ind = T),keep.rownames = T )
rep_coef[,dim:=paste0("Dim",col)]
rep_coef[,coef_main:=betas[,1:3][
  which(p.vals[,1:3]==apply(p.vals[,1:3],1,min))]]
rep_coef[,2:3:=NULL,with=F]

p_rep <- lapply(1:nrow(rep_coef), function(i) {  
  tmp <- data_rep[is.na(case2012)]
  tmp$x  <-tmp[,rep_coef[i,rn],with=F]
  tmp$y <- tmp[,rep_coef[i,dim],with=F]
#if(grepl("PGS",rep[i,1]))
#  mf <- y~gender*birthdate+wave+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+x else 
    mf <- y~gender*birthdate+x 

list(anova= anova(lm(mf, tmp)),coef=coef(lm(mf, tmp)))
})


rep_coef[,N:=unlist(lapply(p_rep,function(x) sum(x[[1]][,1])+1))]
rep_coef[,Fstat:=unlist(lapply(p_rep,function(x) x[[1]]["x",4]))]
rep_coef[,p.val:=unlist(lapply(p_rep,function(x) x[[1]]["x",5]))]
rep_coef[,coef:=unlist(lapply(p_rep,function(x) x[[2]]["x"]))]
colnames(rep_coef)[1] <- ""

kable(rep_coef[,c(1:2,4,3,7), with=F],digits=c(1,1,1,3,3),format="latex",booktabs=T,
      caption="Linear regression in replication cohort") %>%  footnote(general =  c( 
"Linear regression performed for the dimension for the dimension with the stongest association.",
"This was done to test for signconcordance."))
```


Sensitivity
-----------

To test the sensitivity of the associations to different parameter settings in sequence analysis, we selected a subset of subjects that did not require imputation. We did this by only including subjects followed for more than 25 years, and by only looking at the first 25 years of the trajectories:

``` r
#seq25 <- seq[as.numeric(dt_mds[case2012==1,year]) < 1991,]
seq <- seq[rownames(seq) %in% dt[case2012==1,pid],]
seq25 <- seq[as.numeric(dt_mds[case2012==1,year]) < 1991 & 
               as.numeric(dt_mds[case2012==1,SCZ]) < 25,]
seq25 <- seq25[,seq(1,25)]
```

We then performed the sequence analysis with three different substitution cost settings:

Constant costs

``` r
sub.cost <- seqsubm(seqdata = seq25, method = "CONSTANT")
```

Jaccard Distance

``` r
alp <- seqdef(as.data.frame(alphabet(seq25)),stsep="[.]")
rownames(alp) <- alphabet(seq25)
sub.cost_jacc <- matrix(0, nrow(alp),nrow(alp))
alp <- alp[order(seqlength(alp)),]
for(i in 1:nrow(alp))
  for(j in 1:nrow(alp))
    if(i>=j){
      sub.cost_jacc[i,j] <- 1- sum(unlist(alp[i,1:seqlength(alp[i,])]) %in% 
                                     unlist(alp[j,1:seqlength(alp[j,])]))/
        length(unique(c(unlist(alp[i,1:seqlength(alp[i,])]), 
                        unlist(alp[j,1:seqlength(alp[j,])]))))} else {
                          sub.cost_jacc[i,j] <- NA}
sub.cost_jacc <- as.matrix(as.dist(sub.cost_jacc ))
rownames(sub.cost_jacc)  <- colnames(sub.cost_jacc) <- rownames(alp)
sub.cost_jacc <- sub.cost_jacc[order(rownames(sub.cost_jacc)),
                               order(rownames(sub.cost_jacc))]
rownames(sub.cost_jacc) <- rownames(sub.cost)
colnames(sub.cost_jacc) <- colnames(sub.cost)
```

1-Simple matching coefficient (Equivalent to Multichannel Sequence Analysis)

``` r
sub.cost_smc <- matrix(0, nrow(alp),nrow(alp))
for(i in 1:nrow(alp))
  for(j in 1:nrow(alp))
    if(i>=j){
      sub.cost_smc[i,j] <- sum(!unlist(alp[i,1:seqlength(alp[i,])]) %in%
                                 unlist(alp[j,1:seqlength(alp[j,])]))/8
        } else sub.cost_smc[i,j] <- NA
sub.cost_smc <- as.matrix(as.dist(sub.cost_smc ))
rownames(sub.cost_smc)  <- colnames(sub.cost_smc) <- rownames(alp)
sub.cost_smc <- sub.cost_smc[order(rownames(sub.cost_smc)),
                             order(rownames(sub.cost_smc))]
rownames(sub.cost_smc) <- rownames(sub.cost)
colnames(sub.cost_smc) <- colnames(sub.cost)
```

We then computed dissimilarities using three different sequence alignment methods: - Hamming Distance (HAM) - Optimal Matching with indelcosts of \(max(SC)/2\) (OM\_.5) - Optimal Matching with indelcosts of \(max(SC)\) (OM\_1)

Additionally we added six state distribution based dissimilarities [Deville & Saporta 1983](http://cedric.cnam.fr/~saporta/ArticleDevilleSaporta1983.pdf).

``` r
dist_list <- list(
OM_jacc_.5 =  seqdist(seq25, method = "OM", sm = sub.cost_jacc, 
                      indel = .5*max(sub.cost_jacc)),
OM_jacc_1 =  seqdist(seq25, method = "OM", sm = sub.cost_jacc, 
                     indel = max(sub.cost_jacc)),
HAM_jacc = seqdist(seq25, method = "HAM", sm = sub.cost_jacc),

OM_smc_.5 =  seqdist(seq25, method = "OM", sm = sub.cost_smc, 
                     indel = .5*max(sub.cost_smc)),
OM_smc_1 =  seqdist(seq25, method = "OM", sm = sub.cost_smc, 
                    indel = max(sub.cost_smc)),
HAM_smc = seqdist(seq25, method = "HAM", sm = sub.cost_smc),


OM_const_.5 = seqdist(seq25, method = "OM", sm = "CONSTANT", indel = 1),
OM_const_1 = seqdist(seq25, method = "OM", sm = "CONSTANT", indel = 2),
HAM_const = seqdist(seq25, method = "HAM", sm = "CONSTANT"),

eucl_2 = seqdist(seq25[,2:25], method = "EUCLID", step = 2),
chi2_2 = seqdist(seq25[,2:25], method = "CHI2", step =2 ),


eucl_5 = seqdist(seq25, method = "EUCLID", step = 5),
chi2_5 = seqdist(seq25, method = "CHI2", step =5 ),


eucl_12 = seqdist(seq25[,2:25], method = "EUCLID", step = 12),
chi2_12 = seqdist(seq25[,2:25], method = "CHI2", step =12 )
)
```

We then perfomed multidimensional Scaling and MANCOVA:

``` r
MD_list <- mclapply(dist_list, function(x) cmdscale(x, k = 3), 
                    mc.cores = 16)

# add MDS results from main analysis: 
main_dim_25 <- as.matrix(dt_mds[pid %in% rownames(seq25),.(Dim1,Dim2,Dim3)])
MD_list<- c(MD_list,list(main_dim_25))


l <- length(MD_list)
var <- rownames(betas)[-6]
m <- length(var)
pc.wave.cor <- grepl("PGS",var) 
setkey(data,pid)
data25 <-data[pid %in% rownames(seq25)] 


man_list <- lapply(1:m, function(x) lapply(1:l, function(y){
  tmp<- as.data.frame(data25)
    tmp$y1 <- MD_list[[y]]
    tmp$x1 <- tmp[,paste(var)[x]]
   if(pc.wave.cor[x])   summary(manova(y1 ~ birthdate +gender+wave+
                                         C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+x1,
                                       data=tmp )) else
      summary(manova(y1~birthdate +gender+x1,data=tmp ))}))


p.mat <- Reduce("rbind",lapply(man_list, 
                               function(x) sapply(x, function(y) 
                                 y$stats["x1",6])))

Ns <- unlist(lapply(var, function(x) sum(!is.na(data25[,paste(x),with=F]))))
row.names(p.mat) <- paste0(var, " (N=",Ns, ")")
colnames(p.mat) <- c(names(dist_list),"With imputation")
```

``` r
options(digits=2, scipen=T)
kable(p.mat,digits=100,format="latex",booktabs=T, caption="Associations of 
      significant variables across different parameter permutations in 
      25-year subset") %>% 
  kable_styling(latex_options = "scale_down") %>%  footnote(general =  c( 
"OM = Optimal Matching, HAM= Hamming Distance,",
"Eucl = Euclidian distance (state distribition (number=k)),",
"chi2= chi-squared (state distribition(number=k)),",
"jacc = 1-jaccard coefficient, smc= Simple Matching coefficent,const= Constant Costs.",
"MANCOVAs done with first three dimensions of MDS as dependent variable, adjusting for", 
"age and sex.For genetic varibles adjustning additionally for 10 principal components",
"and genotype wave."))
```

``` r
apply(p.mat,1,function(x)mean(x<.05))
```

``` r
# Onset age and number of Diagnoses
 
data25[,n_dia:=rowSums(do.call(cbind,lapply(.SD,function(x) x<25)),na.rm=T),.SDcols=
         c("F10","F30","F40","F50","F60","F70","F84","F90")]  

for(col in c("F10","F30","F40","F50","F60","F70","F84","F90")) set(data25, i=which(data25[[col]]>25), j=col, value=NA)


data25[,age_mean:=apply(cbind(F10,F30,F40,F50,F60,F70,F84,F90),1,mean,na.rm=T)]

an <- anova(lm(age_min~birthdate+gender+factor(n_dia),
            data=data25[n_dia>0]))
options(digits=3)
ggplot(data25[n_dia>0],aes(x=factor(n_dia),y=age_mean,fill=factor(n_dia)))+
  geom_boxplot(outlier.size = NA)+labs(x="Number of diagnoses",y="Mean diagnosis age",
title="Association between age of onset and number of diagnoses
in 25-year follow-up subsample") +
  theme(legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  annotate("text", label = paste0("F= ",round(an$`F value`[3],1),", p= ",
                                  format(an$`Pr(>F)`[3],digits=3)), 
           x = 3.5, y = 8, size = 5)
```

``` r


# Cl5 ~ PGS - Education, Years - comparing to random sample:
mean  <- c(
coef(lm(scale(`PGS - Education, Years`)~as.numeric(!is.na(case))+gender+birthdate+
                C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
           data=dt_ld2[pid %in% dt3[cl5==5,pid]|random==1] ))[2],
coef(lm(scale(`PGS - Education, Years`)~as.numeric(!is.na(case))+gender+birthdate+
                C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
           data=dt_ld2[pid %in% dt3[cl5==1,pid]|random==1] ))[2])

conf <- c(
confint(lm(scale(`PGS - Education, Years`)~as.numeric(!is.na(case))+gender+birthdate+
                C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
           data=dt_ld2[pid %in% dt3[cl5==5,pid]|random==1] ))[2,],
confint(lm(scale(`PGS - Education, Years`)~as.numeric(!is.na(case))+gender+birthdate+
                C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
           data=dt_ld2[pid %in% dt3[cl5==1,pid]|random==1] ))[2,])

pv  <- c(
summary(lm(scale(`PGS - Education, Years`)~as.numeric(!is.na(case))+gender+birthdate+
                C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
           data=dt_ld2[pid %in% dt3[cl5==5,pid]|random==1] ))$coefficients[2,4],
summary(lm(scale(`PGS - Education, Years`)~as.numeric(!is.na(case))+gender+birthdate+
                C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
           data=dt_ld2[pid %in% dt3[cl5==1,pid]|random==1] ))$coefficients[2,4])

label <- paste("Cl",c("5","1"))
lower <- conf[c(1,3)]
upper <- conf[c(2,4)]



df <- data.frame(label, mean, lower, upper, pv=paste("p=",round(pv,4)))

df$label <- factor(df$label, levels=rev(df$label))

fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper,label=pv)) +
        geom_pointrange() + 
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        labs(x="Schizophrenia Subgroup",y="Education PGS (95% CI)",title=
               "Estimated difference in Education PGS 
             between subgroups of SCZ patients 
             and random polulation controls") +
        theme_bw()  +
        geom_text(nudge_x = .1)

knit_hooks$set(plot = function(x, options) {
  paste("\n\\begin{figure}\n",         "\n\\begin{centering}\n",         "\\caption{", options$fig.cap, "}\n",
        "\\includegraphics[width=\\maxwidth]{",
        opts_knit$get("base.url"), paste(x, collapse = "."),
        "}\n",
        "\n\\\n
Estimates based on a linear regression of Education PGS as function of Schizophrenia diagnosis, adjusting for age, sex, 10 principal component and genotype wave. Separate regression were performed for each of the two SCZ subgroups.",
        "\n\\end{centering}\n",
        "\n\\end{figure}\n",
        sep = '')
})

fp
```

``` r
## Browser 

# Summarizing data

library(shiny)
library(ggplot2)
library(TraMineR)

#mds<- read.csv(
#  "/data/projects/IBP/gonthe/sequence_analysis/all/output/mds_trajectories_k7_n31127_181023.csv",
#  header=TRUE)
#mds_scores <- mds[which(mds$case2012==1),]
#mds_norm = as.data.frame(scale(mds_scores[,4:10]))

#mds <- cbind(pid=mds_scores$pid, case=mds_scores$case, case2012=mds_scores$case2012, mds_norm)

load(file = "~/trajectory/mds190531.Rdata")

load("/data/projects/IBP/gonthe/sequence_analysis/all/dist/seq_case-cohort_180314_thin.Rda")
load("/data/projects/IBP/gonthe/sequence_analysis/all/dist/dt_case-cohort_180314_thin.Rda")


dt_mds <- dt[case2012==1]
dt_mds[,c("Dim1","Dim2","Dim3","Dim4","Dim5","Dim6","Dim7"):=as.data.frame(wmd)]

mds <- dt_mds[,c("pid","Dim1","Dim2","Dim3","Dim4","Dim5","Dim6","Dim7")]


seq <- seq[rownames(seq) %in% mds$pid,]

# omit missing
seq <- seqdef(seq,missing="censored")

dt <- dt[dt$pid %in% mds$pid,]


mds[,"pid"] <- rownames(mds) <- rownames(seq) <- sample(1:nrow((mds))) 

# Clustering based on dimension 1-3
clust <- hclust(dist(mds[,c("Dim1","Dim2","Dim3")]),method="ward.D2")
tree_13 <-  cutree(clust,k=1:200)
min_13 <- min(which(sapply(apply(tree_13,2,table),min)<5))-1
tree2_13 <- tree_13[,1:min_13]
tree3_13 <- tree2_13[!duplicated(tree2_13[,ncol(tree2_13)]),]

# Clustering based on combinations of dimension 1-7
hclustlist <- lapply(1:7, function(x) lapply(x:7,function(y) hclust(dist(mds[,c(paste0("Dim",x),paste0("Dim",y)),with=F]),"ward.D2")))
tree <- lapply(hclustlist,function(x) lapply(x, cutree,k=1:200))
min <- lapply(tree,function(x) lapply(x,function(y) min(which(sapply(apply(y,2,table),min)<5))-1))
tree2 <- lapply(1:length(tree), function(x) lapply(1:length(tree[[x]]), function(y)  tree[[x]][[y]][,1:min[[x]][[y]]]))
tree3 <- lapply(1:length(tree2), function(x) lapply(1:length(tree2[[x]]), function(y) tree2[[x]][[y]][!duplicated(tree2[[x]][[y]][,ncol(tree2[[x]][[y]])]),]))


t <- lapply(ncol(tree2_13), function(x)
  {  tmp <-cbind(dt,cl=tree2_13[,x])
  tmp[,.(F10 = mean(!is.na(F10)), F30 = mean(!is.na(F30)), F40= mean(!is.na(F40)), 
         F50 = mean(!is.na(F50)), F60 = mean(!is.na(F60)), F70 =mean(!is.na(F84)),
         F84 = mean(!is.na(F84)), F90 = mean(!is.na(F90))),cl]})


t_all <- lapply(1:length(tree2), function(x) lapply(1:length(tree2[[x]]), function(y)
                                               lapply(ncol(tree2[[x]][[y]]), function(z)
  {  tmp <-cbind(dt,cl=tree2[[x]][[y]][,z])
  tmp[,.(F10 = mean(!is.na(F10)), F30 = mean(!is.na(F30)), F40= mean(!is.na(F40)), 
         F50 = mean(!is.na(F50)), F60 = mean(!is.na(F60)), F70 =mean(!is.na(F84)),
         F84 = mean(!is.na(F84)), F90 = mean(!is.na(F90))),cl]})
))


rare <- lapply(alphabet(seq), function(x) apply(seq,1, function(y) sum(y %in% x)))
nseq <- lapply(rare, function(x) lapply(x, function(y) sum(y>0)))
nseq2 <- lapply(nseq, function(x) sum(unlist(x)))
rare_states <- alphabet(seq)[which(nseq2<5)]
rare_seqs <- rare[which(nseq2<5)]
rare_list <- Reduce('+',rare_seqs)
sum(rare_list>0)


comb <- levels(seq$a1)
comb[levels(seq$a1) %in% rare_states]<- "Other"
for(i in 1:36) levels(seq[[i]]) <- comb

attributes(seq)$alphabet <- unique(comb)
attributes(seq)$labels <- unique(comb)

#for(i in 1:sum(rare_list>0)) 
#  seq[rare_list>0,][1,][as.numeric(seq[rare_list>0,][1,]) %in% which(nseq2<5)] <- "other"
#<-rep("censored",36) 

seq1<- seqdef(seq[,seq(3,36,3)],missing = "%")
seq<- seqdef(seq[,seq(3,36,3)],missing = "*")




rare_times <- lapply(seq, function(x) which(table(x)<4 & table(x)>0))

for(i in 1:length(rare_times))
  seq[which(as.numeric(seq[,i]) %in% unlist(rare_times[i])),i] <- "Other"    


freq_13 <- lapply(1:min_13, function(z)
    seqstatd(seq[tree2_13[,min_13]==z,], weighted=F)$Frequencies  
    )
n_13 <- lapply(1:min_13, function(z)   sum(tree2_13[,min_13]==z))

sm <- sapply(1:12, function(x) seqlength(seq)>=x)
nt_13 <- sapply(1:12, function(z) sapply(1:min_13, function(y) sum(tree2_13[sm[,z],min_13]==y)))


freq_bl <- lapply(1:length(tree), function(x) lapply(1:length(tree[[x]]), function(y)  lapply(1:min[[x]][[y]], function(z) seqstatd(seq[tree2[[x]][[y]][,min[[x]][[y]]]==z,],weighted=F)$Frequencies  )))
n_bl <- lapply(1:length(tree), function(x) lapply(1:length(tree[[x]]),function(y)  lapply(1:min[[x]][[y]], function(z)   sum(tree2[[x]][[y]][,min[[x]][[y]]]==z))))

nt_bl <- lapply(1:length(tree), function(x) lapply(1:length(tree[[x]]),function(y)   
  sapply(1:12, function(z) sapply(1:min[[x]][[y]], function(z1) sum(tree2[[x]][[y]][sm[,z],min[[x]][[y]]]==z1)))))
  


cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
          "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
palette(cols)

alpa <- alphabet(seq)

setwd("shiny")


freq_13 <- lapply(freq_13, function(x) {
         x2<- x
         x2[is.na(x2)] <- 0
         colnames(x2) <- sub("a","",colnames(x2))
         x2}) 
freq_bl <- lapply(1:length(freq_bl), function(x) lapply(1:length(freq_bl[[x]]), function(y) {
         
         y2<- freq_bl[[x]][[y]]
         lapply(y2, function(z){
         y3<- z   
         y3[is.na(y3)] <- 0
         colnames(y3) <- sub("a","",colnames(y3))
         y3})}
         ))


save(tree3,tree2, freq_bl, n_bl, nt_bl,tree3_13,tree2_13, t, t_all, freq_13, n_13, nt_13, cols,mds, alpa, file="data.Rda") 
rm(list=ls())
```

``` r
# App Script:

library(shiny)
library(ggplot2)
library(TraMineR)
library(gridExtra)
library(scales)
load("data.Rda")
source("seqmodst2.R")
library(RColorBrewer)


ui <- fluidPage(
  headerPanel('Schizophrenia Trajectories'),
  sidebarPanel(
    radioButtons('clustering_data', 'Clustering', c("Dimentions 1-3", "Selected Dimensions")),
    selectInput('xcol', 'X Variable', names(mds)[2:8]),
    selectInput('ycol', 'Y Variable',  names(mds)[2:8],
                selected =  names(mds)[2:8][[2]]),
    numericInput('clusters', 'Groups', 5,
                 min = 1, max = 300,  ),
    radioButtons('type', 'Plot type', c("Barplot","Modal state","Modal state - all", "Chronogram"),selected = "Modal state - all"),
  conditionalPanel(condition = "input.type == 'Barplot'",
                     checkboxGroupInput("com", "Comorbidity",c("F10","F30","F40","F50","F60","F70","F84","F90"),selected=c("F10")))),

  mainPanel(
    p("Multidimensional scaling and clustering of the Disease Trajectories:"),
    verbatimTextOutput("info"),
    plotOutput('plot1', click = "plot_click"),
    p("Cluster Characteristics:"),
    plotOutput('plot2',height = 200),
    plotOutput('plot3'),
    p("All states observed in <5 individuals are labeled as", span("other.",style="color:brown"), 
    "All states observed in <4 individuals at a given age are also labeled as",
    span("other.",style="color:brown"))
  )
)

server <- function(input, output) {
  
  selectedData <- reactive({
    mds[, c(input$xcol, input$ycol),with=F]
  })
  
  clusters <- reactive({
  if(input$clustering_data=="Dimentions 1-3")
    {if(input$clusters>ncol(tree3_13)) NULL else tree2_13[,input$clusters]}
  else {
  x <- which( names(mds)[2:8] %in% input$xcol)
  y <- which( names(mds)[2:8] %in% input$ycol)
  ind1 <- min(x,y)
  ind2 <- abs(x-y)+1
  if(input$clusters>ncol(tree3[[ind1]][[ind2]])) NULL  else tree2[[ind1]][[ind2]][,input$clusters]
    }})

  df <-reactive({
    dft <- cbind(selectedData(),factor(clusters()))
  colnames(dft) <- c("x","y","cl")
  clo <- order(tapply(dft$x,dft$cl,mean))
  dft$cl <- factor(clusters(), levels = clo, label = 1:length(clo))
  list(dtf=dft,clo=clo)})

  
  output$info <- renderText({
    x <- which( names(mds)[2:8] %in% input$xcol)
    y <- which( names(mds)[2:8] %in% input$ycol)
    ind1 <- min(x,y)
    ind2 <- abs(x-y)+1
    if(is.null(clusters())) paste0("maximum number of clusters for these dimensions is: ", ncol(tree3[[ind1]][[ind2]])) else NULL
    })

  click_saved <- reactiveValues(singleclick =NULL)
  observeEvent(eventExpr = input$plot_click, handlerExpr = { click_saved$singleclick <- input$plot_click })
  np <-  reactive({
   # if(!is.null(click_saved$singleclick))
    nearPoints(df()$dtf,click_saved$singleclick, threshold = 100,
               maxpoints = 1,
               addDist = TRUE) #else nearPoints(df(),list(x=0,y=0,col=1), threshold = 100,maxpoints = 1,addDist = TRUE)
  })  
  
  output$plot1 <- renderPlot({
    par(mar = c(5.1, 4.1, 0, 1))
    g <- ggplot(df()$dtf,aes(x=x,y= y,col = factor(cl)))+
    geom_point()+
    theme(legend.position="none",
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black")
                      )
    #if(nrow(np())<1) np <- matrix(rep(0,4),nrow = 1) else
    np <- np() 
    #print(np)
    
    if(!input$type %in% c("Barplot","Modal state - all")){
      cla <- df()$dtf$cl
      np1 <- cla[which(rownames(mds) %in% rownames(np))] 
      print(np1)
      if(nrow(np())<1) cl <- cla %in% 1 else cl <- cla %in% np1      
      df <- df()$dtf
      g <- g+stat_ellipse(data=df[cl,])
    } 
    if(input$clusters<10) {pal <- brewer.pal(input$clusters,"Set1")
                           names(pal) <- 1:input$clusters
      g <- g + scale_color_manual(values= pal)}
    g
    })



  cl <- reactive({
    if(nrow(np())<1) df()$dtf$cl[which(df()$dtf$cl %in% 1)][1] else df()$dtf$cl[which(rownames(mds) %in% rownames(np()))][1]
  })
  
  
  
  seq1 <- reactive({
     if(input$type %in% c("Modal state", "Chronogram")){ 
    x <- which( names(mds)[2:8] %in% input$xcol)
    y <- which( names(mds)[2:8] %in% input$ycol)
    ind1 <- min(x,y)
    ind2 <- abs(x-y)+1
     if(input$clustering_data=="Dimentions 1-3"){
       c <-which(factor(tree3_13[,input$clusters], levels = df()$clo, label = 1:length(df()$clo))==cl())
       
       #freq_x <- Reduce('+',lapply(c, function(x) freq_131[[x]]* n_13[[x]]))/Reduce('+',n_13[c])
       freq_x <- t(t(Reduce('+',lapply(c, function(x) t(t(freq_13[[x]])* nt_13[x,]))))/
                     colSums(nt_13[c,,drop=F]))[,1:11]
       n_x <- Reduce('+',n_13[c])    } else {
          c <-which(factor(tree3[[ind1]][[ind2]][,input$clusters], levels = df()$clo, 
                           label = 1:length(df()$clo))==cl())
          freq_x <- t(t(Reduce('+',lapply(c, function(x) t(t(freq_bl[[ind1]][[ind2]][[x]])* 
                                                             nt_bl[[ind1]][[ind2]][x,]))))/
                        colSums(nt_bl[[ind1]][[ind2]][c,,drop=F]))[,1:11]
          n_x <- Reduce('+',n_bl[[ind1]][[ind2]][c])   }
    freq_x[is.na(freq_x)]<-0
    cpal <- rep("grey", length(alpa))
    stat<-seqstatl(seqmodst2(freq_x,n=n_x,weighted=F))
    cpal[alpa %in% stat] <-  c(cols[1:length(stat)])
    cpal[alpa=="None"]<-"lightgray" 
    cpal[alpa=="Other"]<-"black" 
    #cpal[alpa=="censored"]<-"white"
    list(freq=freq_x,n=n_x,stat=stat,cpal=cpal)
    } else 0
    
  })

  
  seq_all <- reactive({
    if(input$type=="Modal state - all"){    
    lapply(1:input$clusters,function(cl){
    x <- which( names(mds)[2:8] %in% input$xcol)
    y <- which( names(mds)[2:8] %in% input$ycol)
    ind1 <- min(x,y)
    ind2 <- abs(x-y)+1
     if(input$clustering_data=="Dimentions 1-3"){
        c <-which(factor(tree3_13[,input$clusters], levels = df()$clo, 
                         label = 1:length(df()$clo))==cl)
       
       
       #freq_x <- Reduce('+',lapply(c, function(x) freq_131[[x]]* n_13[[x]]))/Reduce('+',n_13[c])
       freq_x <- t(t(Reduce('+',lapply(c, function(x) t(t(freq_13[[x]])* nt_13[x,]))))/colSums(nt_13[c,,drop=F]))[,1:11]
       n_x <- Reduce('+',n_13[c])    } else {
           c <-which(factor(tree3[[ind1]][[ind2]][,input$clusters], levels = df()$clo, 
                           label = 1:length(df()$clo))==cl)
          freq_x <- t(t(Reduce('+',lapply(c, function(x) t(t(freq_bl[[ind1]][[ind2]][[x]])* nt_bl[[ind1]][[ind2]][x,]))))/
                        colSums(nt_bl[[ind1]][[ind2]][c,,drop=F]))[,1:11]
          n_x <- Reduce('+',n_bl[[ind1]][[ind2]][c])   }
    freq_x[is.na(freq_x)]<-0
    cpal <- rep("grey", length(alpa))
    stat<-seqstatl(seqmodst2(freq_x,n=n_x,weighted=F))
    cpal[alpa %in% stat] <-  c(cols[(1:length(stat))+c(1,3,4,6,7,8)[cl]])
    cpal[alpa=="None"]<-"lightgray" 
    cpal[alpa=="Other"]<-"black" 
    #cpal[alpa=="censored"]<-"white"
    list(freq=freq_x,n=n_x,stat=stat,cpal=cpal)
    })} else 0
  })

  cpal_all <- reactive({
    if(input$type %in% c("Modal state - all")){ 
      col <-  seq_all()[[1]]$cpal[alpa %in%  seq_all()[[1]]$stat] 
      st <- seq_all()[[1]]$stat
      if(length(seq_all())>1) for(i in 2:length(seq_all())){ 
      #cp[alpa %in%  seq_all()[[i]]$stat] 
      #<- seq_all()[[i]]$cpal[alpa %in%  seq_all()[[i]]$stat]
      st <- c(st, seq_all()[[i]]$stat)
      col <- c(col, seq_all()[[i]]$cpal[alpa %in%  seq_all()[[i]]$stat] )
      }
      col <- col[!duplicated(st)]
      st <- st[!duplicated(st)]
      col <- col[order(st)]
      st <- st[order(st)]
      
      cp <-  seq_all()[[1]]$cpal
      cp[alpa %in% st] <-col
      list(col=col
           ,stat=st,cpal= cp 
           ) 
           } else 0})
  
  
  output$plot2 <- renderPlot({
    par(mar = c(5.1, 4.1, 0, 1))
    if(input$type=="Modal state")   {
      plot(seqmodst2(seq1()$freq,n=seq1()$n,weighted=F,freq.mean=T,with.missing=F),cpal=seq1()$cpal,xlab="age") 
                } else if(input$type=="Modal state - all"){      
                   if(input$clusters > 6)  { plot(0,type='n',axes=FALSE,ann=FALSE) 
          text(.9,.5," <-- Maximum number of clusters in one plot is six, select single cluster plots") } else {
            if(input$clusters %in% 1:3) par(mfrow=c(input$clusters,1),oma=c(0,0,1.5,0))  else #if(input$clusters %in% 4:6) 
                     par(mfrow=c(ceiling(input$clusters/2),2),oma=c(0,0,1.5,0))
                  for(i in 1:input$clusters) { 
                  plot.stslist.modst2(seqmodst2(seq_all()[[i]]$freq,n=seq_all()[[i]]$n,weighted=F,freq.mean=T,with.missing=F),#cpal=seq_all()[[i]]$cpal
                                      cpal=cpal_all()$cpal,ylab="",sub="",cex.axis=2,main="",tx=F,yaxis=F,xlab="age")  
                  par(xpd = T)
                   if(input$clusters<10)  points(x = -1,y = .8,col = brewer.pal(input$clusters,"Set1")[i], pch = 16,cex=3) else
                  points(x = -1,y = .8,col = hue_pal()(input$clusters)[i], pch = 16,cex=3)
                  text(x = -1,y = 0,cex=1.5,labels = paste0("n=",seq_all()[[i]]$n))
            } 
            mtext(paste("Modal state for each cluster"),outer = T,cex=1.2) }
          
          } else if(input$type=="Chronogram"){
      res<- list(seq1()$freq,n=seq1()$n)
    names(res) <- c("Frequencies", "ValidStates")
      class(res) <- c("stslist.statd", "list")
      attr(res, "nbseq") <- seq1()$n
      attr(res, "cpal") <- seq1()$cpal
      attr(res, "xtlab") <- colnames(seq1()$freq)
      attr(res, "weighted") <- F
      plot(res,weighted=F)} else { 
          if(length(input$com)==0) { plot(0,type='n',axes=FALSE,ann=FALSE) 
          text(.9,.5," <-- Select Comorbidities to plot") } else {
          if(input$clustering_data=="Dimentions 1-3")  { 
            tmp <- t[[1]] 
            tmp$cluster <- factor(tree3_13[,input$clusters] ,levels = df()$clo, 
                           label = 1:length(df()$clo))
            tmp$n <- unlist(n_13)
            }else { 
            x <- which( names(mds)[2:8] %in% input$xcol)
            y <- which( names(mds)[2:8] %in% input$ycol)
            ind1 <- min(x,y)
            ind2 <- abs(x-y)+1
            tmp <- t_all[[ind1]][[ind2]][[1]]
            tmp$cluster <- factor(tree3[[ind1]][[ind2]][,c(input$clusters)],levels = df()$clo, 
                           label = 1:length(df()$clo))
            tmp$n <- unlist(n_bl[[ind1]][[ind2]])}
            p <- lapply(input$com, function(z) {
              tmp$c <- tmp[,paste(z),with=F]
              tmp2 <- tmp[,.(Prop =sum(c*n)/sum(n)),cluster]
              p1 <- ggplot(tmp2, aes(x=cluster, y=Prop,fill=cluster))+
                geom_bar(stat="identity")+labs(title= paste(z),x="Cluster")+
                theme(legend.position="none",
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black")
                      )
              if(input$clusters<10)  p1 <- p1 + scale_fill_manual(values = 
                                                                  brewer.pal(input$clusters,"Set1"))
              p1                                                  })

            do.call('grid.arrange',c(p))}}
        })

#   library(RColorBrewer)
# myColors <- brewer.pal(5,"Set1")
# names(myColors) <- levels(dt3[,cl5])
# colScale <- scale_colour_manual(name = "Cluster",values = myColors)
# fillScale <- scale_fill_manual(name = "Cluster",values = myColors)
    

  output$plot3 <- renderPlot({
    par(mar = c(5.1, 4.1, 0, 1))
  if(input$type=="Barplot") return()  else
  if(input$type=="Modal state"){   alp <- seq1()$stat
    col1 <- seq1()$cpal[alpa %in% seq1()$stat] }else
  if(input$type=="Modal state - all"){   
    alp <- cpal_all()$stat
    #col1 <- cpal_all()$cpal[alpa %in% cpal_all()$stat] 
    col1 <- cpal_all()$col
    }else if(input$type=="Chronogram") {
      alp <- alpa[apply(seq1()$freq,1,max)>.1]
      col1 <- seq1()$cpal[apply(seq1()$freq,1,max)>.1]}
    plot(0,type='n',axes=FALSE,ann=FALSE)
     order <- order(alp)
      alp <- alp[order] 
      col1 <- col1[order]
      nn <- c("F10 - Subst. Abuse",
          "F30 - Affective Disord.",
          "F40 - Neurotic Disord.",
          "F50 - Eating Disord.",
          "F60 - Personality Disord.",
    "F70 - Intellect. Disability",
    "F84 - Autism Spectrum Disord.",
    "F90 - Childhood Disord.")
      alp <- replace(alp,which(alp %in% c("F10","F30","F40","F50","F60","F70","F84","F90")),
                     nn[c("F10","F30","F40","F50","F60","F70","F84","F90") %in% alp])
    legend( .8,1, alp ,fill = col1,cex=1)
    
  })


  
}

shinyApp(ui = ui, server = server)
```
