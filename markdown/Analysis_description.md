This is a document providing the code for the analyses in the paper "..." by Krebs et al.

Sequence Analysis in Schizophrenia
==================================

Libraries and Data
------------------

Load required libraries including

``` r
library(devtools)
install_github("MortenKrebs/diagtraject")
library(diagtraject)
```

Load data and create sequence state object:

The Study was performed using the iPSYCH cohort (See [Pedersen et al](https://www.nature.com/articles/mp2017196):

``` r
# 3 data.frames with one row pr individual and collums containing: 
# id, birthday and ##date of first diagnosis for each category of diagnoses

load("diagdates.Rda") # individuals diagnosed with Schizophrenia before Dec 31, 2012  
load("diagdates_random.Rda") # random population sample
load("diagdates2016.Rda") # individuals diagnosed with Schizophrenia 2013-2016
```

``` r
colnames(diagdates)
```

    ##  [1] "pid"       "birthdate" "F1"        "F3"        "F4"       
    ##  [6] "F50"       "F60"       "F70"       "F84"       "F9"       
    ## [11] "SCZ"       "gender"

``` r
nrow(diagdates)
```

    ## [1] 5432

\begin{table}[H]
\centering
\begin{tabular}{>{\raggedright\arraybackslash}p{4.5cm}>{\raggedright\arraybackslash}p{4.5cm}>{\raggedright\arraybackslash}p{4.5cm}}
\toprule
Diagnosis & ICD.10 & ICD.8\\
\midrule
Substance abuse & F10-F19 & 291.x9, 294.39, 303.x9, 303.20, 303.28, 303.90, 304.x9\\
Mood disorders* & F30-39 & 296.x9 (excl 296.89), 298.09, 298.19, 300.49, 301.19\\
Anxiety disorders and Obsessive compulsive disorder & F40.0-F40.2 F41.0-F41.1 F42 F43.0-F43.1 & 300.39\\
Eating disorders & F50 & 305.60, 306.50, 306.58, 306.59\\
Personality disorder & F60 & 305.60, 306.50, 306.58, 306.59\\
\addlinespace
Mental Retardation & F70-79 & 301.x9 (excl 301.19), 301.80, 301.81, 301.82, 301.84\\
Pervasive developmental disorders & F84 & 299.00, 299.01, 299.02, 299.03\\
Behavioural and emotional disorders with onset usually occurring in childhood and adolescence & F90-98 & 306.x9 308.0x\\
\bottomrule
\end{tabular}
\end{table}
\*For recurrent depression onset was defined as the second admission that occurred at least 8 weeks after last discharge with these ICD-8 codes

Diagnosis State Sequences:
--------------------------

``` r
dt <- data.table(diagdates)

summary(xtabs(N~.,dt[,.N,.(F1=!is.na(F1),F3=!is.na(F3),
                           F4=!is.na(F4),F50=!is.na(F50),
                           F60=!is.na(F60),F70=!is.na(F70),
                           F84=!is.na(F84),F9=!is.na(F9))]))
## Call: xtabs(formula = N ~ ., data = dt[, .N, .(F1 = !is.na(F1), F3 = !is.na(F3), 
##     F4 = !is.na(F4), F50 = !is.na(F50), F60 = !is.na(F60), F70 = !is.na(F70), 
##     F84 = !is.na(F84), F9 = !is.na(F9))])
## Number of cases in table: 5432 
## Number of factors: 8 
## Test for independence of all factors:
##  Chisq = 2046, df = 247, p-value = 2e-280
##  Chi-squared approximation may be incorrect
```

``` r
fdate <- function(x) as.Date(x,"%d/%m/%Y" )

dt[,"birthdate":= fdate(birthdate) ]
dt[,c("F1","F3", "F4","F50","F60", "F70", "F84", "F9","SCZ") :=
     lapply(list(F1,F3,F4,F50,F60,F70,F84,F9,SCZ),function(x){  
       (as.numeric(fdate(x))- as.numeric(birthdate))/365})]
#Calculate age in 31/12/2016
dt[,censored:= (as.numeric(fdate("31/12/2016"))-as.numeric(birthdate))/365] 

dt.m <-melt(dt, id.vars = 'pid', direction = "long", 
            measure.vars  = list(c(3:10,13)), value.name = "time",
            variable.name = "event")
dt.m<- dt.m[order(pid)]
dt.m[,time:=time1]
dt.m <- dt.m[!is.na(time)]

events <- levels(dt.m$event)
drop <-  matrix(FALSE, nrow = length(events), ncol = length(events),
                dimnames = list(events,events))
drop['censored',] <- T
drop[,'censored'] <- T
diag(drop) <- F

e2sm <-seqe2stm(events = events, dropMatrix = drop)
```

Creating sequence object with three different settings:

### Varying sequences lengths

``` r
seq_mis <- seqdef(TSE_to_STS(dt.m,
                             id = "pid", timestamp = "time", event= "event",
                             tmin=1, tmax=ceiling(max(dt.m$time)),stm = e2sm)
                  , firstState ="None", missing="censored")
```

``` r
range(seqlength(seq_mis))
## [1] 16 35
length(alphabet(seq_mis))
## [1] 193
```

### Plotting most frequent states:

Top-20 most frequent sequences in first 25 years for individuals followed more than 25 years

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
                     c("Substance Abuse", "Affect. Disord.",
                       "Neurot. Disord.", "Personal. Disord."), 
                     sep= " - ")


attributes(seq25)$cpal <- rep("grey", length(attributes(seq25)$labels) )
attributes(seq25)$cpal[which(
  lapply(attributes(seq25)$labels, function(x) sum(
    grepl(x, rownames(attributes(tab)$weights))))>0)]<- c(
      cols[c(1,2,5,3,7,4,6)], "white")

par(mfrow=c(1,2))
seqfplot(seq25,tlim=20:1, pbarw=F, with.legend=F, yaxis="pct")
plot(0,type='n',axes=FALSE,ann=FALSE)
legend(x=0.7,y=1, lab[o],
      fill = c(cols[c(1,2,5,3,7,4,6)][o], "white"),cex=.8)
```

![](Analysis_description_files/figure-markdown_github/seqfplot-1.pdf)

### "censored" as separate state

For later use a version keeping `"censored"` as a seperate state was created:

``` r
seq <- seqdef(TSE_to_STS(dt.m,id = "pid", timestamp = "time", 
                         event="event",
                         tmin=1, tmax=ceiling(max(dt.m$time)),
                         stm = e2sm),
              firstState ="None")
```

### Additional last state

Identifying the state at time of censoring:

``` r
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

To obtain population-wide estimates of transition probabilities, the random population cohort of 30000 individuals was included (See [Pedersen et al](https://www.nature.com/articles/mp2017196) for details) and estimates were computed using inverse sampling probability weighting (Bowman II):

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

Distribution of age and gender in the Danish population publicly available from [Statistics Denmark](www.statistikbanken.dk):

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

Inverse sampling probability weighting:

``` r
sample_dist <- sample_dist[order(gender,year)]
sample_dist[,pop_weights := DK_pop[,n]/sample_dist[,N]]
dt_pop <- merge(dt_pop, sample_dist, by=c("year","gender"))
dt_pop[!is.na(SCZ), pop_weights := 1]
setkey(dt_pop,pid)
```

Formatted as [above](#Last)

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
dt.m[,time:=time1]
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

Time-varying Transition Rates computed with `TraMineR::seqtrate:`:

``` r
seq_mis_last_uni_pop <- unique(seq_mis_last_pop)
attributes(seq_mis_last_uni_pop)$weights <- match(
  seqconc(seq_mis_last_pop),seqconc(seq_mis_last_uni_pop))

tr_t <- seqtrate(seq_mis_last_uni_pop, time.varying = T,weighted = T)
```

### Substitution Costs:

Substitution Cost Matrix, using Jaccard distance:

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
```

    ## [1] 202 202

``` r
kable(sub.cost_jacc[1:5,1:5],format="latex",booktabs=T)
```

\begin{tabular}{lrrrrr}
\toprule
  & censored-> & F1-> & F1.F3-> & F1.F3.F4-> & F1.F3.F4.F50->\\
\midrule
censored-> & 0 & 0.00 & 0.00 & 0.00 & 0.00\\
F1-> & 0 & 0.00 & 0.50 & 0.67 & 0.75\\
F1.F3-> & 0 & 0.50 & 0.00 & 0.33 & 0.50\\
F1.F3.F4-> & 0 & 0.67 & 0.33 & 0.00 & 0.25\\
F1.F3.F4.F50-> & 0 & 0.75 & 0.50 & 0.25 & 0.00\\
\bottomrule
\end{tabular}
### Handling right-censoring:

To calculate dissimilarities between sequences with right censoring, we used inferred states weighted by the probabilities of that state, given the last observed state in the sequence. When calculating the dissimilarity between two sequences \(i\) and \(j\) of unequal length, the dissimilarity, \(D(i,j)\)
\[ D(i,j)  = d_{obs} + d_{inf} \] where \(d_{obs}(i,j)\) is the dissimilarity between the sequences before right censoring occurs and \(d_{inf}\) is the sum of the substitution cost matrix weighted by the inferred probabilities.

Dissimilarities for observed states using optimal matching (OM):

``` r
  sms <- which(rownames(sub.cost_jacc) %in% 
                 rownames(seqsubm(seq,"CONSTANT")))
  d_OM <- seqdist(seq, method = "OM", indel=.5, 
                  sm = sub.cost_jacc[sms,sms])
```

### Imputation:

\[ d_{inf}(i,j) = \sum\limits_{t=1}^{t_{max}} \sum Pr(i)_t Pr(j)_t^T \circ  SC \] And for imputed states using `diagtraject::mis.cost()` :

``` r
ms <- mis.cost(seq_mis, last, tr_t, sm = sub.cost_jacc, 
               cens.type = "right", imp.length = "max", 
               diag=F, sum_to_1=T, resol.comp = resol.comp, 
               resol.ratio = resol.reduc, mc.cores=27)
```

### Overall dissimilarities

Obtained by \(d_{obs} + d_{inf}\)

``` r
dist_OM <- d_OM + as.matrix(ms$dist)
```

Multidimentional scaling
------------------------

Metric multidimentional scaling using `vegan::wcmdscale()` :

``` r
# finding unique sequences:
mcor <- match(seqconc(seq_mis_last),seqconc(unique(seq_mis_last))) 
uni <- (!duplicated(mcor))
dist <- list("dist"=dist_OM[uni,uni],"weight"= table(mcor),"mcor"= mcor)
```

``` r
library(vegan)
library(ggplot2)
wmd <- wcmdscale(dist$dist,w = table(dist$mcor), eig = T, k=13)
```

### Proportion of variance explained for `k=1:13`:

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

![](Analysis_description_files/figure-markdown_github/R2-1.pdf)

### Stressplots for `k=2:13`:

``` r
wmd_all <-wcmdscale(dist_unique_OM$dist,w = table(mcor), eig = T)
par(mfrow=c(3,4))
for(i in 2:13)
stressplot(wmd_all, k=i, p.col="blue", l.col="red", lwd=2)
```

![](Analysis_description_files/figure-markdown_github/Stress-1.pdf)

### Bootstrap Stability of MDS:

Bootstrap stability of MDS over \(n=100\) bootstrap samples of the dissimilarity matrix. Where stability is measured as the stability coefficient:

\[ ST = 1- \frac{\sum_{i=1}^{n}  \parallel X_{i}^{*}-\bar{X}^{*} \parallel^2}{\sum_{i=1}^{n}  \parallel X_{i}^{*} \parallel ^2} \]

where \(X_{i}^{*}\) is the bootstrap MDS solutions after a Procrustes transformation and \(\bar{X}^{*}\) is the average result across all bootstap solutions.

It can be interperpreted as the ratio of between and total variance (See [de Leeuw](https://link.springer.com/content/pdf/10.1007%2FBF01896814.pdf))

Computed using `diagtraject::bootmds()`:

``` r
bootmds(dist_unique_OM$dist,k = 7, nrep = 100, w = table(mcor))
```

    ## [1] 1

Correlation with gender and year of birth
-----------------------------------------

``` r
dt_mds <- cbind(dt, wmd[[6]]$points[mcor,])
dt_mds[,n_dia := dt.m[,.N,pid]$N-1]
```

Dimension 1+2 and number of diagnoses:

``` r
ggplot(dt_mds, aes(x=Dim1, y=Dim2, color=factor(n_dia)))+  
  geom_jitter(h=.05,w=.05)+ 
  scale_colour_brewer(palette = "Spectral")+
  labs(x="Dimension 1", y="Dimension 2",
       color= "Number of comorbid diagnoses")
```

![](Analysis_description_files/figure-markdown_github/Comorb_count-1.pdf)

Dimension 1 and year of birth

``` r
first_dia <- melt(dt_mds, id.vars = 'pid', 
                  direction = "long", measure.vars =  list(c(3:12)),
                  value.name = "time")
setkey(first_dia,pid)
setkey(dt_mds,pid)
dt_mds[,age_min := first_dia[,min(time1,na.rm=T),pid]$V1]
```

``` r
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
  annotate("text",x=as.Date("1996",format = "%Y"), 
           y=0, parse=T,color= "blue",
           label = as.character(expression(
             Dim%~%birthdate+onset_age+gender)))+
  annotate("text",x=as.Date("1998",format = "%Y"), 
           y=5, parse=T,color= "red",
           label = as.character(expression(Dim%~%birthdate)))
p
```

![](Analysis_description_files/figure-markdown_github/DimYear-1.pdf)

Dimension 1 and 2 ~ gender

``` r
ggplot(dt_mds, aes(x=Dim1, y=Dim2, 
                   color=factor(gender,labels = c("female","male"))))+
  geom_jitter(h=.05,w=.05)+ 
  scale_colour_manual(values = c("red","blue"))+
  labs(x="Dimension 1", y="Dimension 2",color= "Sex")
```

![](Analysis_description_files/figure-markdown_github/GenderDim-1.pdf)

Dimension Characteristics:
--------------------------

Characteristics for Dimension 1-7:

``` r
dt_mds[,paste0("Dim",1:7,"r") := lapply(list(
  Dim1,Dim2,Dim3,Dim4,Dim5,Dim6,Dim7), rank,ties.method="random")]
dt_mds[,paste0("Dim",1:7,"f") := lapply(list(
  Dim1r,Dim2r,Dim3r,Dim4r,Dim5r,Dim6r,Dim7r), 
  function(x) cut(x, breaks=50))]

dim <- paste0("Dim",1:7,"f")

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
    scale_fill_brewer(palette = "YlGn")+
    scale_y_continuous(limits =c(0,100))+
    labs(title="SA", x=sub("_f","",x),y="%")+ 
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
    scale_fill_brewer(palette = "YlGn")+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Aff", x= sub("_f","",x),y="%")+ 
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
    scale_fill_brewer(palette = "YlGn")+
    scale_y_continuous(limits =c(0,100))+
    labs(title="Neuro", x= sub("_f","",x),y="%")+ 
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
    scale_fill_brewer(palette = "YlGn")+
    scale_y_continuous(limits =c(0,100))+
    labs(title="ED", x= sub("_f","",x),y="%")+ 
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
    scale_fill_brewer(palette = "YlGn")+
    scale_y_continuous(limits =c(0,100))+
    labs(title="PD", x= sub("_f","",x),y="%")+ 
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
    scale_fill_brewer(palette = "YlGn")+
    scale_y_continuous(limits =c(0,100))+
    labs(title="MR", x= sub("_f","",x),y="%")+
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
    scale_fill_brewer(palette = "YlGn")+
    scale_y_continuous(limits =c(0,100))+
    labs(title="ASD", x= sub("_f","",x),y="%")+
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
    scale_fill_brewer(palette = "YlGn")+
    scale_y_continuous(limits =c(0,100))+
    labs(title="ChD", x= sub("_f","",x),y="%", 
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
    labs(title="", x= sub("_f","",x),y="%", fill="Number of diagnoses")+ 
    theme(legend.direction="vertical",
          #legend.title=element_blank(),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
  guides(fill=guide_legend(ncol=2)) })

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

library(grid)
library(gridExtra)

emp <-textGrob("")

# All: 
plist_all<- c(p_n,list(emp),p_10,list(leg),p_30,list(emp),
              p_40,list(leg2),p_50,list(emp),p_60,list(emp),
              p_70,list(emp),p_84,list(emp),p_90,list(emp))
#Dim1
plist1<- list(p_n[[1]], leg,leg2, p_30[[1]],p_60[[1]],p_90[[1]])
#Dim2
plist2<- list(p_n[[2]], leg,leg2,p_30[[2]],p_84[[2]],p_90[[2]])
#Dim3
plist3<- list(p_n[[3]], leg,leg2,p_10[[3]],p_30[[3]],p_84[[3]])

do.call("grid.arrange", c(plist1, ncol=3))
```

![](Analysis_description_files/figure-markdown_github/Dim1_7-1.pdf)

``` r
do.call("grid.arrange", c(plist2, ncol=3))
```

![](Analysis_description_files/figure-markdown_github/Dim1_7-2.pdf)

``` r
do.call("grid.arrange", c(plist3, ncol=3))
```

![](Analysis_description_files/figure-markdown_github/Dim1_7-3.pdf)

``` r
do.call("grid.arrange", c(plist_all, ncol=8))
```

![](Analysis_description_files/figure-markdown_github/Dim1_7-4.pdf)

Clustering
----------

Hierarchical clustering based on Sequence Analysis:

``` r
library(WeightedCluster)

c <- dist_unique_OM$mcor[which(dt[,!is.na(case2012)])]
clust <- hclust(as.dist(dist_unique_OM$dist[c,c]),method="ward.D2")
dt3 <- dt[!is.na(case2012)]
dt3[,cl2:= factor(cutree(clust,2))]
dt3[,cl3:= factor(cutree(clust,3))]
dt3[,cl4:= factor(cutree(clust,4))]
dt3[,cl5:= factor(cutree(clust,5))]
dt3[,cl6:= factor(cutree(clust,6))]

dt3 <- merge(dt_mds, dt3[,.(pid,cl2,cl3,cl4,cl5,cl6)],by="pid")

ggplot(dt3, aes(x=Dim1, y=Dim2, color=cl6))+  geom_jitter(h=.05,w=.05)
```

![](Analysis_description_files/figure-markdown_github/hclust-1.pdf) Cluster Statistics using `WeightedCluster::wcClusterQuality`

``` r
wcClusterQuality(as.dist(dist_unique_OM$dist[c,c]),dt3[,cl6])
```

    ## $stats
    ##    PBC     HG   HGSD    ASW   ASWw     CH     R2   CHsq   R2sq     HC 
    ##   0.43   0.59   0.59   0.17   0.17 415.24   0.28 795.05   0.42   0.18 
    ## 
    ## $ASW
    ##      ASW   ASWw
    ## 1  0.450  0.451
    ## 2  0.066  0.067
    ## 3  0.028  0.029
    ## 4  0.263  0.264
    ## 5 -0.050 -0.049
    ## 6  0.352  0.352

Assesment of clusterstability in bootstrap test using 100 permutations of data and computing the mean jaccard index using a `fpc::clusterboot`:

``` r
library(fpc)
boot<-clusterbootw(as.dist(dist_unique_OM$dist[c,c]),
                   clustermethod=disthclustCBI, k= 6, 
                   B=100, method="ward.D2", members=NULL, mc.cores = 27)
```

``` r
boot$bootmean
```

    ## [1] 0.53 0.59 0.61 0.65 0.53 0.83

Visualizing stability using `WGCNA`:

``` r
library(WGCNA)
labels <- matrix(NA, ncol=length(boot$partition), nrow=101)
labels[1,] <-  boot$partition
for(i in 2:101) labels[i,] <- matchLabels(boot$bootpartition[,(i-1)],
                                          boot$partition)

colors <- labels2colors(labels)
plotDendroAndColors(boot$result$result,t(colors), main="",
                    c("Full dataset", paste("Resamp.", 
                                            c(1:(dim(labels)[2]-1)))),
                    dendroLabels =F, hang=.03,autoColorHeight = T)
```

![](Analysis_description_files/figure-markdown_github/clusterboot-1.pdf)

Association with putative risk variables
----------------------------------------

Load additional genetic and registry variables

``` r
load("additional.Rda")
```

### Polygenic Scores:

\begin{table}

\caption{\label{tab:unnamed-chunk-39}Summary Statistics used for
      Polygenic Score Calculations}
\centering
\begin{tabular}[t]{lrrrl}
\toprule
Phenotype & n.cases & n.controls & PMID & X.\\
\midrule
Anorexia Nervosa & 3495 & 10982 & 28494655 & \\
Anxiety & 7016 & 14745 & 26754954 & \\
Bipolar Affective Disorder & 9412 & 137760 & 173062 & Biorxiv\\
Birth Weight & 153781 & NA & 27680694 & \\
Body mass index & 339224 & NA & 25673413 & \\
\addlinespace
Cannabis use, Lifetime & 32330 & NA & 27023175 & \\
Depressive Symptoms & 161460 & NA & 27089181 & \\
Education, Years & 293723 & NA & 27225129 & \\
Extraversion & 160713 & NA & 24828478 & \\
Neuroticism & 170911 & NA & 27089181 & \\
\addlinespace
Schizophrenia & 36989 & 113075 & 25056061 & \\
Subjective Well-being & 298420 & NA & 27089181 & \\
\bottomrule
\end{tabular}
\end{table}
n cases indicate number of cases in the discovery samples

Polygenic Scores were computed using ...

Association with schizophrenia is based on a binomial logistic regression of 2861 cases and 18843 controls - adjusting for age, gender, 10 principal components of genetic similarity and genotype wave. Pseudo\(R^{2}\) is calculated using Nagelkerke formula (`rcompanion::nargelkerke`s) - comparing the full model to the model without the PGS. p-values printed in indicate significance after correction for multiple testing (p 0.05/12).

``` r
source("~/trajectory/single/all/script/Nagelkerke.R")

dim(additional)
dt_pg <- merge(additional,dt,by = "pid")

dt_pg[!is.na(C1),.N,case]

pval <- lapply(12:23, function(x) {
  tmp = dt_pg 
  tmp$y = tmp[,x,with=F]
  summary(glm(!is.na(case)~y+gender+
                C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave, 
              data=tmp, family="binomial" ))[[13]][2,4]})

null <- glm(!is.na(case)~birthdate+gender+
              C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave, 
            data=dt_pg, family="binomial" )

nag <-  lapply(12:23, function(x) {
  tmp = dt_pg 
  tmp$y = tmp[,x,with=F]
  nagelkerke(fit = glm(!is.na(case)~y+birthdate+gender+
                         C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave, 
                       data=tmp, family="binomial" ),
             null = null)$Pseudo.R.squared.for.model.vs.null[3,]})
```

``` r
options(digits=3, scipen=T)
kable(data.frame(Phenotype= name, p= unlist(pval), R2= unlist(nag)),digits=c(1,32,4),,format="latex",booktabs=T)
```

\begin{tabular}{lrr}
\toprule
Phenotype & p & R2\\
\midrule
Anorexia Nervosa & 2.11e-02 & 0.0005\\
Anxiety & 5.39e-01 & 0.0000\\
Bipolar Affective Disorder & 7.61e-04 & 0.0011\\
Birth Weight & 4.29e-01 & 0.0001\\
Body mass index & 3.34e-02 & 0.0004\\
\addlinespace
Cannabis use, Lifetime & 3.69e-03 & 0.0008\\
Depressive Symptoms & 2.40e-12 & 0.0046\\
Education, Years & 7.66e-05 & 0.0014\\
Extraversion & 6.89e-03 & 0.0007\\
Neuroticism & 1.80e-12 & 0.0046\\
\addlinespace
Schizophrenia & 2.24e-29 & 0.0116\\
Subjective Well-being & 4.38e-11 & 0.0040\\
\bottomrule
\end{tabular}
### Rare Variants:

Computed by Andrea Ganna ...

### Merge MDS with risk variables:

``` r
data <- merge(dt_mds[,1:33,with=F],additional,by="pid")
dim(data)
```

    ## [1] 5432   91

``` r
ind_birth <- c(34:41,86:87) 
ind_mat_infek <- 42:43
ind_prs <- (44:55)[pval<.05/12]
ind_pc_wave <- 56:66
ind_sev <- 67:73
ind_infek <- 74:77
ind_rare <- 80:81
ind_parental <- 88:91
```

### Registry Variables:

#### Medical Birth Registry (MBR)

(see [documentation](https://sundhedsdatastyrelsen.dk/-/media/sds/filer/registre-og-services/nationale-sundhedsregistre/graviditet-foedsler-og-boern/foedselsregisteret/dokumentation-af-fdselsregisteret-1973-1996.xls?la=da):

APGAR 5: Appearance, Pulse, Grimace, Activity, Respiration at 5 minutes

``` r
table(data$apgar5)[table(data$apgar5)>4]
## 
##    0   10    5    6    7    8    9    A 
##   15 4978   10   18   34   85  221   51
data$apgar5=factor(data$apgar5,
                   levels=c("0","1","2","3","4","5",
                            "6","7","8","9","10","A"))
data$apgar5=as.character(data$apgar5)
data$apgar5[data$apgar5=="A"]=NA
data$apgar5=as.numeric(data$apgar5)
table(data$apgar5)[table(data$apgar5)>4]
## 
##    0    5    6    7    8    9   10 
##   15   10   18   34   85  221 4978
```

Birth Length:

``` r
#names(table(data$"Birth Length")[table(data$"Birth Length")>0])
leng=data$"Birth Length"
leng[leng %in%c("A","1","10")]=NA
leng=as.numeric(as.character(leng))
data$"Birth Length"=leng;rm(leng)
table(data$"Birth Length")[table(data$"Birth Length")>4]
## 
##  36  37  38  39  41  42  43  44  45  46  47  48  49  50  51  52  53  54 
##   5   5   7   5  11  14  14  23  35  66 122 252 397 807 799 905 752 540 
##  55  56  57  58  59 
## 298 154  62  27   5
```

Birth Weight: When calculating the value of bw\_score all singletons in the MBR born between 1981 and 2005 and having information of birth weight as well as gestational age were divided according to gender, gestational age in weeks, and for persons with gestational age of 28 weeks or more also into the following groups of calendar year at date of birth: 1981-1985, 1986- 1990, 1991-1995, 1996-2000, 2001-2005. The variable bw\_score contains the percentage of persons with same gender, gestational age, and calendar period (only if gestational age is 28 weeks or larger) who have the same or a smaller birth weight than the index person.

``` r
hist(data$"Birth Weight",na.rm = T,main = "Birth Weight Score")
## Warning in plot.window(xlim, ylim, "", ...): "na.rm" is not a graphical
## parameter
## Warning in title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...):
## "na.rm" is not a graphical parameter
## Warning in axis(1, ...): "na.rm" is not a graphical parameter
## Warning in axis(2, ...): "na.rm" is not a graphical parameter
```

![](Analysis_description_files/figure-markdown_github/BWS-1.pdf)

Maternal Smoking during pregnancy: From 1991 to 1996 it was registered as a binary (smoker of non-smoker). In the period 1997 to 2005 the register had more refined data, but this was transformed in a binary (non-smoker=non-smoker. smoker= smoker, stopped smoking during the first trimester, stopped smoking after the first trimester, smokes at most 5 cigarettes daily, smokes 6-10 cigarettes daily, smokes 11-20 cigarettes daily, smokes more than 20 cigarettes daily or smoker quantity unspecified).

``` r
table(data$B_RYGER)
## 
##   0   1 
## 565 468
names(table(data$C_RYGER))
## [1] "0"  "10" "20" "21" "22" "99"
data$C_RYGER[data$C_RYGER==99]=NA
table(data$C_RYGER)[table(data$C_RYGER)>4]
## 
##  0 21 22 
## 35  7 10
data$"Maternal Smoking in Pregnancy" <- data$B_RYGER
data$"Maternal Smoking in Pregnancy"[!is.na(data$C_RYGER)] <-1
table(data$"Maternal Smoking in Pregnancy")
## 
##   0   1 
## 565 525

data$wave[is.na(data$wave)]="Other"
data <- data[,wave:=factor(wave)]

ind_birth <- c(ind_birth[-c(7:8)],92); names(data)[ind_birth]
## [1] "Maternal Age"                  "Paternal Age"                 
## [3] "Birth Length"                  "apgar5"                       
## [5] "Gestational Age"               "Birth Weight"                 
## [7] "B_SECTIOU"                     "B_I11"                        
## [9] "Maternal Smoking in Pregnancy"
```

Number of Hospitalizations: All psychiatric hospital contacts with Schizophrenia as either main diagnosis (aktionsdiagnose/hoveddiagnose) or basic diagnosis (grundmorbus) in the Danish Psychiatric Central Research Register complete until December 31, 2016.

``` r
nms = grep("d2100",names(data))
nms <- names(data)[nms];nms
## [1] "d2100_ptype0_contacts" "d2100_ptype0_days"     "d2100_ptype1_contacts"
## [4] "d2100_ptype1_visits"   "d2100_ptype2_visits"   "d2100_ptype2_novisits"
## [7] "d2100_ptype3_contacts"
data <- data[,"Number of hospitalizations":= 
               Reduce('+',.SD), .SDcols=nms[c(1,3,7)]]
table(data$"Number of hospitalizations")[table(
  data$"Number of hospitalizations")>4]
## 
##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
## 1503  900  545  359  278  222  169  161  120  125   92   79   77   54   46 
##   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
##   60   39   42   36   31   30   25   30   22   19   23   23    7   13   11 
##   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44 
##   13    9   12   21   10   13    8   12    9    8    8    9    5    8    8 
##   45   49   51   54   60 
##   11    5    5    5    5
```

Total time hospitalized: Contains the total number of days admitted as 24 hour inpatient with a Schizophrenia diagnosis. For an admission date of the admission (unfinished admission) the admission is presumed ongoing until December 31, 2016. An exception for this rule is when the person has another contact in the Danish Psychiatric Central Research Register starting after the first unfinished admission (ptype = 0). In this case all unfinished admissions with ptype = 0 for this person are counted as having a duration of 1 day only. For persons with overlapping admissions with patient type 0 each day only counts once.

``` r
data <- data[,"Total time hospitalized" := d2100_ptype0_days]
table(data$"Total time hospitalized")[table(
  data$"Total time hospitalized")>4]
## 
##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
## 1919   44  117   57   51   31   28   31   36   21   15   27   22   23   25 
##   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
##   21   22   21   19   11   13   20   16   14   14   11   22   18   11   12 
##   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44 
##   18   10   17   12   12   16   15   20    5   19   16   13    9   12   19 
##   45   46   47   48   49   50   51   52   53   54   55   56   57   58   59 
##   19   13   16   11   21   15   13    7    8   11   12   17   17   10    9 
##   60   61   62   63   64   65   66   67   68   69   70   71   72   73   74 
##   12   11   19   11   14   13   12   12   15    7    7   10    8   10    8 
##   75   76   77   78   79   80   81   82   83   84   85   86   87   88   89 
##   12   13   17    7   11   11   13   10   10    9   18    9   14   12    9 
##   90   91   92   93   94   95   96   97   98   99  100  101  102  103  104 
##    8   10    9    6    9    6   12    6    5    8    9   15    5    9    7 
##  105  106  107  108  109  110  111  112  113  114  115  116  117  118  119 
##    8   13    8    9   10   10    8   11   10    9    7    7    5    9    9 
##  120  121  122  123  124  126  127  128  129  130  132  133  134  135  137 
##   10    6   12    9    9   10    8    9    6   11    7   11    8   10    6 
##  138  139  140  141  143  145  147  149  152  155  156  157  159  160  163 
##    9    5   10    6    7   10    9    7    9    5    8    5    5    6    7 
##  164  166  167  168  170  171  172  173  174  175  177  179  181  183  186 
##    5    8    8    5    5    5    7    7    8    6    9    8    6    5    8 
##  187  188  190  192  193  194  195  196  200  201  204  206  212  213  215 
##    6    7    6    9    8   10    7    6    5    5    6    6    5    7    6 
##  217  220  224  225  227  228  231  232  233  234  238  240  242  243  245 
##    6    9    5    5    6    5    6    5    6    6    6    6    5    5    7 
##  247  250  251  252  253  258  261  266  273  274  285  289  295  307  312 
##    5    5    7    7    6    5    5    5    5    7    6    7    5    6    7 
##  315  323  325  337  342  344  357  375  404  423  531  543 
##    7    5    5    5    5    5    5    5    5    5    6    5
ind_sev <- 93:94; names(data)[ind_sev]
## [1] "Number of hospitalizations" "Total time hospitalized"
```

From the National Patient Register, information an maternal infections during pregnancy was obtained. It was defined as a maternal diagnosis of infection (bacterial and viral) within a period of time prior to the date of birth that corresponded to the gestational age of the child.

``` r
table(data$"Maternal Infection during Pregnancy, Viral")
## 
##    0    1 
## 5412   20
table(data$"Maternal Infection during Pregnancy, Bacterial")
## 
##    0    1 
## 5270  162
```

### Standardize PRS scores

``` r
prs = names(data)[ind_prs]  
data[!is.na(C1),c(prs):= lapply(.SD,scale),.SDcols=prs]
```

### Association Analysis

Wrapper functions for `stats::manova` and `stats::anova`:

``` r
mancova.func <- function(tmp,ind_x,ind_y,mf){
tmp<- as.data.frame(tmp)
p.val <- Fstat <- rep(NA,length(ind_x))
for(p in (1:length(ind_x))){
  tmp$x=tmp[,ind_x[p]]
  tmp$y = as.matrix(tmp[,ind_y])
  manova.y=manova(mf, data = tmp[!is.na(tmp$x),]);
  p.val[p] = summary(manova.y)$stats[1,6]
  Fstat[p] = summary(manova.y)$stats[1,3]
  } 
names(p.val)=names(tmp)[ind_x]
names(Fstat)=names(tmp)[ind_x]
list(p.val=p.val,Fstat=Fstat) }

ancova.func <- function(tmp,ind_x,ind_y,mf){
tmp<- as.data.frame(tmp)
p.val <- betas <- array(NA,dim=c(length(ind_y),length(ind_x)))
for(p in (1:length(ind_x))){
  for(q in (1:length(ind_y))){
    tmp$x=tmp[,ind_x[p]]
    tmp$y = as.matrix(tmp[,ind_y[q]])
    lm = lm(mf, data = tmp[!is.na(tmp$x),])
    anova.y=anova(lm)
    p.val[q,p] = anova.y[1,5]
    betas[q,p] = coef(lm)[2]}} 
colnames(p.val)=names(tmp)[ind_x]
colnames(betas)=names(tmp)[ind_x]
list(p.val=p.val,beta=betas) }
```

Number of tests in genetic risk analysis:

``` r
n_test_gen <- length(c(ind_prs, ind_rare[1]));n_test_gen
```

    ## [1] 8

Number of tests in registry variable analysis:

``` r
n_test_regist <- length(c(ind_sev,ind_birth,ind_infek,
                          ind_mat_infek, ind_parental));n_test_regist
```

    ## [1] 21

Polygenic Scores:

``` r
ind_y=23:29; # Selecting Dim 1-3
ind_x=ind_prs;names(data)[ind_x]
## [1] "PGS - Bipolar Affective Disorder" "PGS - Cannabis use, Lifetime"    
## [3] "PGS - Depressive Symptoms"        "PGS - Education, Years"          
## [5] "PGS - Neuroticism"                "PGS - Schizophrenia"             
## [7] "PGS - Subjective Well-being"
mf = formula(y ~ x +birthdate+gender 
             +C1 +C2 +C3 +C4 +C5 +C6 +C7 +C8 +C9 +C10+wave)

manova.results_prs<-mancova.func(data[!is.na(data$C1),], ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_prs$p.val <= .05/n_test_gen]
anova.results_prs<-ancova.func(data[!is.na(data$C1),], ind_x_sig,ind_y,mf)
```

Association with rare mutations:

``` r
ind_x=ind_rare;names(data)[ind_x]
```

    ## [1] "Disruptive or Damaging Mutations" "Synonymous Mutations"

``` r
mf = formula(y ~ x + birthdate*gender+C1+C2+C3+C4+C5+C6+C7+C8+C8+C9+C10)

manova.results_rare<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_rare$p.val <= .05/n_test_gen]
anova.results_rare<-ancova.func(data, ind_x_sig,ind_y,mf)
```

Contacts and visits

``` r

ind_x=ind_sev;names(data)[ind_x]
## [1] "Number of hospitalizations" "Total time hospitalized"
mf = formula(y ~ x + birthdate*gender)

manova.results_sev<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_sev$p.val <= .05/n_test_regist]
anova.results_sev<-ancova.func(data, ind_x_sig,ind_y,mf)
```

Association with birth & pregnancy variables:

``` r

ind_x=c(ind_birth);names(data)[ind_x]
## [1] "Maternal Age"                  "Paternal Age"                 
## [3] "Birth Length"                  "apgar5"                       
## [5] "Gestational Age"               "Birth Weight"                 
## [7] "B_SECTIOU"                     "B_I11"                        
## [9] "Maternal Smoking in Pregnancy"
mf = formula(y ~ x + birthdate*gender)

manova.results_birth<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_birth$p.val <= .05/n_test_regist]
anova.results_birth<-ancova.func(data, ind_x_sig,ind_y,mf)
```

Association with infection diagnoses during pregancy

``` r
ind_x=ind_mat_infek;names(data)[ind_x]
## [1] "Maternal Infection during Pregnancy, Bacterial"
## [2] "Maternal Infection during Pregnancy, Viral"
summary(data[,ind_x])
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    42.0    42.2    42.5    42.5    42.8    43.0
mf = formula(y ~ x + birthdate*gender)

manova.results_mat_infek<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_mat_infek$p.val <= .05/n_test_regist]
#anova.results_mat_infek<-ancova.func(data, ind_x_sig,ind_y,mf)
```

Parental diagnosis of Schizophrenia and all Fxx diagnosis

``` r

ind_x=ind_parental;names(data)[ind_x]
## [1] "Maternal Diagnosis, Any Psychiatric"
## [2] "Maternal Diagnosis, Schizophrenia"  
## [3] "Paternal Diagnosis, Any Psychiatric"
## [4] "Paternal Diagnosis, Schizophrenia"
mf = formula(y ~ x + birthdate*gender)

manova.results_fam<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_fam$p.val <= .05/n_test_regist]
anova.results_fam<-ancova.func(data, ind_x_sig,ind_y,mf)
```

Association with infection diagnoses

``` r

ind_x=ind_infek;names(data)[ind_x]
## [1] "Otitis Infection"    "Bacterial Infection" "CNS Infection"      
## [4] "Viral Infection"
mf = formula(y ~ x + birthdate*gender)

manova.results_infek<-mancova.func(data, ind_x,ind_y,mf)

ind_x_sig = ind_x[manova.results_infek$p.val <= .05/n_test_regist]
anova.results_infek<-ancova.func(data, ind_x_sig,ind_y,mf)
```

Heatmap of associations using twerk of `corrplot::corrplot` included in `diagtraject::corrplot`:

``` r
anovas <- list(
  anova.results_prs,
  anova.results_rare,
  anova.results_birth,
  anova.results_infek,
  #anova.results_mat_infek,
  anova.results_sev,
  anova.results_fam)
betas <- t(Reduce('cbind',lapply(anovas,function(x) x$beta)))
p.vals <- t(Reduce('cbind',lapply(anovas,function(x) x$p.val)))

colnames(betas) <- colnames(p.vals) <- paste("Dimension",1:ncol(betas))  

library(corrplot)
load_all("~/trajectory/schizotracks190123/diagtraject0.6/")
corr.plot <- corrplot(betas[,1:3], is.corr=FALSE, p.mat=p.vals[,1:3],
                      method="p.square", insig="pch", 
                      sig.level=0.05/nrow(betas), pch="*", pch.cex=1.2, 
                      cl.ratio=0.25, cl.align.text="l", cl.lim=c(-6,6), 
                      tl.col="black", bg="#fffffc",ssl.pos = "r")
```

![](Analysis_description_files/figure-markdown_github/corrplot-1.pdf)

``` r
kable(
rbind(
as.data.frame(manova.results_prs),
as.data.frame(manova.results_rare),
as.data.frame(manova.results_prs),
as.data.frame(manova.results_birth),
as.data.frame(manova.results_infek),
as.data.frame(manova.results_mat_infek),
as.data.frame(manova.results_sev),
as.data.frame(manova.results_fam)
),format="latex",booktabs=T)
```

\begin{tabular}{lrr}
\toprule
  & p.val & Fstat\\
\midrule
PGS - Bipolar Affective Disorder & 0.327 & 1.154\\
PGS - Cannabis use, Lifetime & 0.560 & 0.833\\
PGS - Depressive Symptoms & 0.020 & 2.375\\
PGS - Education, Years & 0.005 & 2.921\\
PGS - Neuroticism & 0.023 & 2.333\\
\addlinespace
PGS - Schizophrenia & 0.000 & 4.322\\
PGS - Subjective Well-being & 0.391 & 1.054\\
Disruptive or Damaging Mutations & 0.000 & 7.005\\
Synonymous Mutations & 0.000 & 3.978\\
PGS - Bipolar Affective Disorder1 & 0.327 & 1.154\\
\addlinespace
PGS - Cannabis use, Lifetime1 & 0.560 & 0.833\\
PGS - Depressive Symptoms1 & 0.020 & 2.375\\
PGS - Education, Years1 & 0.005 & 2.921\\
PGS - Neuroticism1 & 0.023 & 2.333\\
PGS - Schizophrenia1 & 0.000 & 4.322\\
\addlinespace
PGS - Subjective Well-being1 & 0.391 & 1.054\\
Maternal Age & 0.000 & 8.090\\
Paternal Age & 0.001 & 3.712\\
Birth Length & 0.000 & 5.422\\
apgar5 & 0.187 & 1.433\\
\addlinespace
Gestational Age & 0.059 & 1.941\\
Birth Weight & 0.000 & 4.637\\
B\_SECTIOU & 0.626 & 0.754\\
B\_I11 & 0.534 & 0.864\\
Maternal Smoking in Pregnancy & 0.010 & 2.662\\
\addlinespace
Otitis Infection & 0.004 & 2.950\\
Bacterial Infection & 0.000 & 10.469\\
CNS Infection & 0.710 & 0.656\\
Viral Infection & 0.001 & 3.436\\
Maternal Infection during Pregnancy, Bacterial & 0.026 & 2.270\\
\addlinespace
Maternal Infection during Pregnancy, Viral & 0.302 & 1.194\\
Number of hospitalizations & 0.000 & 32.478\\
Total time hospitalized & 0.000 & 33.209\\
Maternal Diagnosis, Any Psychiatric & 0.000 & 10.267\\
Maternal Diagnosis, Schizophrenia & 0.001 & 3.652\\
\addlinespace
Paternal Diagnosis, Any Psychiatric & 0.000 & 5.932\\
Paternal Diagnosis, Schizophrenia & 0.000 & 4.301\\
\bottomrule
\end{tabular}
Replication
-----------

### Replication dataset:

Sequence Analysis is performed in non-overlapping dataset consisting of patients diagnosed January 1 2013 - December 31 2016.

``` r
load("diagdates2016.Rda")  
diagdates_rep <- data.table(rbind(diagdates,diagdates2016))
```

Sequence object is created as [above](#method) repeated for dataset of replication samples and primary cohort:

``` r
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
dt.m[,time:=time1]
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

Weighted MDS is performed assigning weights of `10^-20` to replication sample:

``` r
dt_rep[,mcor2:=dist_unique_OM_rep$mcor]
case_cor <- dt_rep[case==1,sum(case2012,na.rm = T)+1e-20,mcor2]
setkey(case_cor, mcor2)
w <- case_cor$V1
```

Showing that MDS of Case2012 in unaltered

``` r
wmd2012 <- wcmdscale(dist_unique_OM_rep$dist[w>=1,w>=1],w = w[w>=1], k=7) 
wmd_rep <- wcmdscale(dist_unique_OM_rep$dist,w = case_cor$V1, k=7,) 

#test
#wmd2012t<- wcmdscale(dist_OM, k=7) 
#head(wmd_rep[dt_rep$mcor2,][!is.na(dt_rep$case2012),])
#cor(wmd_rep[dt_rep$mcor2,][!is.na(dt_rep$case2012),],wmd2012t)

sum(abs(scale(wmd2012[,1])-scale(wmd_rep[w>=1,1])))
## [1] 6.01e-13

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
```

Thus producing a projection of Case2016 unto the MDS og Case2012:

### Merging with putative risk variable data

Testing which significant associations replicate:

``` r
data_rep <- merge(dt_rep,additional,by="pid")
```

Modify as [above](#Recode):

``` r
ind_birth <- c(34:41,86:87) 
ind_mat_infek <- 42:43
ind_prs <- (44:55)[pval<.05/12]
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

ind_birth <- c(ind_birth[-c(7:8)],92); names(data_rep)[ind_birth]
nms = grep("d2100",names(data_rep))
nms <- names(data_rep)[nms];nms
data_rep <- data_rep[,"Number of hospitalizations":= 
                       Reduce('+',.SD), .SDcols=nms[c(1,3,7)]]
table(data_rep$"Number of hospitalizations")[
  table(data_rep$"Number of hospitalizations")>4]
data_rep <- data_rep[,"Total time hospitalized" := d2100_ptype0_days]
table(data_rep$"Total time hospitalized")[
  table(data_rep$"Total time hospitalized")>4]
ind_sev <- 93:94; names(data_rep)[ind_sev]
table(data_rep$"Maternal Infection during Pregnancy, Viral")
table(data_rep$"Maternal Infection during Pregnancy, Bacterial")
```

### Standardize PRS scores

``` r
prs = names(data_rep)[grep("PGS",names(data_rep))]  
data_rep[!is.na(C1),c(prs):= lapply(.SD,scale),.SDcols=prs]
```

### Association analysis in replication data:

Test replication of associations with \(p < 0.05/17\):

``` r
rep_geno <- data.table(which(p.vals[1:3,1:3]<0.05/nrow(betas),
                             arr.ind = T),keep.rownames = T )
rep_geno[,dim:=paste0("Dim",col)]
rep_geno[,coef_main:=betas[1:3,1:3][
  which(p.vals[1:3,1:3]<0.05/nrow(betas))]]
rep_geno[,2:3:=NULL,with=F]


p_rep <- lapply(1:nrow(rep_geno), function(i) {  
  tmp <- data_rep[is.na(case2012)&!is.na(C1)]
  tmp$x  <-tmp[,rep_geno[i,rn],with=F]
  tmp$y <- tmp[,rep_geno[i,dim],with=F]
mf <- y~x+gender*birthdate+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave
list(anova= anova(lm(mf, tmp)),coef=coef(lm(mf, tmp)))
})

rep_geno[,N:=unlist(lapply(p_rep,function(x) sum(x[[1]][,1])+1))]
rep_geno[,Fstat:=unlist(lapply(p_rep,function(x) x[[1]][1,4]))]
rep_geno[,p.val:=unlist(lapply(p_rep,function(x) x[[1]][1,5]))]
rep_geno[,coef:=unlist(lapply(p_rep,function(x) x[[2]][2]))]
```

``` r
options(digits=3, scipen=T)
kable(rep_geno,digits=c(1,1,2,1,1,32,2),format="latex",booktabs=T)
```

\begin{tabular}{llrrrrr}
\toprule
rn & dim & coef\_main & N & Fstat & p.val & coef\\
\midrule
Disruptive or Damaging Mutations & Dim1 & -0.14 & 154 & 2.3 & 0.1320 & -0.46\\
PGS - Education, Years & Dim3 & 0.12 & 650 & 4.8 & 0.0284 & 0.16\\
PGS - Schizophrenia & Dim3 & -0.08 & 650 & 6.6 & 0.0102 & -0.20\\
\bottomrule
\end{tabular}
``` r
rep_regist <- data.table(which(p.vals[5:16,1:3]<0.05/nrow(betas),
                               arr.ind = T),keep.rownames = T )
rep_regist[,dim:=paste0("Dim",col)]
rep_regist[,coef_main:=betas[5:16,1:3][
  which(p.vals[5:16,1:3]<0.05/nrow(betas))]]
rep_regist[,2:3:=NULL,with=F]

p_rep <- lapply(1:nrow(rep_regist), function(i) {  
  tmp <- data_rep[is.na(case2012)]
  tmp$x  <-tmp[,rep_regist[i,rn],with=F]
  tmp$y <- tmp[,rep_regist[i,dim],with=F]
mf <- y~x+gender*birthdate
list(anova= anova(lm(mf, tmp)),coef=coef(lm(mf, tmp)))
})


rep_regist[,N:=unlist(lapply(p_rep,function(x) sum(x[[1]][,1])+1))]
rep_regist[,Fstat:=unlist(lapply(p_rep,function(x) x[[1]][1,4]))]
rep_regist[,p.val:=unlist(lapply(p_rep,function(x) x[[1]][1,5]))]
rep_regist[,coef:=unlist(lapply(p_rep,function(x) x[[2]][2]))]
```

``` r
kable(rep_regist,digits=c(1,1,3,1,1,32,3),format="latex",booktabs=T)
```

\begin{tabular}{llrrrrr}
\toprule
rn & dim & coef\_main & N & Fstat & p.val & coef\\
\midrule
Birth Length & Dim1 & -0.050 & 851 & 1.4 & 2.33e-01 & -0.011\\
Bacterial Infection & Dim1 & 0.920 & 870 & 3.2 & 7.25e-02 & 0.521\\
Viral Infection & Dim1 & 0.538 & 870 & 6.0 & 1.45e-02 & 0.505\\
Number of hospitalizations & Dim1 & 0.039 & 857 & 30.8 & 3.84e-08 & 0.085\\
Total time hospitalized & Dim1 & 0.001 & 857 & 2.9 & 8.77e-02 & 0.001\\
\addlinespace
Maternal Diagnosis, Any Psychiatric & Dim1 & 0.845 & 858 & 3.1 & 7.92e-02 & 0.307\\
Paternal Diagnosis, Any Psychiatric & Dim1 & 0.539 & 858 & 1.5 & 2.23e-01 & 0.146\\
Number of hospitalizations & Dim2 & 0.006 & 857 & 3.9 & 4.94e-02 & 0.032\\
Maternal Diagnosis, Schizophrenia & Dim2 & -0.739 & 858 & 1.0 & 3.25e-01 & 0.678\\
Paternal Diagnosis, Schizophrenia & Dim2 & -0.624 & 858 & 0.6 & 4.37e-01 & -0.429\\
\addlinespace
Maternal Age & Dim3 & 0.033 & 856 & 20.2 & 8.13e-06 & 0.040\\
Paternal Age & Dim3 & 0.021 & 843 & 2.5 & 1.12e-01 & 0.003\\
Birth Weight & Dim3 & 0.004 & 834 & 3.4 & 6.37e-02 & 0.003\\
Number of hospitalizations & Dim3 & -0.016 & 857 & 0.1 & 8.03e-01 & -0.009\\
Total time hospitalized & Dim3 & -0.001 & 857 & 0.6 & 4.25e-01 & -0.001\\
\addlinespace
Maternal Diagnosis, Any Psychiatric & Dim3 & -0.257 & 858 & 1.3 & 2.62e-01 & -0.120\\
Paternal Diagnosis, Any Psychiatric & Dim3 & -0.211 & 858 & 0.2 & 6.41e-01 & 0.081\\
\bottomrule
\end{tabular}
``` r
sum((rep_regist$coef_main<0) + (rep_regist$coef<0)!=1)
```

[1] 15

Robustness
----------

Selecting subset of data that does not require imputation:

``` r
#seq25 <- seq[as.numeric(dt_mds[case2012==1,year]) < 1991,]
seq <- seq[rownames(seq) %in% dt[case2012==1,pid],]
seq25 <- seq[as.numeric(dt_mds[case2012==1,year]) < 1991 & 
               as.numeric(dt_mds[case2012==1,SCZ]) < 25,]
seq25 <- seq25[,seq(1,25)]
```

Sequence Analysis with 3 different substitution cost settings: Jaccard Distance:

``` r
sub.cost <- seqsubm(seqdata = seq25, method = "CONSTANT")
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

1- Simple matching coefficient (Equivalent to Multichannel Sequence Analysis)

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

Computing dissimilarities

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

Multidimensional Scaling and MANCOVA:

``` r
MD_list <- mclapply(dist_list, function(x) cmdscale(x, k = 7), 
                    mc.cores = 16)

l <- length(MD_list)
m <- length(rownames(betas))
pc.wave.cor <- rownames(betas) %in% rep_geno$rn
setkey(data,pid)
data25 <-data[pid %in% rownames(seq25)] 


man_list <- lapply(1:m, function(x) lapply(1:l, function(y){
  tmp<- as.data.frame(data25)
    tmp$y1 <- MD_list[[y]]
    tmp$x1 <- tmp[,paste(rownames(betas))[x]]
   if(pc.wave.cor[x])   summary(manova(y1 ~ x1 + birthdate +gender+wave+
                                         C1+C2+C3+C4+C5+C6+C7+C8+C9+C10,
                                       data=tmp )) else
      summary(manova(y1 ~ x1+ birthdate +gender,data=tmp ))}))


p.mat <- Reduce("rbind",lapply(man_list, 
                               function(x) sapply(x, function(y) 
                                 y$stats[1,6])))

row.names(p.mat) <- rownames(betas)
colnames(p.mat) <- names(dist_list)
```

``` r
options(digits=2, scipen=T)
kable(p.mat,digits=40,format="latex",booktabs=T) %>% 
  kable_styling(latex_options = "scale_down")
```

\begin{table}[H]
\centering
\resizebox{\linewidth}{!}{
\begin{tabular}{lrrrrrrrrrrrrrrr}
\toprule
  & OM\_jacc\_.5 & OM\_jacc\_1 & HAM\_jacc & OM\_smc\_.5 & OM\_smc\_1 & HAM\_smc & OM\_const\_.5 & OM\_const\_1 & HAM\_const & eucl\_2 & chi2\_2 & eucl\_5 & chi2\_5 & eucl\_12 & chi2\_12\\
\midrule
PGS - Education, Years & 7.1e-02 & 5.6e-02 & 5.6e-02 & 6.7e-02 & 6.7e-02 & 6.7e-02 & 2.1e-01 & 1.5e-01 & 1.5e-01 & 2.2e-01 & 2.8e-02 & 6.9e-01 & 3.1e-02 & 1.0e-01 & 0.013643\\
PGS - Schizophrenia & 1.6e-02 & 1.7e-02 & 1.7e-02 & 1.8e-02 & 1.8e-02 & 1.8e-02 & 4.7e-02 & 7.2e-02 & 7.3e-02 & 1.9e-01 & 7.8e-01 & 8.4e-02 & 8.1e-01 & 3.8e-02 & 0.739486\\
Disruptive or Damaging Mutations & 6.5e-03 & 9.6e-03 & 9.6e-03 & 5.9e-02 & 5.9e-02 & 5.9e-02 & 7.1e-02 & 6.5e-02 & 6.5e-02 & 6.0e-03 & 2.5e-01 & 9.1e-03 & 2.1e-01 & 6.5e-02 & 0.297921\\
Synonymous Mutations & 8.2e-02 & 9.9e-02 & 9.9e-02 & 4.5e-02 & 4.5e-02 & 4.5e-02 & 2.0e-01 & 2.3e-01 & 2.3e-01 & 2.1e-01 & 7.3e-01 & 7.4e-02 & 6.1e-01 & 2.9e-01 & 0.556170\\
Maternal Age & 8.3e-05 & 1.2e-04 & 1.2e-04 & 1.2e-03 & 1.2e-03 & 1.2e-03 & 1.3e-04 & 6.3e-04 & 6.8e-04 & 1.4e-03 & 1.4e-01 & 3.2e-03 & 1.1e-01 & 1.1e-03 & 0.087149\\
\addlinespace
Paternal Age & 1.1e-01 & 1.9e-01 & 1.9e-01 & 1.3e-01 & 1.3e-01 & 1.3e-01 & 1.3e-01 & 2.3e-01 & 2.3e-01 & 2.9e-01 & 7.2e-01 & 2.2e-01 & 7.0e-01 & 2.6e-01 & 0.804683\\
Birth Length & 1.5e-03 & 1.8e-03 & 1.8e-03 & 6.1e-06 & 6.1e-06 & 6.1e-06 & 1.2e-01 & 2.5e-01 & 2.4e-01 & 3.5e-01 & 2.2e-01 & 3.3e-01 & 2.3e-01 & 5.8e-01 & 0.223223\\
Birth Weight & 4.7e-03 & 6.4e-03 & 6.4e-03 & 1.2e-01 & 1.2e-01 & 1.2e-01 & 9.3e-02 & 1.8e-01 & 1.8e-01 & 6.8e-02 & 1.2e-01 & 1.2e-01 & 8.8e-02 & 5.3e-02 & 0.115862\\
Bacterial Infection & 5.2e-13 & 1.3e-12 & 1.3e-12 & 2.8e-10 & 2.8e-10 & 2.8e-10 & 1.6e-10 & 1.4e-09 & 1.3e-09 & 1.7e-10 & 4.4e-01 & 1.5e-10 & 4.1e-01 & 3.8e-11 & 0.463429\\
Viral Infection & 1.1e-02 & 1.5e-02 & 1.5e-02 & 4.2e-02 & 4.2e-02 & 4.2e-02 & 1.9e-02 & 2.7e-02 & 2.7e-02 & 2.4e-02 & 5.9e-01 & 2.1e-02 & 6.0e-01 & 1.4e-02 & 0.509768\\
\addlinespace
Number of hospitalizations & 1.6e-15 & 3.8e-15 & 3.8e-15 & 6.8e-19 & 6.8e-19 & 6.8e-19 & 1.9e-12 & 1.6e-11 & 1.7e-11 & 1.3e-12 & 3.1e-01 & 1.0e-12 & 3.1e-01 & 5.4e-13 & 0.331074\\
Total time hospitalized & 1.3e-34 & 3.0e-33 & 3.0e-33 & 8.6e-25 & 8.6e-25 & 8.6e-25 & 5.2e-31 & 4.4e-28 & 5.6e-28 & 2.4e-27 & 8.4e-01 & 8.6e-26 & 8.5e-01 & 5.9e-29 & 0.795509\\
Maternal Diagnosis, Any Psychiatric & 5.6e-08 & 1.1e-08 & 1.1e-08 & 3.3e-09 & 3.3e-09 & 3.3e-09 & 2.4e-05 & 4.7e-06 & 4.5e-06 & 3.7e-06 & 5.7e-06 & 2.7e-06 & 4.7e-06 & 1.5e-05 & 0.000013\\
Maternal Diagnosis, Schizophrenia & 5.6e-05 & 4.6e-05 & 4.6e-05 & 2.0e-06 & 2.0e-06 & 2.0e-06 & 6.6e-04 & 8.7e-04 & 8.9e-04 & 3.6e-04 & 2.9e-01 & 2.0e-04 & 2.9e-01 & 6.3e-03 & 0.273661\\
Paternal Diagnosis, Any Psychiatric & 5.5e-05 & 2.2e-05 & 2.2e-05 & 1.1e-02 & 1.1e-02 & 1.1e-02 & 3.5e-05 & 3.6e-05 & 3.4e-05 & 6.7e-05 & 3.3e-01 & 6.2e-05 & 3.8e-01 & 1.6e-04 & 0.189496\\
Paternal Diagnosis, Schizophrenia & 2.3e-02 & 1.7e-02 & 1.7e-02 & 9.7e-02 & 9.7e-02 & 9.7e-02 & 2.4e-02 & 2.0e-02 & 1.9e-02 & 1.3e-02 & 5.5e-01 & 1.1e-02 & 5.1e-01 & 1.0e-02 & 0.445770\\
\bottomrule
\end{tabular}}
\end{table}
Post-hoc Regressions
--------------------

``` r
#data2 <- merge(data, dt_mds[,.(pid,age_min,n_dia)],by="pid")
data2 <- data
dim(data2)
```

    ## [1] 5432   94

Dimension 1 and Rare Mutations:

Adjusting for count of synonomous mutations:

``` r
anova(lm(Dim1~`Disruptive or Damaging Mutations`+`Synonymous Mutations`+
           gender+birthdate+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
         data=data2))$"Pr(>F)"[1]
```

    ## [1] 3.6e-09

Adjusting for Schizophrenia diagnosis before/after Jan 1, 1995.

``` r
anova(lm(Dim1~`Disruptive or Damaging Mutations`+`Synonymous Mutations`+
           gender+birthdate+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave+
           (format(age_min+birthdate,"%Y")>1995),data=data2))$"Pr(>F)"[1]
```

    ## [1] 3.6e-09

``` r
anova(lm(Dim1~`Disruptive or Damaging Mutations`+`Synonymous Mutations`+
           gender+birthdate+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave+
           (format(SCZ+birthdate,"%Y")>1995),data=data2))$"Pr(>F)"[1]
```

    ## [1] 3.6e-09

Adjusting for family history

``` r
anova(lm(Dim1~`Disruptive or Damaging Mutations`+gender+
           `Synonymous Mutations`+birthdate+
           `Maternal Diagnosis, Any Psychiatric`+
           `Maternal Diagnosis, Schizophrenia`+
           `Paternal Diagnosis, Any Psychiatric`+
           `Paternal Diagnosis, Schizophrenia`+
           C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,data=data2))$"Pr(>F)"[1]
```

    ## [1] 2.7e-09

Adjusting for age at first diagnosis and total number of diagnoses:

``` r
summary(lm(Dim1~`Disruptive or Damaging Mutations`+`Synonymous Mutations`+
             gender+birthdate+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
           data=data2))$coefficients[2,4]
```

    ## [1] 0.0056

``` r
summary(lm(Dim1~`Disruptive or Damaging Mutations`+age_min+
             `Synonymous Mutations`+gender+birthdate+
             C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
           data=data2))$coefficients[2,4]
```

    ## [1] 0.0066

``` r
summary(lm(Dim1~`Disruptive or Damaging Mutations`+n_dia+
             `Synonymous Mutations`+gender+birthdate+
             C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
           data=data2))$coefficients[2,4]
```

    ## [1] 0.23

Rare mutations ~ number of diagnoses:

``` r
summary(glm(`Disruptive or Damaging Mutations`~n_dia+gender+
              `Synonymous Mutations`+birthdate+
              C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
            data=data2,family = "poisson"))$coefficients[2,4]
```

    ## [1] 0.0054

``` r
exp(coef(glm(`Disruptive or Damaging Mutations`~n_dia+gender+
               `Synonymous Mutations`+birthdate+
               C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
             data=data2,family = "poisson")))[2]
```

    ## n_dia 
    ##  0.95

``` r
exp(confint(glm(`Disruptive or Damaging Mutations`~n_dia+
                  gender+`Synonymous Mutations`+birthdate+
                  C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
                data=data2,family = "poisson")))[2,]
```

    ##  2.5 % 97.5 % 
    ##   0.92   0.99

Dimension 3 and Polygenic Scores for Schizophrenia and Educational Attainment:

``` r
anova(lm(Dim3~`PGS - Education, Years`+gender+birthdate
         +C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
         data=data2))$"Pr(>F)"[1]
```

    ## [1] 0.0022

``` r
anova(lm(Dim3~`PGS - Schizophrenia`+gender+birthdate+
           C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
         data=data2))$"Pr(>F)"[1]
```

    ## [1] 0.0022

Substance Abuse ~ Education Polygenic scores:

``` r
summary(glm(!is.na(F10)~scale(`PGS - Education, Years`)+
              gender+birthdate+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
            data=data2,family = "binomial"))$coefficients[2,4]
```

    ## [1] 0.016

``` r
exp(coef(glm(!is.na(F10)~scale(`PGS - Education, Years`)+gender+birthdate+
               C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
             data=data2,family = "binomial")))[2]
```

    ## scale(`PGS - Education, Years`) 
    ##                            0.89

``` r
exp(confint(glm(!is.na(F10)~scale(`PGS - Education, Years`)+
                  gender+birthdate+C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
                data=data2,family = "binomial")))[2,]
```

    ##  2.5 % 97.5 % 
    ##   0.81   0.98

Substance Abuse ~ Schizophrenia Polygenic scores:

``` r
summary(glm(!is.na(F10)~`PGS - Schizophrenia`+gender+birthdate+
              C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
            data=data2,family = "binomial"))$coefficients[2,4]
```

    ## [1] 0.094

``` r
exp(coef(glm(!is.na(F10)~scale(`PGS - Schizophrenia`)+gender+birthdate+
               C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
             data=data2,family = "binomial")))[2]
```

    ## scale(`PGS - Schizophrenia`) 
    ##                          1.1

``` r
exp(confint(glm(!is.na(F10)~scale(`PGS - Schizophrenia`)+gender+birthdate+
                  C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
                data=data2,family = "binomial")))[2,]
```

    ##  2.5 % 97.5 % 
    ##   0.99   1.19

Early Substance Abuse ~ Schizophrenia Polygenic scores

``` r
data2[,median(F10,na.rm = T)]
data2[,earl_F10 := 0]
data2[F10<median(F10,na.rm = T) ,earl_F10 := 1]
```

``` r
summary(glm(earl_F10~`PGS - Schizophrenia`+gender+birthdate+
              C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
            data=data2,family = "binomial"))$coefficients[2,4]
```

    ## [1] 0.0083

``` r
exp(confint(glm(earl_F10~scale(`PGS - Schizophrenia`)+gender+birthdate+
                  C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+wave,
                data=data2,family = "binomial")))[2,]
```

    ##  2.5 % 97.5 % 
    ##    1.0    1.3

Dim3 ~ Maternal and Paternal Age

``` r
summary(lm(Dim3~`Maternal Age`+`Paternal Age`+gender+birthdate,
           data=data2))$coefficients[2:3,4]
```

    ## `Maternal Age` `Paternal Age` 
    ##        0.00076        0.26220

Dim3 ~ Maternal comparing to random sample:

``` r
summary(lm(`Maternal Age`~as.numeric(!is.na(case))+gender+birthdate,
           data=dt_pg[pid %in% dt_mds[Dim3<median(Dim3),pid]|random==1] ))
## 
## Call:
## lm(formula = `Maternal Age` ~ as.numeric(!is.na(case)) + gender + 
##     birthdate, data = dt_pg[pid %in% dt_mds[Dim3 < median(Dim3), 
##     pid] | random == 1])
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -14.553  -3.381  -0.288   3.093  22.421 
## 
## Coefficients:
##                            Estimate Std. Error t value Pr(>|t|)    
## (Intercept)              24.7252244  0.1205287  205.14   <2e-16 ***
## as.numeric(!is.na(case)) -0.8962732  0.1062991   -8.43   <2e-16 ***
## genderM                  -0.0564889  0.0583207   -0.97     0.33    
## birthdate                 0.0004220  0.0000139   30.33   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4.8 on 27132 degrees of freedom
##   (538 observations deleted due to missingness)
## Multiple R-squared:  0.04,   Adjusted R-squared:  0.0399 
## F-statistic:  377 on 3 and 27132 DF,  p-value: <2e-16

summary(lm(`Maternal Age`~as.numeric(!is.na(case))+gender+birthdate,
           data=dt_pg[pid %in% dt_mds[Dim3>quantile(Dim3,.75),pid]|
                        random==1] ))
## 
## Call:
## lm(formula = `Maternal Age` ~ as.numeric(!is.na(case)) + gender + 
##     birthdate, data = dt_pg[pid %in% dt_mds[Dim3 > quantile(Dim3, 
##     0.75), pid] | random == 1])
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -14.543  -3.364  -0.275   3.095  20.145 
## 
## Coefficients:
##                            Estimate Std. Error t value Pr(>|t|)    
## (Intercept)              24.7567024  0.1204704  205.50   <2e-16 ***
## as.numeric(!is.na(case)) -0.1423132  0.1310627   -1.09     0.28    
## genderM                  -0.0436566  0.0588447   -0.74     0.46    
## birthdate                 0.0004173  0.0000139   29.97   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4.8 on 26247 degrees of freedom
##   (104 observations deleted due to missingness)
## Multiple R-squared:  0.0341, Adjusted R-squared:  0.034 
## F-statistic:  309 on 3 and 26247 DF,  p-value: <2e-16
```

Dim1 ~ Infections

``` r
anova(lm(Dim1~`Bacterial Infection`+gender+birthdate,
         data=data2))$"Pr(>F)"[1]
```

    ## [1] 6.1e-07

``` r
anova(lm(Dim1~`Bacterial Infection`+gender+birthdate+
           `Number of hospitalizations`+`Total time hospitalized`,
         data=data2))$"Pr(>F)"[1]
```

    ## [1] 4.5e-07

``` r
anova(lm(Dim1~`Viral Infection`+gender+birthdate,
         data=data2))$"Pr(>F)"[1]
```

    ## [1] 0.000012

``` r
anova(lm(Dim1~`Viral Infection`+gender+birthdate+
           `Number of hospitalizations`+`Total time hospitalized`,
         data=data2))$"Pr(>F)"[1]
```

    ## [1] 9.6e-06

Excluding infections after onset of first psychiatric:

``` r
data2[,Bac_after := 0]
data2[`Bacterial Infection`==1, Bac_after := 
        as.numeric(age_Bacterial_infection>= age_min)]

data2[,Vir_after := 0]
data2[`Viral Infection`==1, Vir_after := 
        as.numeric(age_Viral_infection>= age_min)]
```

``` r
anova(lm(Dim1~`Bacterial Infection`+gender+birthdate,
         data=data2 [Bac_after!=1]))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dim1
    ##                         Df Sum Sq Mean Sq F value Pr(>F)    
    ## `Bacterial Infection`    1      6       6    0.38   0.54    
    ## gender                   1   1421    1421   96.92 <2e-16 ***
    ## birthdate                1  10030   10030  684.27 <2e-16 ***
    ## Residuals             4966  72789      15                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(lm(Dim1~`Viral Infection`+gender+birthdate,data=data2[Vir_after!=1]))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Dim1
    ##                     Df Sum Sq Mean Sq F value Pr(>F)    
    ## `Viral Infection`    1     90      90    6.03  0.014 *  
    ## gender               1   1696    1696  113.42 <2e-16 ***
    ## birthdate            1  10030   10030  670.75 <2e-16 ***
    ## Residuals         5221  78071      15                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
coef(lm(Dim1~`Viral Infection`+gender+birthdate,data=data2[Vir_after!=1]))
```

    ##       (Intercept) `Viral Infection`           genderM         birthdate 
    ##          -6.28210           0.19644          -0.92048           0.00097

Family history:

Browser
-------

Summarizing data

``` r
library(shiny)
library(ggplot2)
library(TraMineR)

mds<- read.csv(
"/data/projects/IBP/gonthe/sequence_analysis/all/output/mds_trajectories_k7_n31127_181023.csv",
               header=TRUE)
mds_scores <- mds[which(mds$case2012==1),]
mds_norm = as.data.frame(scale(mds_scores[,4:10]))
mds <- cbind(pid=mds_scores$pid, case=mds_scores$case, case2012=mds_scores$case2012, mds_norm)

load(
"/data/projects/IBP/gonthe/sequence_analysis/all/dist/seq_case-cohort_180314_thin.Rda")
seq <- seq[rownames(seq) %in% mds$pid,]

mds[,"pid"] <- rownames(mds) <- rownames(seq) <- sample(1:nrow((mds))) 

hclustlist <- lapply(1:7, function(x) lapply(
  x:7,function(y) hclust(dist(mds[,c(x+3,y+3)]),"ward.D2")))
tree <- lapply(hclustlist,function(x) lapply(x, cutree,k=1:100))
min <- lapply(tree,function(x) lapply(
  x,function(y) min(which(sapply(apply(y,2,table),min)<5))-1))
tree2 <- lapply(1:length(tree), function(x) lapply(
  1:length(tree[[x]]), function(y)  tree[[x]][[y]][,1:min[[x]][[y]]]))
tree3 <- lapply(1:length(tree2), function(x) lapply(
  1:length(tree2[[x]]), function(y) 
    tree2[[x]][[y]][!duplicated(
      tree2[[x]][[y]][,ncol(tree2[[x]][[y]])]),]))


rare <- lapply(alphabet(seq), function(x) apply(seq,1, function(y) 
  sum(y %in% x)))
nseq <- lapply(rare, function(x) lapply(x, function(y) sum(y>0)))
nseq2 <- lapply(nseq, function(x) sum(unlist(x)))
rare_states <- alphabet(seq)[which(nseq2<5)]
rare_seqs <- rare[which(nseq2<5)]
rare_list <- Reduce('+',rare_seqs)


comb <- levels(seq$a1)
comb[levels(seq$a1) %in% rare_states]<- "Other"
for(i in 1:36) levels(seq[[i]]) <- comb

attributes(seq)$alphabet <- unique(comb)
attributes(seq)$labels <- unique(comb)

seq<- seqdef(seq[,seq(3,36,3)])

rare_times <- lapply(seq, function(x) which(table(x)<4 & table(x)>0))

for(i in 1:length(rare_times))
seq[which(as.numeric(seq[,i]) %in% unlist(rare_times[i])),i] <- "Other"    

freq_bl <- lapply(1:length(tree), function(x) lapply(
  1:length(tree[[x]]), function(y)  lapply(1:min[[x]][[y]], function(z)
    seqstatd(seq[tree2[[x]][[y]][,min[[x]][[y]]]==z,],
             weighted=F)$Frequencies  )))

n_bl <- lapply(1:length(tree), function(x) lapply(
  1:length(tree[[x]]),function(y)  lapply(
    1:min[[x]][[y]], function(z)   sum(
      tree2[[x]][[y]][,min[[x]][[y]]]==z))))


cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
palette(cols)

alpa <- alphabet(seq)

setwd("shiny")

save(tree3,tree2, freq_bl, n_bl, cols,mds, alpa, file="data.Rda") 
rm(list=ls())
```

App Script:

``` r
library(shiny)
library(ggplot2)
library(TraMineR)
load("data.Rda")
source("seqmodst2.R")

ui <- fluidPage(
  headerPanel('Schizophrenia Trajectories'),
  sidebarPanel(
    selectInput('xcol', 'X Variable', names(mds[4:10])),
    selectInput('ycol', 'Y Variable', names(mds[4:10]),
                selected = names(mds[4:10])[[2]]),
    numericInput('clusters', 'Groups', 30,
                 min = 1, max = 300,  ),
    selectInput('type', 'Plot type', c("Modal state", "Chronogram"))
    
  ),
  mainPanel(
    p("All states observed in <5 individuals are labeled as", 
      span("other.",style="color:brown"), 
    "All states observed in <4 individuals at a 
    given age are also labeled as",
    span("other.",style="color:brown")),
    verbatimTextOutput("info"),
    plotOutput('plot1', click = "plot_click"),
    plotOutput('plot2',height = 200),
    plotOutput('plot3')
  )
)

server <- function(input, output) {
  
  selectedData <- reactive({
    mds[, c(input$xcol, input$ycol)]
  })
  
  clusters <- reactive({
  x <- which(names(mds[4:10]) %in% input$xcol)
  y <- which(names(mds[4:10]) %in% input$ycol)
  ind1 <- min(x,y)
  ind2 <- abs(x-y)+1
  if(input$clusters>ncol(tree3[[ind1]][[ind2]])) F else 
    tree2[[ind1]][[ind2]][,input$clusters]
    })

  df <-reactive({
    dft <- cbind(selectedData(),factor(clusters()))
  colnames(dft) <- c("x","y","cl")
  dft})

  
  output$info <- renderPrint({
    x <- which(names(mds[4:10]) %in% input$xcol)
    y <- which(names(mds[4:10]) %in% input$ycol)
    ind1 <- min(x,y)
    ind2 <- abs(x-y)+1
    if(clusters()==F) cat("maximum number of clusters for 
                          these dimensions is: ", 
                          paste(ncol(tree3[[ind1]][[ind2]])))
    })

  click_saved <- reactiveValues(singleclick =NULL)
  observeEvent(eventExpr = input$plot_click, 
               handlerExpr = { click_saved$singleclick <- 
                                 input$plot_click })
  np <-  reactive({
   # if(!is.null(click_saved$singleclick))
    nearPoints(df(),click_saved$singleclick, threshold = 100,
               maxpoints = 1,
               addDist = TRUE) #else nearPoints(df(),list(x=0,y=0,col=1),
    #threshold = 100,maxpoints = 1,addDist = TRUE)
  })  
  
  output$plot1 <- renderPlot({
    par(mar = c(5.1, 4.1, 0, 1))
    g <- ggplot(df(),aes(x=x,y= y,col = cl))+
    geom_point()+
    theme(legend.position="none")
    #if(nrow(np())<1) np <- matrix(rep(0,4),nrow = 1) else
    np <- np() 
    print(np)
    
    #if(nrow(np)>0){
      cla <-factor(clusters())
      np1 <- cla[which(rownames(mds) %in% rownames(np))] 
      print(np1)
      if(nrow(np())<1) cl <- cla %in% 1 else cl <- cla %in% np1      
      df <- df()
      g <- g+stat_ellipse(data=df[cl,])
    #}
    g
    })



  cl <- reactive({
    if(nrow(np())<1) df()$cl[which(df()$cl %in% 1)][1] else 
      df()$cl[which(rownames(mds) %in% rownames(np()))][1]
  })
  
  seq1 <- reactive({
    x <- which(names(mds[4:10]) %in% input$xcol)
    y <- which(names(mds[4:10]) %in% input$ycol)
    ind1 <- min(x,y)
    ind2 <- abs(x-y)+1
    c <-which(tree3[[ind1]][[ind2]][,input$clusters]==cl())
    #freq_x <- Reduce('+',freq_bl[[ind1]][[ind2]][c])/sum(c)
    freq_x <- Reduce('+',lapply(c, function(x) 
      freq_bl[[ind1]][[ind2]][[x]]* n_bl[[ind1]][[ind2]][[x]]))/
      Reduce('+',n_bl[[ind1]][[ind2]][c])
    n_x <- Reduce('+',n_bl[[ind1]][[ind2]][c])
    print("cp")
    cpal <- rep("grey", length(alpa))
    stat<-seqstatl(seqmodst2(freq_x,n=n_x,weighted=F))
    cpal[alpa %in% stat] <-  c(cols[1:length(stat)])
    cpal[alpa=="None"]<-"lightgray" 
    cpal[alpa=="Other"]<-"brown" 
    cpal[alpa=="censored"]<-"white"
    list(freq=freq_x,n=n_x,stat=stat,cpal=cpal)
    
  })
  
  output$plot2 <- renderPlot({
    par(mar = c(5.1, 4.1, 0, 1))
    if(input$type!="Chronogram") plot(seqmodst2(seq1()$freq,
                                                n=seq1()$n,weighted=F,
                                                freq.mean=T),
                                         cpal=seq1()$cpal) else{
      res<- list(seq1()$freq,n=seq1()$n)
    names(res) <- c("Frequencies", "ValidStates")
      class(res) <- c("stslist.statd", "list")
      attr(res, "nbseq") <- seq1()$n
      attr(res, "cpal") <- seq1()$cpal
      attr(res, "xtlab") <- colnames(seq1()$freq)
      attr(res, "weighted") <- F
      plot(res,weighted=F)}

      
  })
  output$plot3 <- renderPlot({
    par(mar = c(5.1, 4.1, 0, 1))
  if(input$type!="Chronogram"){   alp <- seq1()$stat
    col1 <- seq1()$cpal[alpa %in% seq1()$stat] }else{
      alp <- alpa[apply(seq1()$freq,1,max)>.1]
      col1 <- seq1()$cpal[apply(seq1()$freq,1,max)>.1]}
    plot(0,type='n',axes=FALSE,ann=FALSE)
    legend( 1,1, alp ,fill = col1,cex=2)
    
  })
  
}

shinyApp(ui = ui, server = server)
```
