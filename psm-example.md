Inferring Causality using Propensity Score Matching!
================
Cesar Yahia
August 12, 2020

### Overview

This markdown file is used for evaluating the **causal impact** of an intervention (treatment) on subjects. In particular, the objective is to determine the impact of labor training on post-intervention earning levels (income).

To estimate the causal effect, we have to control for **confounders** that may bias the output. Specifically, we must control for variables such as age, education, race, marital status that influence the earnings regardless of the intervention.

To control for confounders, we implement **propensity score matching** (psm). The psm process invloves fitting a logistic regression (where the output is binary treatment status) such that the logistic regression gives us the probability of treatment for each subject. The propensity score is defined as the probability of treatment given the covariates. Thus, subjects with the same propensity score are equally likely to have been treated. This implies that if we isolate a specific propensity score level, then we are effectively considering subjects that have the same rate of treatment. In other words, for a fixed propensity score level, the distribution of covariates in the treated group is the same as the distribution of covariates in the control group (everyone has the same probability of treatment and so the control group will *not* be skewed towards subjects with specific covariates). Equivalently, matching on propensity scores achieves *balance* between covariate distributions across treated and control groups. Since treated and control groups have the same covariate distribution (i.e., same distribution of confounders), we can effectively infer the causal impact by analyzing outcome differences across the two groups. In contrast, if the distribution of confounders is different between treated and control groups, we would not know if the observed differences result from the intervention or from other confounding variables.

### Load Data

The lalonde data is in the MatchIt library

``` r
data(lalonde, package = 'MatchIt')
View(lalonde)
```

### Data Analysis

First get covariates and create table one

``` r
xvars<-c("age","educ","black","hispan","married","nodegree","re74","re75")
table1<- CreateTableOne(vars=xvars,strata="treat", data=lalonde, test=FALSE)
```

Now get the smd

``` r
print(table1,smd=TRUE)
```

    ##                       Stratified by treat
    ##                        0                 1                 SMD   
    ##   n                        429               185                 
    ##   age (mean (SD))        28.03 (10.79)     25.82 (7.16)     0.242
    ##   educ (mean (SD))       10.24 (2.86)      10.35 (2.01)     0.045
    ##   black (mean (SD))       0.20 (0.40)       0.84 (0.36)     1.668
    ##   hispan (mean (SD))      0.14 (0.35)       0.06 (0.24)     0.277
    ##   married (mean (SD))     0.51 (0.50)       0.19 (0.39)     0.719
    ##   nodegree (mean (SD))    0.60 (0.49)       0.71 (0.46)     0.235
    ##   re74 (mean (SD))     5619.24 (6788.75) 2095.57 (4886.62)  0.596
    ##   re75 (mean (SD))     2466.48 (3292.00) 1532.06 (3219.25)  0.287

Get the income for treated and untreated subjects

``` r
inc_trt <- lalonde$re78[lalonde$treat==1]
inc_con <- lalonde$re78[lalonde$treat==0]
```

Get the difference in mean income between treated and untreated subjects

``` r
mean_diff <- mean(inc_trt) - mean(inc_con)
print(mean_diff)
```

    ## [1] -635.0262

### Propensity Score Matching!

#### Estimate propensity scores using logistic regression

First start by fitting a logistic regression where the outcome is treatment (binary); this will give us the propensity scores!

``` r
psmodel<-glm(treat~age+educ+black+hispan+married+nodegree+re74+ re75,
             family=binomial(),data=lalonde)
#show coefficients
summary(psmodel)
```

    ## 
    ## Call:
    ## glm(formula = treat ~ age + educ + black + hispan + married + 
    ##     nodegree + re74 + re75, family = binomial(), data = lalonde)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.7645  -0.4736  -0.2862   0.7508   2.7169  
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -4.729e+00  1.017e+00  -4.649 3.33e-06 ***
    ## age          1.578e-02  1.358e-02   1.162  0.24521    
    ## educ         1.613e-01  6.513e-02   2.477  0.01325 *  
    ## black        3.065e+00  2.865e-01  10.699  < 2e-16 ***
    ## hispan       9.836e-01  4.257e-01   2.311  0.02084 *  
    ## married     -8.321e-01  2.903e-01  -2.866  0.00415 ** 
    ## nodegree     7.073e-01  3.377e-01   2.095  0.03620 *  
    ## re74        -7.178e-05  2.875e-05  -2.497  0.01253 *  
    ## re75         5.345e-05  4.635e-05   1.153  0.24884    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 751.49  on 613  degrees of freedom
    ## Residual deviance: 487.84  on 605  degrees of freedom
    ## AIC: 505.84
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
#get propensity score
pscore<-psmodel$fitted.values
```

get the minimum and maximum propensity scores

``` r
print(min(pscore))
```

    ## [1] 0.009080193

``` r
print(max(pscore))
```

    ## [1] 0.8531528

#### Match based on propensity scores

Next, we carry out the matching using the Match function. The matching process matches treated subjects with untreated subjects that have the same propensity scores. In the end, the treated group will have the same covariate distribution as the control group.

First, start with matching without a caliper

``` r
set.seed(931139)
psmatch<-Match(Tr=lalonde$treat,M=1,X=pscore,replace=FALSE)
```

Now get the matched data and create standardized differences (compare with previous standardized differences!)

``` r
matched<-lalonde[unlist(psmatch[c("index.treated","index.control")]), ]
xvars<-c("age","educ","black","hispan","married","nodegree","re74","re75")
matchedtable1<- CreateTableOne(vars=xvars,strata="treat", data=matched, test=FALSE)
print(matchedtable1, smd = TRUE)
```

    ##                       Stratified by treat
    ##                        0                 1                 SMD   
    ##   n                        185               185                 
    ##   age (mean (SD))        25.29 (10.65)     25.82 (7.16)     0.058
    ##   educ (mean (SD))       10.55 (2.71)      10.35 (2.01)     0.084
    ##   black (mean (SD))       0.47 (0.50)       0.84 (0.36)     0.852
    ##   hispan (mean (SD))      0.21 (0.41)       0.06 (0.24)     0.453
    ##   married (mean (SD))     0.20 (0.40)       0.19 (0.39)     0.027
    ##   nodegree (mean (SD))    0.65 (0.48)       0.71 (0.46)     0.127
    ##   re74 (mean (SD))     2351.12 (4192.62) 2095.57 (4886.62)  0.056
    ##   re75 (mean (SD))     1605.02 (2601.68) 1532.06 (3219.25)  0.025

Now redo the matching using a caliper!

``` r
set.seed(931139)
psmatchC<-Match(Tr=lalonde$treat,M=1,X=pscore,replace=FALSE, caliper=.1)
```

Now get the matched data and create standardized differences (compare with previous standardized differences!)

``` r
matchedC<-lalonde[unlist(psmatchC[c("index.treated","index.control")]), ]
xvars<-c("age","educ","black","hispan","married","nodegree","re74","re75")
matchedtable1C<- CreateTableOne(vars=xvars,strata="treat", data=matchedC, test=FALSE)
print(matchedtable1C, smd = TRUE)
```

    ##                       Stratified by treat
    ##                        0                 1                 SMD   
    ##   n                        111               111                 
    ##   age (mean (SD))        26.27 (11.10)     26.22 (7.18)     0.006
    ##   educ (mean (SD))       10.37 (2.66)      10.25 (2.31)     0.047
    ##   black (mean (SD))       0.72 (0.45)       0.74 (0.44)     0.040
    ##   hispan (mean (SD))      0.11 (0.31)       0.10 (0.30)     0.029
    ##   married (mean (SD))     0.24 (0.43)       0.24 (0.43)    <0.001
    ##   nodegree (mean (SD))    0.66 (0.48)       0.65 (0.48)     0.019
    ##   re74 (mean (SD))     2704.56 (4759.89) 2250.49 (5746.14)  0.086
    ##   re75 (mean (SD))     1969.10 (3169.08) 1222.25 (3081.19)  0.239

### Causal Impact Analysis

Now we use the matched data (with caliper=0.1) to determine the causal impact of the intervention.

``` r
inc_trt_matched <- matchedC$re78[matchedC$treat==1]
inc_con_matched <- matchedC$re78[matchedC$treat==0]
mean_diff_matched<- mean(inc_trt_matched) - mean(inc_con_matched)
print(mean_diff_matched)
```

    ## [1] 1246.806

*observe that for matched pairs the impact of the labor training intervention is a net gain of $1246 but without matching we obtained that the trained group had a net loss of income of 635 dollars*

*this illustrates the importance of matching in controlling for confounders!*

carry out a paired t-test

``` r
diff_inc <- inc_trt_matched - inc_con_matched
t.test(diff_inc)
```

    ## 
    ##  One Sample t-test
    ## 
    ## data:  diff_inc
    ## t = 1.4824, df = 110, p-value = 0.1411
    ## alternative hypothesis: true mean is not equal to 0
    ## 95 percent confidence interval:
    ##  -420.0273 2913.6398
    ## sample estimates:
    ## mean of x 
    ##  1246.806
