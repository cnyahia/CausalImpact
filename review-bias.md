Measuring the Difference in Scores between Single and Double-Blind Reviews
================
Cesar Yahia
August 15, 2020

The single-blind review process is biased. Reviews are influenced by the author’s gender, institution, academic status (rank/positions held), reputation, etc.

-   Do double-blind reviews make the review process more equitable?

-   How can we test that empirically?

### Randomized Experiment Approach

To assess the efficacy of double-blind reviews, we can run a *randomized* experiment where randomly selected papers are given double-blind reviews (the remaining reviews are single-blind). Then, we can use a *two sample t-test (independent samples)* to evaluate the difference in review scores between the two groups (double and single). We can also stratify the authors in each group according to factors such as the author's country of origin/academic rank/gender and compare the scores between double/single within each stratified class.

#### Assumptions:

1.  Both groups (double and single) receive reviewers with similar features (where a reviewer's features refers to the tendency to give high/low scores). If the number of papers analyzed is large enough, this assumption is justified; precisely, due to randomization, we will have a balanced distribution of reviewers between the two groups.

#### Challenges:

1.  Who gets assigned single or double-blind reviews? Authors might not like being “randomly” assigned into one category or the other.

### Observational Study using Author’s Preferences for Single or Double-Blind Reviews

To avoid challenges resulting from the randomized assignment of double-blind reviews, we could use an *opt-in process* where authors select single or double-blind reviews for their papers.

#### Challenges:

1.  To perform an appropriate causal analysis we need to control for confounding variables that make it difficult to determine whether there is a bias in the review process or the papers are indeed of lower quality. For example, authors that select double-blind may have other characteristics (lack of mentorship, insufficient funding, other inequities earlier in the pipeline) that render their papers of lower quality. To isolate the impact of biased reviews, we need to control for such confounders.

#### Approach

We can first implement *matching* to control for confounders and then infer the causal impact of the double-blind review process.

In brief, for every subject in the double-blind group, matching finds a corresponding subject in the single-blind group that has similar characteristics. We refer to those characteristics as covariates/confounders (these covariates can be the author’s country of origin/academic rank/gender). Then, we create a list of matched pairs, and we infer the causal impact of double-blind reviews by applying a *paired t-test* between the double-blind review score and the single-blind review score.

To find the impact of double-blind reviews on subsets of authors with the same feature (e.g., all authors with the same country of origin), we can isolate the matched pairs that have the desired feature and then infer the causal impact for that specific group.

#### Assumptions:

1.  The tendency of reviewers to give high/low scores now has a greater impact on the results! Since we are implementing a paired t-test between matched pairs, the reviewers must have similar tendencies for every pair!

------------------------------------------------------------------------

------------------------------------------------------------------------

------------------------------------------------------------------------

In what follows, additional description on implementing matching is provided. The analysis below uses information from Coursera's course on causal inference.

The lalonde data in the Matchit library is used. In this case, the objective is to analyze the impact of labor training on income, where some confounding variables (such as age, education, marital status) influence both who gets training and income. We can control for those confounding variables by matching to isolate the impact of labor training.

#### Load Data

Load the lalonde data in the MatchIt library. Treat is a binary variable that takes the value 1 if the subject received labor training. The variable re78 is the outcome of interest representing income. All other variables are confounding variables.

``` r
data(lalonde, package = 'MatchIt')
# View(lalonde)
print(head(lalonde, 10))
```

    ##       treat age educ black hispan married nodegree re74 re75       re78
    ## NSW1      1  37   11     1      0       1        1    0    0  9930.0460
    ## NSW2      1  22    9     0      1       0        1    0    0  3595.8940
    ## NSW3      1  30   12     1      0       0        0    0    0 24909.4500
    ## NSW4      1  27   11     1      0       0        1    0    0  7506.1460
    ## NSW5      1  33    8     1      0       0        1    0    0   289.7899
    ## NSW6      1  22    9     1      0       0        1    0    0  4056.4940
    ## NSW7      1  23   12     1      0       0        0    0    0     0.0000
    ## NSW8      1  32   11     1      0       0        1    0    0  8472.1580
    ## NSW9      1  22   16     1      0       0        0    0    0  2164.0220
    ## NSW10     1  33   12     0      0       1        0    0    0 12418.0700

#### Data Analysis

First get covariates and create table one

``` r
xvars<-c("age","educ","black","hispan","married","nodegree","re74","re75")
table1<- CreateTableOne(vars=xvars,strata="treat", data=lalonde, test=FALSE)
```

Now get the standardize differences (smd) between covariates in treated and control group. Observe high smd values, greater than 0.2, that indicate significant differences in characteristics between the treated and control groups

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

#### Greedy Matching

Implement greedy matching using the Mahalanobis distance. The Mahalanobis distance between covariates *X*<sub>*i*</sub> for subject *i* and covariates *X*<sub>*j*</sub> for subject *j* is given by $d=\\sqrt{(X\_{i} - X\_{j})^{T} S^{-1} (X\_{i} - X\_{j})}$, where *S* is the covariance matrix associated with the covariates.

``` r
greedymatch<-Match(Tr=lalonde$treat,M=1,X=lalonde[xvars],replace=FALSE)
matched<-lalonde[unlist(greedymatch[c("index.treated","index.control")]), ]
```

Now let's look at the matched pairs, observe that we have much lower smd values indicating good balance between covariates in treated and control groups! We have 185 matched pairs.

``` r
matchedtab1<-CreateTableOne(vars=xvars, strata ="treat", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)
```

    ##                       Stratified by treat
    ##                        0                 1                 SMD   
    ##   n                        185               185                 
    ##   age (mean (SD))        24.21 (9.55)      25.82 (7.16)     0.190
    ##   educ (mean (SD))       10.23 (2.37)      10.35 (2.01)     0.052
    ##   black (mean (SD))       0.43 (0.50)       0.84 (0.36)     0.943
    ##   hispan (mean (SD))      0.06 (0.24)       0.06 (0.24)    <0.001
    ##   married (mean (SD))     0.20 (0.40)       0.19 (0.39)     0.027
    ##   nodegree (mean (SD))    0.69 (0.46)       0.71 (0.46)     0.035
    ##   re74 (mean (SD))     2681.77 (4754.79) 2095.57 (4886.62)  0.122
    ##   re75 (mean (SD))     1523.69 (2810.24) 1532.06 (3219.25)  0.003

Let's do some outcome analysis

``` r
inc_trt_m <- matched$re78[matched$treat==1]
inc_con_m <- matched$re78[matched$treat==0]
# diffinc_m <- inc_trt_m - inc_con_m
# paired t-test
t.test(inc_trt_m, inc_con_m, paired=TRUE)
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  inc_trt_m and inc_con_m
    ## t = 1.2523, df = 184, p-value = 0.212
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -473.2129 2117.9775
    ## sample estimates:
    ## mean of the differences 
    ##                822.3823

*observe that, for matched pairs, the impact of the labor training intervention is a net gain of $822.*

*previously (without matching) we got that the trained group had a net income loss of $635!!*

*this illustrates the importance of matching in controlling for confounders!*

Note however that the p-value is high indicating that the difference in means between treated and control groups is not siginificant.
