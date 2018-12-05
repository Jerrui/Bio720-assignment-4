---
title: "12.4 assignment"
author: "Jerry"
date: "2018Äê12ÔÂ4ÈÕ"
output: 
  html_document: 
    keep_md: yes
---

2.

```r
p_t <- function(wAA,wAa,waa,p0,n) {
    p <- rep(NA, n)  
    w_bar <- rep(NA, n)
   	p[1] <- p0 
	w_bar[1] <- p[1]*p[1]*wAA + 2*p[1]*(1-p[1])*wAa + (1-p[1])*(1-p[1])*waa
	for ( i in 2:n) {
		w_bar[i - 1] <- p[i-1]*p[i-1]*wAA + 2*p[i-1]*(1-p[i-1])*wAa + (1-p[i-1])*(1-p[i-1])*waa
		p[i] <- p[i-1]*p[i-1]*wAA/w_bar[i - 1] + p[i-1]*(1-p[i-1])*wAa/w_bar[i - 1]
	}
if (any(p > 0.9999)) {
    p.fixation <- min(which.max(p > 0.9999))
    cat("fixation for A occurs approximately at generation:", p.fixation )	
} else if (any((1-p) > 0.9999)) {
    q.fixation <- min(which.max((1-p) > 0.9999))
    cat("fixation for a occurs approximately at generation:", q.fixation )	
    }
  else {
        maxAlleleFreq <- max(p)
    	cat("fixation of A1 does not occur, max. allele frequency is:", print(maxAlleleFreq, digits = 2))
    }
plot(x = 1:n, y = p, 
     xlab="generations", 
	 ylab="Allele frequency (p)", main = "generations VS Allele frequency (p)",
	 pch = 20, col = "red", cex.lab = 1.5)   
    return(p)
}
p_t(1,0.9,0.8,0.3,100)
```

```
## fixation for A occurs approximately at generation: 96
```

![](12345_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```
##   [1] 0.3000000 0.3244186 0.3497597 0.3759023 0.4027082 0.4300248 0.4576887
##   [8] 0.4855293 0.5133734 0.5410490 0.5683902 0.5952402 0.6214553 0.6469071
##  [15] 0.6714845 0.6950951 0.7176653 0.7391401 0.7594826 0.7786726 0.7967049
##  [22] 0.8135880 0.8293416 0.8439952 0.8575859 0.8701572 0.8817568 0.8924355
##  [29] 0.9022460 0.9112417 0.9194759 0.9270011 0.9338684 0.9401270 0.9458240
##  [36] 0.9510042 0.9557099 0.9599805 0.9638533 0.9673627 0.9705407 0.9734168
##  [43] 0.9760182 0.9783702 0.9804956 0.9824155 0.9841491 0.9857140 0.9871262
##  [50] 0.9884003 0.9895495 0.9905858 0.9915201 0.9923623 0.9931214 0.9938055
##  [57] 0.9944219 0.9949772 0.9954775 0.9959281 0.9963339 0.9966995 0.9970287
##  [64] 0.9973251 0.9975920 0.9978323 0.9980487 0.9982435 0.9984189 0.9985769
##  [71] 0.9987190 0.9988470 0.9989622 0.9990659 0.9991592 0.9992432 0.9993189
##  [78] 0.9993869 0.9994482 0.9995034 0.9995530 0.9995977 0.9996379 0.9996741
##  [85] 0.9997067 0.9997360 0.9997624 0.9997862 0.9998075 0.9998268 0.9998441
##  [92] 0.9998597 0.9998737 0.9998863 0.9998977 0.9999079 0.9999171 0.9999254
##  [99] 0.9999329 0.9999396
```

3.

```r
allele_drift <- function(fA,alleles,n){
A_freq <- rep(NA,n)
a_freq <- rep(NA,n)

A_freq[1] <- fA
a_freq[1] <- (1-fA)
for (i in 2:n) {
allele_counts <- sample(c("a", "A"), size = alleles, replace = TRUE, prob = c(a_freq[i-1],A_freq[i-1]))
a_freq[i] <- as.vector(table(allele_counts)/length(allele_counts))[1]
A_freq[i] <- (1-(a_freq[i]))
}
plot(a_freq~c(1:n), type="l",
     xlab="generations", 
	 ylab="Allele frequency(A)", main = "genetic drift",
	 xlim = c(1,100),
   ylim = c(0,1), 
	 pch = 20, col = "red", cex.lab = 1.5)
lines(A_freq~c(1:n),col = "blue")
}
allele_drift(0.5,400,1000)
```

![](12345_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

4.

```r
allele_drift2 <- function(fA,alleles,n){
A_freq <- rep(NA,n)
a_freq <- rep(NA,n)
A_freq[1] <- fA
a_freq[1] <- (1-fA)
for (i in 2:n) {
allele_counts <- sample(c("a", "A"), size = alleles, replace = TRUE, prob = c(a_freq[i-1],A_freq[i-1]))
a_freq[i] <- as.vector(table(allele_counts)/length(allele_counts))[1]
A_freq[i] <- 1-(a_freq[i])
}
return(A_freq)
}
allele_counts_4 <- function(fA,alleles,n,x){
  repl <- replicate(x, allele_drift2(fA,alleles,n))
  return(sum(repl[100,] == 0)/x)
}
allele_counts_4(0.5,400,100,1000)
```

```
## [1] 0.007
```

```r
allele_counts_4(0.25,400,100,1000)
```

```
## [1] 0.104
```

```r
allele_counts_4(0.1,400,100,1000)
```

```
## [1] 0.428
```

5.
![](12345_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
## 
## [[6]]
## NULL
## 
## [[7]]
## NULL
## 
## [[8]]
## NULL
## 
## [[9]]
## NULL
## 
## [[10]]
## NULL
## 
## [[11]]
## NULL
## 
## [[12]]
## NULL
## 
## [[13]]
## NULL
## 
## [[14]]
## NULL
## 
## [[15]]
## NULL
## 
## [[16]]
## NULL
## 
## [[17]]
## NULL
## 
## [[18]]
## NULL
## 
## [[19]]
## NULL
## 
## [[20]]
## NULL
## 
## [[21]]
## NULL
## 
## [[22]]
## NULL
## 
## [[23]]
## NULL
## 
## [[24]]
## NULL
## 
## [[25]]
## NULL
## 
## [[26]]
## NULL
## 
## [[27]]
## NULL
## 
## [[28]]
## NULL
## 
## [[29]]
## NULL
## 
## [[30]]
## NULL
## 
## [[31]]
## NULL
## 
## [[32]]
## NULL
## 
## [[33]]
## NULL
## 
## [[34]]
## NULL
## 
## [[35]]
## NULL
## 
## [[36]]
## NULL
## 
## [[37]]
## NULL
## 
## [[38]]
## NULL
## 
## [[39]]
## NULL
## 
## [[40]]
## NULL
## 
## [[41]]
## NULL
## 
## [[42]]
## NULL
## 
## [[43]]
## NULL
## 
## [[44]]
## NULL
## 
## [[45]]
## NULL
## 
## [[46]]
## NULL
## 
## [[47]]
## NULL
## 
## [[48]]
## NULL
## 
## [[49]]
## NULL
## 
## [[50]]
## NULL
## 
## [[51]]
## NULL
## 
## [[52]]
## NULL
## 
## [[53]]
## NULL
## 
## [[54]]
## NULL
## 
## [[55]]
## NULL
## 
## [[56]]
## NULL
## 
## [[57]]
## NULL
## 
## [[58]]
## NULL
## 
## [[59]]
## NULL
## 
## [[60]]
## NULL
## 
## [[61]]
## NULL
## 
## [[62]]
## NULL
## 
## [[63]]
## NULL
## 
## [[64]]
## NULL
## 
## [[65]]
## NULL
## 
## [[66]]
## NULL
## 
## [[67]]
## NULL
## 
## [[68]]
## NULL
## 
## [[69]]
## NULL
## 
## [[70]]
## NULL
## 
## [[71]]
## NULL
## 
## [[72]]
## NULL
## 
## [[73]]
## NULL
## 
## [[74]]
## NULL
## 
## [[75]]
## NULL
## 
## [[76]]
## NULL
## 
## [[77]]
## NULL
## 
## [[78]]
## NULL
## 
## [[79]]
## NULL
## 
## [[80]]
## NULL
## 
## [[81]]
## NULL
## 
## [[82]]
## NULL
## 
## [[83]]
## NULL
## 
## [[84]]
## NULL
## 
## [[85]]
## NULL
## 
## [[86]]
## NULL
## 
## [[87]]
## NULL
## 
## [[88]]
## NULL
## 
## [[89]]
## NULL
## 
## [[90]]
## NULL
## 
## [[91]]
## NULL
## 
## [[92]]
## NULL
## 
## [[93]]
## NULL
## 
## [[94]]
## NULL
## 
## [[95]]
## NULL
## 
## [[96]]
## NULL
## 
## [[97]]
## NULL
## 
## [[98]]
## NULL
## 
## [[99]]
## NULL
## 
## [[100]]
## NULL
```
6(1).

```r
set.seed(720)
pval <- function(intercept, slope, size, sd) {
x <- seq(from =1, to = 10, length.out = size) 
y_deterministic <- intercept + slope*x
y_simulated <- rnorm(length(x), mean = y_deterministic, sd)
mod_sim <- lm(y_simulated ~ x)
p_val_slope <- summary(mod_sim)$coef[2,4]
p_val_slope
}
replicate(n=100,pval(0.5,0.1,20,2))
```

```
##   [1] 0.362062526 0.217974167 0.325034673 0.677325090 0.157603696
##   [6] 0.544150051 0.001007433 0.220102889 0.697565361 0.660627728
##  [11] 0.483512912 0.811256356 0.390050687 0.034661125 0.041481784
##  [16] 0.673833634 0.163332876 0.263695337 0.058104974 0.785159791
##  [21] 0.184864751 0.548399577 0.708047360 0.721201379 0.472971946
##  [26] 0.476600043 0.681806670 0.439674193 0.007749069 0.101583283
##  [31] 0.570335114 0.363456484 0.981378755 0.279368058 0.059317198
##  [36] 0.666428397 0.760305358 0.212748656 0.112076261 0.993768181
##  [41] 0.028560124 0.199968574 0.987352386 0.570840828 0.113020343
##  [46] 0.170249123 0.332478686 0.379061906 0.013045769 0.631443481
##  [51] 0.351408960 0.061309577 0.429809936 0.906217448 0.987918167
##  [56] 0.575060544 0.399641993 0.560329435 0.162860275 0.449256724
##  [61] 0.681368936 0.813694047 0.231167243 0.304744707 0.774859478
##  [66] 0.516522132 0.215979399 0.889648720 0.024917109 0.387601015
##  [71] 0.898891052 0.367253606 0.093460586 0.654086377 0.730831050
##  [76] 0.643812027 0.004181689 0.914436258 0.619802264 0.502260319
##  [81] 0.655293738 0.643036472 0.892473791 0.344221021 0.811120088
##  [86] 0.466932196 0.386107926 0.834578225 0.014378627 0.020349513
##  [91] 0.562775147 0.103054224 0.213018061 0.809534823 0.030115352
##  [96] 0.053753390 0.822394494 0.109312920 0.709089451 0.263872850
```

6(2)

```r
hist(replicate(n=1000,pval(0.5,0.1,20,2)))
```

![](12345_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
table((replicate(n=1000,pval(0.5,0.1,20,2)))<0.05)/length((replicate(n=1000,pval(0.5,0.1,20,2)))<0.05)
```

```
## 
## FALSE  TRUE 
## 0.924 0.076
```

6(3)

The result of (slope = 0) is smaller than it of (slope = 0.1). Because when slope = 0, mean = y_deterministic = a. There is no linear relationship so you won't have a high p value expectation.

```r
hist(replicate(n=1000,pval(0.5,0,20,2)))
```

![](12345_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
table((replicate(n=1000,pval(0.5,0,20,2)))<0.05)/length((replicate(n=1000,pval(0.5,0,20,2)))<0.05)
```

```
## 
## FALSE  TRUE 
## 0.947 0.053
```

6(4)

```r
replicates <- rep(NA,100)
samples <- seq(10, 100, 5)
proportion <- rep(NA,length(samples))
for(i in 1:length(samples)) {
replicates <- replicate(n=100,pval(0.5,0.1,samples[i],1.5))
  proportion[i] <- sum(replicates<0.05)/length(replicates)
} 
plot(x= samples, y = proportion,  xlab="sample size", 
	 ylab="proportion of p-values less than 0.05", main = "sample size VS proportion of p-values less than 0.05")
```

![](12345_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

