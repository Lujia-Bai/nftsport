# nftsport
Bootstrap-based portmanteau tests for multivariate  non-stationary functional time series with an increasing number of lags with a data example in energy. The paper is in https://arxiv.org/abs/2501.00118. 

## Quick start 

```{r}
# install.packages("devtools")
devtools::install_github("Lujia-Bai/nftsport")
```

## Data 
We analyze a three-dimensional time series of hourly energy consumption, measured in megawatts (MW), for American Electric Power, The Dayton Power and Light Company, and Northern Illinois Hub. Each day’s consumption profile is treated as a continuous function over 24 hours.
Given that 2006 ranked as the sixth warmest year globally, we focus on the period from February 16, 2006, to July 4, 2007. This interval comprises 501 days with complete 24-hour observations.

```{r}
data("pjmfts")
```

## Visuallization

You can also embed plots, for example:

```{r pressure, echo=FALSE}
normalize <- function(x) {
  t(apply(x, 1, function(row) (row - mean(row)) / sd(row)))
}
fts = pjmfts
c1 <- normalize(fts[seq(1,500,length.out=6),,1])
c2 <- normalize(fts[seq(1,500,length.out=6),,2])
c3 <- normalize(fts[seq(1,500,length.out=6),,3])

# 5 panels (one per day)
par(mfrow = c(2, 3))

for (i in 1:nrow(c1)) {
  matplot(t(c1[i,,drop=FALSE]),
          type = "l",
          col = "blue",
          lty = 1,
          ylim = range(c(c1[i,], c2[i,], c3[i,])),
          main =   unique_dates[seq(1,500,length.out= 6)+499][i],
          xlab = "Time",
          ylab = "Value")
  
  matlines(t(c2[i,,drop=FALSE]), col = "red", lty = 1)
  matlines(t(c3[i,,drop=FALSE]), col = "green", lty = 1)
  
  legend("bottomright",
         legend = c("Company 1","Company 2","Company 3"),
         col = c("blue","red","green"),
         lty = 1,
         cex = 0.8)
}
```

## Functional portmanteau test

concurrent functional observations 

```{r}
# original test
result = fport(fts, 3, lb=1, ub=3)
print(result$pvalue)
```

including non-concurrent functional observations  
```{r}
# extended test
resulte = fport_extend(fts, 3, lb=1, ub=3)
print(resulte$pvalue)
```
