# Simulate data:

status <- c(rep(1,500),rep(2,500))
score.null <- rnorm(1000,mean=50,sd=25) # null data; no difference in group means
null.data <- cbind(status,score.null)

score.altern <- c(rnorm(500,mean=50,sd=5),rnorm(500,mean=51,sd=5)) # alternative data; one group has mean of 50, other has 51
altern.data <- cbind(status,score.altern)

null.betas <- vector(length=1000)
null.p <- vector(length=1000)

a.betas <- vector(length=1000)
a.p <- vector(length=1000)

# Examine difference in null vs alternative:

null.model <- lm(score.null~as.factor(status)) # no sig difference in true data
a.model <- lm(score.altern~as.factor(status)) # sig difference in true data

# Now bootstrap:

for(i in 1:1000){ 
  sample.index <- sample(1:1000, 100,replace=TRUE)
  # resample from null data 1k times and save beta and b-values from models
  null.perm <- lm(score.null[sample.index]~as.factor(status[sample.index]))
  null.betas[i] <- as.numeric(null.perm$coefficients[2])
  null.p[i]<- anova(null.perm)$'Pr(>F)'[1]
  # resample from alternative data 1k times and save beta and b-values from models
  a.perm <- lm(score.altern[sample.index]~as.factor(status[sample.index]))
  a.betas[i] <- as.numeric(a.perm$coefficients[2])
  a.p[i]<- anova(a.perm)$'Pr(>F)'[1]
}

plot(a.p,a.betas) # when a true effect exists in the data, greater number positive slopes than negative, more so at low p-values

plot(null.p,null.betas) # when no true effect, distribution of sign of slop is distributed equally (approx)
