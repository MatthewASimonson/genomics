###############################
# great plotting examples 1:  #
###############################


options(PlosApiKey = "<insert your API key here!>")
#install_github("rplos", "ropensci")
library("rplos")
library("ggplot2")
require("dplyr")

# Convert author strings to counts
countAuths <- function(cell)
  length(unlist(strsplit(cell, ";")))

countAuths <- Vectorize(countAuths)

# Query PLoS API for 1k papers per journal per year,
# count the number of authors and return a data.frame
getAuths <- function(j, lim=1000, start.year=2006){
  cat("Getting results for journal: ", j, "\n")
  # seem to be in reverse order by year?
  results <- sapply(start.year:2013, function(i) data.frame(year = i,
                                                            auths = searchplos(
                                                              q  = paste0('publication_date:[', i,
                                                                          '-01-01T00:00:00Z TO ', i,
                                                                          '-12-31T23:59:59Z]'),
                                                              fl = "author",
                                                              fq = list("doc_type:full",
                                                                        paste0("cross_published_journal_key:", j)),
                                                              start=0, limit=lim, sleep=6),
                                                            year=i), simplify=F)
  results <- do.call(rbind, results)
  results$counts <- countAuths(results$author)
  results$journal <- j
  results
}

journals <- journalnamekey()
plos.all <- sapply(journals[c(1:5, 7)], getAuths, simplify=F)
plos <- do.call(rbind, plos.all)

# Fig. 1: Bean plot showing distribution of author counts
#         per journal overall
ggplot(plos, aes(x=journal, y=counts, fill=journal)) +
  geom_violin(scale="width") +
  geom_boxplot(width=.12, fill=I("black"), notch=T,
               outlier.size=NA, col="grey40") +
  stat_summary(fun.y="median", geom="point", shape=20, col="white") +
  scale_y_log10(breaks=c(1:5, seq(10, 50, by=10), 100, 200, 300)) +
  coord_flip() + labs(x="", y="Number of authors per paper") +
  theme_classic() + theme(legend.position="none") +
  scale_fill_brewer()

# Fig 2. ECDFs of the author count distributions
# 5in x 5in
ggplot(plos, aes(x=counts, col=journal)) +
  stat_ecdf(geom="smooth", se=F, size=1.2) + theme_bw() +
  scale_x_log10(breaks=c(1:5, seq(10, 50, by=10), 100, 200, 300)) +
  theme(legend.position=c(.75,.33)) +
  labs(x="Number of authors per paper", y="ECDF",
       col="") + coord_cartesian(xlim=c(1,300)) +
  scale_color_brewer(type="qual", palette=6)

# Fig 3. Trends in author counts over time with
#        confidence limits on the means
# 7 x 7
ggplot(plos, aes(x=year, y=counts, col=journal, fill=journal)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               width=.2, alpha=I(.5)) +
  stat_summary(fun.y="mean", geom="line") +
  labs(list(x="Year", y="Mean number of authors per paper")) +
  theme_bw() + theme(legend.position=c(.2,.85)) +
  scale_fill_brewer(type="qual", palette=2,
                    guide=guide_legend(direction="vertical",
                                       label.position="bottom",
                                       title=NULL, ncol=2,
                                       label.hjust=0.5)) +
  scale_color_brewer(type="qual", palette=2, guide="none")


# from http://stackoverflow.com/a/17024184/1274516
# show regression equation on each graph facet
lm_eqn  <-  function(df){
  m  <- summary(lm(counts ~ year, df))
  eq <- substitute(~~y~"="~beta*x+i~(R^2==r2),
                   list(beta = format(m$coefficients[2,"Estimate"],
                                      digits = 3),
                        i = format(m$coefficients[1,"Estimate"], digits=3),
                        r2 = format(m$r.squared, digits=2)))
  as.character(as.expression(eq))
}

means <- group_by(plos, journal, year) %.% summarise(counts=mean(counts))
b <- by(means, means$journal, lm_eqn)
df <- data.frame(beta=unclass(b), journal=names(b))
summary(lm(counts ~ year + journal, data=means))

means <- group_by(means, journal) %.% summarise(m=max(counts))
df$top <- means$m * 1.2

# Fig 4. Facetted linear regression of author inflation per journal
# 6 x 8.5
ggplot(plos, aes(x=year, y=counts, col=journal, fill=journal)) +
  stat_summary(fun.data="mean_cl_boot", geom="errorbar",
               width=.2, alpha=I(.5)) +
  stat_summary(fun.y="mean", geom="point") +
  #stat_summary(fun.y="median", geom="point", shape=4) +
  facet_wrap(~journal, scales="free_y") +
  geom_smooth(method="lm") +
  scale_x_continuous(breaks=2006:2013) +
  labs(list(x="", y="Mean number of authors per paper")) +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_brewer(type="qual", palette=2, guide="none") +
  scale_color_brewer(type="qual", palette=2, guide="none") +
  geom_text(data=df, aes(x=2009.5, y=top, label=beta), size=3, parse=T)

# Overall estimate of author inflation?
# .21 extra authors per paper per year, on average
s <- summary(lm(counts ~ year + journal, data=plos))


# Summary barchart data:
bc <- data.frame(journal = unique(means$journal),
                 trend   = c(0.2490979,
                             0.1211823,
                             0.5201688,
                             0.4088536,
                             0.05894102,
                             0.1828939),
                 std.err = c(0.08224567,
                             0.02213142,
                             0.1493662,
                             0.06361849,
                             0.03891493,
                             0.03798822),
                 IF      = c(12.690,
                             4.867,
                             8.517,
                             15.253,
                             3.730,
                             8.136))

bc$journal <- factor(bc$journal, levels=bc$journal[order(bc$trend)])

# Fig 5. Barchart of author inflation estimate per journal.
# 7 x 5
ggplot(bc, aes(x=journal, y=trend, fill=journal, ymin=trend-std.err,
               ymax=trend+std.err)) +
  geom_bar(stat="identity") +
  geom_errorbar(width=.2) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  labs(x="",
       y="Estimate of annual author inflation (additional mean authors per paper)") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_brewer(palette="Blues", guide="none")

pcc <- cor(bc$trend, bc$IF)
# Fig 6. Correlation of author inflation and journal impact factors.
# 5 x 5
ggplot(bc, aes(x=trend, y=IF, col=journal)) +
  geom_text(aes(label=journal)) + xlim(0,.6) +
  labs(x="Author inflation estimate",
       y="Journal impact factor (2012)") +
  scale_color_brewer(type="qual", palette=2, guide="none") +
  annotate("text", x=.05, y=15,
           label=paste0("rho == ", format(pcc, digits=2)), parse=T)

# N.S. (p = 0.18)
cor.test(bc$trend, bc$IF)



###############################
# great plotting examples 2:  #
###############################

library("ggplot2")
library("scales")

## Data: http://www.theguardian.com/news/datablog/2009/oct/21/icm-poll-data-labour-conservatives#data (State of the parties spreadsheet)
# ICM poll resultsof ~1000 people taken every month or so
# (more often before an election) as well as election
# results dating from June 1984.
sop <- read.csv("StateOfTheParties.csv", stringsAsFactors=F)

## Data cleanup
sop[,2:5] <- apply(sop[,2:5], 2, function(x) as.numeric(gsub("%", "", x)))
sop[,1] <- as.Date(sop[,1], format="%d-%m-%Y")
colnames(sop)[1] <- "Date"

# correct for some rounding errors leading to 101/99 %
sop$rsum <- apply(sop[,2:5], 1, sum)
table(sop$rsum)
sop[,2:5] <- sop[,2:5] / sop$rsum

# Melt for ggplot
melt.fun <- function(x)
  melt(x, measure.vars=c("CON", "LAB", "LIB.DEM", "OTHER"), id.vars="Date")

props <- melt.fun(sop)

# Grab election data for additional layer
elections <- sop[which(grepl("ELECTION", sop$Sample)),]
elections <- melt.fun(elections)

# Party colours (ish)
cols <- brewer.pal(8, "Pastel1")[c(2,1,6,4)]

# Area plot overview of public support for main parties over 30 years
ggplot(props, aes(x=Date, y=value, fill=variable, group=variable)) +
  geom_area() + scale_fill_manual(values=cols) +
  geom_bar(data=elections, stat="identity", position=position_stack(.1),
           width=I(100), col="grey40") +
  scale_y_continuous(expand=c(0,0), labels=percent_format()) +
  scale_x_date(expand=c(0,0)) +
  labs(list(y="Public support", fill="Party", x="")) #+
ggtitle("UK general elections 1984-2014")


## Look at the time between elections
e.dates <- c(unique(elections$Date), as.Date("2015-05-07"))
results <- c(rep("CON",2), rep("LAB", 3), "CON.LIB", "?")

sop.rows <- which(grepl("ELECTION", sop$Sample))
e.years <- c(sub(".*(\\d{4})$", "\\1", sop[sop.rows, "Sample"], perl=T), "2015")

sop$runup <- rep(c(1:7), c(sop.rows[1], diff(c(sop.rows, nrow(sop)))))
sop$win <- rep(results, c(sop.rows[1], diff(c(sop.rows, nrow(sop)))))
e <- sop

e <- melt(e, measure.vars=c("CON", "LAB", "LIB.DEM", "OTHER"),
          id.vars=c("runup", "Date", "win"))

e$prior <- with(e, as.numeric(Date - e.dates[runup]))
e$header <- with(e, paste0(e.years[e$runup], " (", win, ")"))

# Fn to make colours legible on white bg
darken <- function(hex.col, amount=.15)
  return(rgb(t(col2rgb(hex.col)*(1-amount)), max=255))
darken <- Vectorize(darken)
dcols <- unname(darken(cols))

# Previous 6 election support history
ggplot(subset(e, header != "2015 (?)"),
       aes(x=prior, y=value, col=variable, fill=variable)) +
  facet_wrap(~header, scales="free_x") +
  scale_color_manual(values=dcols) +
  scale_fill_manual(values=dcols) +
  geom_smooth(method="loess") +
  geom_line() + geom_point() + theme_bw() +
  scale_y_continuous(labels=percent_format()) +
  labs(list(y="Public support", x="Days relative to election",
            col="Party", fill="Party"))

# boxlots of election runups
ggplot(e, aes(x=header, y=value, fill=variable)) +
  theme_bw() + geom_boxplot(position=position_dodge(.8)) +
  scale_fill_manual(values=cols) +
  labs(list(y="Public support", x="Days relative to election",
            col="Party", fill="Party")) +
  ggtitle("UK general elections 1984-2014") +
  theme(legend.position=c(.9,.85), legend.background=element_blank())

## Which prior results are most predictive of end result

ggplot(subset(e, variable %in% c("LAB", "CON")),
       aes(x=prior, y=value, col=variable, fill=variable)) +
  facet_grid(header~., scales="free_y") +
  scale_color_manual(values=dcols) +
  scale_fill_manual(values=dcols) +
  geom_smooth(method="loess") +
  scale_x_continuous(expand=c(0.01,0.01)) +
  theme_bw() + scale_y_continuous(labels=percent_format()) +
  labs(list(y="Public support", x="Days relative to election",
            col="Party", fill="Party"))
