##  Random nucleotides for testing:
writeLines(paste(sample(nuc, 120000, replace = T), collapse=""), "/Users/matthew/PycharmProjects/BigBoiCode/nucs.fasta")
## for plotting:
library(tidyverse)
par(mfrow=c(2,2))
data <- read_csv("/Users/matthew/PycharmProjects/BigBoiCode/zipf", col_names = F)
`Frequency Ranking` <-  (1:length(data[[1]]))
`Zipf Enumeration` <- (data[[1]])
plot(`Frequency Ranking`, `Zipf Enumeration`, type="l")
title("Evaluating for Zipf Distribution")

`log(Frequency Ranking)` <-  log(1:length(data[[1]]))
`log(Zipf Enumeration)` <- log(data[[1]])
plot(`log(Frequency Ranking)`, `log(Zipf Enumeration)`, type="l")
title("Evaluating for Zipf Distribution (logarithmic)")

data <- read_csv("/Users/matthew/PycharmProjects/BigBoiCode/zipf.freq.txt", col_names = F)
`Rank` <-  (1:length(data[[1]]))
`Frequency` <- (data[[1]])
plot(`Rank`, `Frequency`, type="l")
title("Evaluating Frequency vs Rank");data <- read_csv("/Users/matthew/PycharmProjects/BigBoiCode/zipf.freq.txt", col_names = F)

`log(Rank)` <-  log(1:length(data[[1]]))
`log(Frequency)` <- log(data[[1]])
plot(`log(Rank)`, `log(Frequency)`, type="l")
title("Evaluating Frequency vs Rank (logarithmic)")
