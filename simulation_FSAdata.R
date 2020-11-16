library(FSA)
library(FSAdata)
library(magrittr)



data(CreekChub)

CreekChub %<>% mutate(lcat10 = lencat(len, w=20))

sample2 <- filter(CreekChub, !is.na(age))
alk_freq <- xtabs(~lcat10 + age, data = sample2)
alk <- prop.table(alk_freq, margin = 1)
length_n <- xtabs(~lcat10, data = CreekChub)
alkAgeDist(alk, lenA.n = rowSums(alk_freq), len.n = length_n)



