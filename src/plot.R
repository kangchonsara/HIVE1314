data = read.table("HIVE1314_similarities.csv", sep=",", header=TRUE)
data = na.omit(data)
data
plot (data$yob, data$Asimilarity)
plot (data$yob, data$Bsimilarity)
plot (data$yob, data$HA2similarity)

Asim.lm = lm(data$yob ~ data$Asimilarity)
Asim.res = resid(Asim.lm)

Asim.res

plot(data$post, Asim.res)
Asim.fit = lm(data$post ~ Asim.res)
summary(Asim.fit)