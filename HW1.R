#get car data
CarData <- read.table("carc.dat")
#extract mileage and headqtr from table
x <- CarData[c(3,14)]

#PART 1: fivenum summaries of US, Japan, EU
fivenum(x[x$V14 == 1, 1])
fivenum(x[x$V14 == 2, 1])
fivenum(x[x$V14 == 3, 1])

#PART 2: factor headqtrs and draw boxplot
headquarter <- factor(x[,2])
levels(headquarter) <- c("US","JAPAN", "EU")
x <- cbind(x, headquarter)
colnames(x)=c("Mileage", "HeadInt", "HeadFac")
boxplot(Mileage ~ HeadFac, data=x, at=1:3, main="Car Mileage Data", col='green')

#PART 3: draw histogram for total mileage & country groups
hist(x[,1], main="All Mileage Data", xlab="Mileage")

par(mfrow = c(3, 1))
hist(x[headquarter=="US",1], main="US", xlab="Mileage")
hist(x[headquarter=="JAPAN",1], main="JAPAN", xlab="Mileage")
hist(x[headquarter=="EU",1], main="EU", xlab="Mileage")







