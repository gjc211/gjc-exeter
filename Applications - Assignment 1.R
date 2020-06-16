# 1. Converting to matrix, then solving
A <- array(c(1,2,-3,1,1,3,4,2,1,0,1,-1,1,-1,2,1), dim=c(4,4))
b <- c(13,-1,10,1)
solve(A,b)
# x=2, y=0, z=6, w=5

# 2.a) Converting to form ready for PCA 
food_ <- read.csv("food.csv", header = TRUE)
head(food_)

# Making sure observations in rows and variables in columns by transposing the data: 
# Making sure headers are not recognised as part of data before transposing
rownames(food_) <- food_[,1]
food_ <- food_[,-1]
food_ = t(food_)

# b) PCA
# Centralising data without changing dimensions and finding covariance matrix
food.centred <- scale(food_,center=TRUE,scale=FALSE)
food.cov <- cov(food.centred)

# Finding eigenvectors and eigenvalues and determining how many principal components needed
food.eig <- eigen(food.cov)
eigenvectors <- food.eig$vectors
eigenvalues <- food.eig$values
eigenvalues
sum(eigenvalues[1:2])/sum(eigenvalues)
# First two principal components contains 96% of information

# c) 
# Transforming dataset onto first principal component and turing into dataframe
food.reduced <- food.centred %*% eigenvectors[,1]
colnames(food.reduced) <- c("PC1")
food.reduced <- as.data.frame(food.reduced)
Country = c("England", "Wales", "Scotland", "N Ireland")

# Finally plotting onto graph with the one principal component 
library(ggplot2)
ggplot(data=food.reduced, aes(x = PC1, y = 0)) + geom_point(aes(color=Country)) + theme(axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank(),legend.position = "bottom", panel.grid.minor = element_blank()) + theme(aspect.ratio = 0.3)

# Now using two principal components
food.reduced1 <- food.centred %*% eigenvectors[,1:2]
colnames(food.reduced1) <- c("PC1", "PC2")
food.reduced1 <- as.data.frame(food.reduced1)
Country = c("England", "Wales", "Scotland", "N Ireland")
ggplot(data=food.reduced1, aes(x = PC1, y = PC2)) + geom_point(aes(color=Country)) + theme(aspect.ratio = 1)

# N Ireland is clear outlier as seen by using just the first principle component
# The addition of the second pricipal component shows the variation of the other three countries 
# along that axis, in terms of eating habits.

# d)
# N Ireland consume much less of a lot of foods such as fish, but for example consume more potatoes.


# 3. Image compression
# a) Separating the matrices formed from the colour image for PCA to be performed on each 
install.packages("jpeg")
library(jpeg)
sheep <- readJPEG("sheep.jpg")
array1 <- sheep[,,1]
array2 <- sheep[,,2]
array3 <- sheep[,,3]
sheep1.pca <- prcomp(array1,center=FALSE)
sheep2.pca <- prcomp(array2,center=FALSE)
sheep3.pca <- prcomp(array3,center=FALSE)

# b)
# Constructing a compressed version of each colour channel, using a chosen number of principal components
sheep.compressed1 <- sheep1.pca$x[,1:200] %*% t(sheep1.pca$rotation[,1:200])
sheep.compressed2 <- sheep2.pca$x[,1:200] %*% t(sheep2.pca$rotation[,1:200])
sheep.compressed3 <- sheep3.pca$x[,1:200] %*% t(sheep3.pca$rotation[,1:200])

# c)
# Compiling the separately compressed matrices back into an array 
finalarray1 <- array(c(sheep.compressed1,sheep.compressed2,sheep.compressed3),dim = c(2248,4000,3)) 
finalarray2 <- array(c(sheep.compressed1,sheep.compressed2,sheep.compressed3),dim = c(2248,4000,3))
finalarray3 <- array(c(sheep.compressed1,sheep.compressed2,sheep.compressed3),dim = c(2248,4000,3))
finalarray4 <- array(c(sheep.compressed1,sheep.compressed2,sheep.compressed3),dim = c(2248,4000,3))
finalarray5 <- array(c(sheep.compressed1,sheep.compressed2,sheep.compressed3),dim = c(2248,4000,3))

# Exporting the final jpg image 
writeJPEG(finalarray1, paste("sheep_compressed_20.jpg"))
writeJPEG(finalarray2, paste("sheep_compressed_5.jpg"))
writeJPEG(finalarray3, paste("sheep_compressed_2.jpg"))
writeJPEG(finalarray4, paste("sheep_compressed_60.jpg"))
writeJPEG(finalarray5, paste("sheep_compressed_200.jpg"))

# Finding compression ratios
file.info("sheep_compressed_60.jpg")$size/file.info("sheep.jpg")$size*100
file.info("sheep_compressed_2.jpg")$size/file.info("sheep.jpg")$size*100
file.info("sheep_compressed_5.jpg")$size/file.info("sheep.jpg")$size*100
file.info("sheep_compressed_20.jpg")$size/file.info("sheep.jpg")$size*100
file.info("sheep_compressed_200.jpg")$size/file.info("sheep.jpg")$size*100

# Compression ratio %:
# 20 Principal components = 7.68
# 5 Principal components = 5.96
# 2 Principal components = 5.32
# 60 Principal components = 9.80
# 200 Principal components = 12.00






