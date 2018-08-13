library(dplyr)
library(ggplot2)
library(Rtsne)
library(plotly)

#load data
df <- read.csv("outdata11.csv")
df3 <- df %>% select(-X,-SubjectID,-ALSFRS_slope,-onset_delta.y)

#Kmeans
set.seed(2017)
# K-Means Cluster Analysis 1
fit0 <- kmeans(df3, 4) # 4 cluster solution
# append cluster assignment 1
df4 <- data.frame(df3, fit0$cluster)

# 2D tsne
tsne_out_2d <- as.data.frame(Rtsne(as.matrix(df4[,-1]),dims=2)$Y)
tsne_out_2d$kmeans_cluster_label <- as.factor(df4$fit0.cluster)

ggplot(tsne_out_2d, aes(x=V1, y=V2,color=kmeans_cluster_label)) + 
    geom_point(size=2,alpha=0.8)+ggtitle("2D-tSNE")

# 3D tsne
tsne_out_3d <- as.data.frame(Rtsne(as.matrix(df4[,-1]),dims=3)$Y)
tsne_out_3d$kmeans_cluster_label <- as.factor(df4$fit0.cluster)
p <- plot_ly(tsne_out_3d, x = ~V1, y = ~V2, z = ~V3, color = ~kmeans_cluster_label, marker=list(size=2)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'V1'),
                        yaxis = list(title = 'V2'),
                        zaxis = list(title = 'V3')))
p
