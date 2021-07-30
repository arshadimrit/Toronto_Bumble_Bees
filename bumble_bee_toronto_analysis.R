library(dplyr)
library(Hmisc)
library(vegan)

library(ggfortify)
library(factoextra)
library(corrplot)
library(rpart)
library(ggrepel)

library(ade4)

###########################################################################################################
# Making correlations plot

df = read.csv('/filepath/toronto_bumble_bees_grid.csv')

df$aveMeanFD = log1p(df$aveMeanFD)  # log transform mean foraging distance

data_subset = df[ , c('aveMeanFD', 'colNe', 'buildPerc', 'roofPerc', 'meadPerc', 'treeCount', 
                    'water_perc', 
                    'elevat', 'houseDensity', 'succPerc', 'green1Perc', 
                    'tree_canopy_perc', 'beachPerc', 'wetPerc', 'grass_shrub_perc', 
                    'road_perc', 'slope', 'forestPerc', 'bare_earth_perc', 
                    'other_paved_perc', 
                    'popDensity', 'popFemale', 'popMale', 'popTotalCom', 'pop_60plus', 
                    'pop_less20', 'pop_20.39', 'pop_40.59', 'indTI', 'famTI')]

corr_data = subset(data_subset) # making a copy 

corr_data = corr_data %>%
  rename(
    "Building %" = buildPerc, "Green Roof %" = roofPerc, "Meadow %" = meadPerc, 
    "Tree Count" = treeCount, "Water %" = water_perc, "Average Elevation" = elevat, 
    "House Density" = houseDensity, "Successional %" = succPerc, "City Park %" = green1Perc, 
    "Tree Canopy %" = tree_canopy_perc, "Beach-Bluff %" = beachPerc, "Wetland %" = wetPerc, 
    "Grass & Shrub %" = grass_shrub_perc, "Road %" = road_perc, "Average Slope" = slope, 
    "Forest %" = forestPerc, "Bare Earth %" = bare_earth_perc, "Other Paved %" = other_paved_perc, 
    "log (Foraging Distance)" = aveMeanFD, "Number Colonies" = colNe, "Human Population Density" = popDensity, 
    "Number Human Females" = popFemale, "Number Human Males" = popMale, "Total Human Population" = popTotalCom, 
    "Population > 60 yrs" = pop_60plus, " Population < 20 yrs" = pop_less20, 
    "Population 20-39 yrs" = pop_20.39, "Population 40-59" = pop_40.59, "Total Individual Income" = indTI, 
    "Total Family Income" = famTI
  )

res = rcorr(as.matrix(corr_data), type = c("spearman")) # run spearman correlations
p_benjamini = matrix(p.adjust(res$P, 'BH'), nrow=30, ncol=30) # adjust p-values using Benjamini-Hochberg

corrplot(res$r, type="upper", order="original", 
         p.mat = p_benjamini, sig.level = 0.05, insig = "blank", tl.col="black")

############################### RDA ################################


df_pred_response = subset(data_subset)

df_pred_response = df_pred_response %>%
  rename(
    "Buildings" = buildPerc, "Meadows" = meadPerc, 
    "Tree_Count" = treeCount, "Elevation" = elevat, 
    "House_Density" = houseDensity, "City_Parks" = green1Perc, 
    "Tree_Canopy" = tree_canopy_perc, 
    "Grass_Shrubs" = grass_shrub_perc, "Roads" = road_perc,
    "Forests" = forestPerc, "Bare_Earth" = bare_earth_perc, "Other_Paved_Surfaces" = other_paved_perc, 
    "Foraging_Distance" = aveMeanFD, "Number_Colonies" = colNe, "Human_pop_Density" = popDensity, 
    "Pop_over_60" = pop_60plus, "Pop_40_59" = pop_40.59, "Total_Individual_Income" = indTI, 
    "Human_total_Pop" = popTotalCom)

vif_pruned_rda = rda(df_pred_response[, c(1,2)] ~ 
                       Buildings + Meadows + Forests + House_Density + City_Parks +
                       Roads + Bare_Earth + Other_Paved_Surfaces + Human_pop_Density + 
                       Total_Individual_Income + Tree_Count + Elevation + Tree_Canopy + 
                       Human_total_Pop + Pop_over_60 + + Pop_40_59 + Grass_Shrubs, 
                     data = df_pred_response)

RsquareAdj(vif_pruned_rda)
summary(vif_pruned_rda)

anova_global = anova(vif_pruned_rda, step=1000) # testing the significance of the RDA model
anova_axis = anova.cca(vif_pruned_rda, by='axis', step=1000)

anova_global
anova_axis

# Trying to make beautiful rda plot

rda_summary = summary(vif_pruned_rda)
sp = as.data.frame(rda_summary$species[, 1:2])
st = as.data.frame(rda_summary$sites[, 1:2])
yz = as.data.frame(rda_summary$biplot[, 1:2])

ggplot() +
  geom_segment(data = sp,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "red")+
  geom_text_repel(data = sp,aes(RDA1,RDA2,label=row.names(sp)))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=0.6,colour = "blue")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)), box.padding = 0.8)+
  labs(x=paste("RDA 1 (60.9 %)", sep=""),
       y=paste("RDA 2 (30.2 %)", sep=""))+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank())

### Variation partitioning ###

varpart(Y = df_pred_response[, c(1,2)], X = ~ Buildings + Other_Paved_Surfaces + House_Density +
          Roads, 
        ~ Meadows + ~ Forests  + City_Parks  + Bare_Earth + Tree_Count + Elevation + Tree_Canopy + 
          Grass_Shrubs, 
        ~ Human_pop_Density + Total_Individual_Income + Human_total_Pop + Pop_over_60 + + Pop_40_59, 
        data = df_pred_response)
