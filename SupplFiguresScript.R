library(ggplot2)
library(viridis)
library(stringr)
library(gridExtra)
library(dplyr)
library(ape)
library(cowplot)
library(lme4)
library(lmerTest)
library(data.table)
library(ggrepel)


colourlist <- c('#543005',
                '#8c510a',
                '#bf812d',
                '#dfc27d',
                '#c7eae5',
                '#80cdc1',
                '#01665e')


source('/Applications/polyDFE-master/postprocessing.R')


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}




################################################ Suppl. Fig. 1: Species DFEs, all models ################################################################

species_list = c('Quercus_petraea','Betula_pendula','Populus_nigra', 'Pinus_pinaster','Fagus_sylvatica', 'Picea_abies')
model_average_statistics = data.frame()
discretised_plot_list = data.frame()
discretised_plot_list_model_avg = data.frame()
samplenumber = 40
num = 1

for (species in species_list){
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  print(sp_short)
  print(num)
  speciesstr = str_replace(species, '_', ' ')
  plotstr = paste(species, samplenumber, sep = '.')
  print(plotstr)
  
  
  path = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/'
  path = paste(path, sp_short, sep = '')
  print(path)
  
  est_mod1 = paste(path, 'PolyDFE_4model1_all40.out', sep = '/')
  est_mod2 = paste(path, 'PolyDFE_4model2_all40.out', sep = '/')
  est_mod3 = paste(path, 'PolyDFE_4model3_all40.out', sep = '/')
  est_mod4 = paste(path, 'PolyDFE_4model4_all40.out', sep = '/')
  
  
  
  est1 <- parseOutput(est_mod1)
  dfe1 <- data.frame(getDiscretizedDFE(est1[[1]], c(-100, -10, -1, 0, 1)))
  est2 <- parseOutput(est_mod2)
  dfe2 <- data.frame(getDiscretizedDFE(est2[[1]], c(-100, -10, -1, 0, 1)))
  est3 <- parseOutput(est_mod3)
  dfe3 <- data.frame(getDiscretizedDFE(est3[[1]], c(-100, -10, -1, 0, 1)))
  est4 <- parseOutput(est_mod4)
  dfe4 <- data.frame(getDiscretizedDFE(est4[[1]], c(-100, -10, -1, 0, 1)))
  
  #  if (species %in% c('Picea_abies', 'Pinus_sylvestris')) {
  #    overall <- c(parseOutput(est_mod3), parseOutput(est_mod4))
  #  }
  #  else {
  overall <- c(parseOutput(est_mod1), parseOutput(est_mod2), parseOutput(est_mod3), parseOutput(est_mod4))
  #  }
  
  aic = getAICweights(overall)
  disc_dfe <- sapply(1:length(overall), function(i) getDiscretizedDFE(overall[[i]], c(-100, -10, -1, 0, 1)))
  disc_dfe_model_avg = sum(sapply(1:length(overall), function(i) aic[i, "weight"] * disc_dfe[1,][i]))
  
  disc_dfe_model_avg = apply(disc_dfe, 1, function(x) 
    sum(sapply(1:length(overall), function(i) aic[i, "weight"] * x[i])) )
  
  disc_dfe_model_avg <- cbind(disc_dfe_model_avg, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  disc_dfe_model_avg <- cbind(disc_dfe_model_avg, species)
  colnames(disc_dfe_model_avg) <- c('fraction', 'categories', 'Species')
  discretised_plot_list_model_avg = rbind(discretised_plot_list_model_avg, disc_dfe_model_avg)
  
  
  eps_an <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["eps_an"]))
  theta_bar <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["theta_bar"]))
  Sd <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_d"]))
  b <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["b"]))
  pb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["p_b"]))
  Sb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_b"]))
  alpha_dfe <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * as.numeric(overall[[i]]$alpha[1])))
  
  (print(b))
  prime_values <- data.frame(species, plotstr, t(disc_dfe_model_avg[,1]), b, Sd, pb, Sb, alpha_dfe, eps_an)
  colnames(prime_values) <- c('Species', 'name', "<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1", 'b', 'S_d', 'p_b', 'S_b', 'alpha', 'eps_an')
  model_average_statistics <- rbind(prime_values, model_average_statistics)
  
  
  
  dfe1 <- data.frame(t(dfe1))
  dfe1 <- cbind(dfe1, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  dfe1 <- cbind(dfe1, c(species), c('mod1'))
  
  dfe2 <- data.frame(t(dfe2))
  dfe2 <- cbind(dfe2, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  dfe2 <- cbind(dfe2, c(species), c('mod2'))
  
  dfe3 <- data.frame(t(dfe3))
  dfe3 <- cbind(dfe3, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  dfe3 <- cbind(dfe3, c(species), c('mod3'))
  
  dfe4 <- data.frame(t(dfe4))
  dfe4 <- cbind(dfe4, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  dfe4 <- cbind(dfe4, c(species), c('mod4'))
  
  colnames(dfe1) <- c('fraction', 'categories', 'Species', 'Model')
  dfe1$categories <- factor(dfe1$categories, levels = dfe1$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe1)
  
  colnames(dfe2) <- c('fraction', 'categories', 'Species', 'Model')
  dfe2$categories <- factor(dfe2$categories, levels = dfe2$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe2)
  
  colnames(dfe3) <- c('fraction', 'categories', 'Species', 'Model')
  dfe3$categories <- factor(dfe3$categories, levels = dfe3$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe3)
  
  colnames(dfe4) <- c('fraction', 'categories', 'Species', 'Model')
  dfe4$categories <- factor(dfe2$categories, levels = dfe4$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe4)
  
  num = num + 1
}



discretised_plot_list_model_avg$fraction <- as.numeric(discretised_plot_list_model_avg$fraction)
discretised_plot_list_model_avg$categories <- factor(discretised_plot_list_model_avg$categories, levels = unique(discretised_plot_list_model_avg$categories))

discretised_plot_list$Species <- sub("_", " ", discretised_plot_list$Species)
discretised_plot_list_model_avg$Species <- sub("_", " ", discretised_plot_list_model_avg$Species)
discretised_plot_list$Species <- factor(discretised_plot_list$Species, levels = c('Fagus sylvatica','Quercus petraea','Betula pendula','Populus nigra','Picea abies','Pinus pinaster'))
discretised_plot_list_model_avg$Species <- factor(discretised_plot_list_model_avg$Species, levels = c('Fagus sylvatica','Quercus petraea','Betula pendula','Populus nigra','Picea abies','Pinus pinaster'))


ggplot(discretised_plot_list, aes(fill = Species, colour = Model, x=categories, y=fraction)) +
  scale_fill_manual(values = colourlist) +
  #  scale_color_manual(values = c('black', 'gold')) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() + xlab(expression(paste("Fitness effects (", italic("N"['e']), italic("s"), ")" ) )) + ylab('Fraction of new mutations') + 
  theme(legend.text = element_text(face = c("italic")))+
  coord_cartesian(ylim = c(0, 0.8))


################################################ Suppl. Fig. 2 and 3: Species DFEs, all sites and GC-conservative only ################################################################

species_list = c('Quercus_petraea','Betula_pendula','Populus_nigra', 'Pinus_pinaster','Fagus_sylvatica', 'Picea_abies')

model_average_statistics = data.frame()
discretised_plot_list = data.frame()
discretised_plot_list_model_avg_GC = data.frame()
samplenumber = 40
num = 1

for (species in species_list){
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  print(sp_short)
  print(num)
  speciesstr = str_replace(species, '_', ' ')
  plotstr = paste(species, samplenumber, sep = '.')
  print(plotstr)
  
  
  path = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/'
  path = paste(path, sp_short, sep = '')
  print(path)
  
  est_mod1 = paste(path, 'PolyDFE_4model1_all40GC.out', sep = '/')
  est_mod2 = paste(path, 'PolyDFE_4model2_all40GC.out', sep = '/')
  est_mod3 = paste(path, 'PolyDFE_4model3_all40GC.out', sep = '/')
  est_mod4 = paste(path, 'PolyDFE_4model4_all40GC.out', sep = '/')
  
  
  
  est1 <- parseOutput(est_mod1)
  dfe1 <- data.frame(getDiscretizedDFE(est1[[1]], c(-100, -10, -1, 0, 1)))
  est2 <- parseOutput(est_mod2)
  dfe2 <- data.frame(getDiscretizedDFE(est2[[1]], c(-100, -10, -1, 0, 1)))
  est3 <- parseOutput(est_mod3)
  dfe3 <- data.frame(getDiscretizedDFE(est3[[1]], c(-100, -10, -1, 0, 1)))
  est4 <- parseOutput(est_mod4)
  dfe4 <- data.frame(getDiscretizedDFE(est4[[1]], c(-100, -10, -1, 0, 1)))
  
 
  overall <- c(parseOutput(est_mod1), parseOutput(est_mod2), parseOutput(est_mod3), parseOutput(est_mod4))
  
  aic = getAICweights(overall)
  disc_dfe <- sapply(1:length(overall), function(i) getDiscretizedDFE(overall[[i]], c(-100, -10, -1, 0, 1)))
  disc_dfe_model_avg = sum(sapply(1:length(overall), function(i) aic[i, "weight"] * disc_dfe[1,][i]))
  
  disc_dfe_model_avg = apply(disc_dfe, 1, function(x) 
    sum(sapply(1:length(overall), function(i) aic[i, "weight"] * x[i])) )
  
  disc_dfe_model_avg <- cbind(disc_dfe_model_avg, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  disc_dfe_model_avg <- cbind(disc_dfe_model_avg, species)
  colnames(disc_dfe_model_avg) <- c('fraction', 'categories', 'Species')
  discretised_plot_list_model_avg_GC = rbind(discretised_plot_list_model_avg_GC, disc_dfe_model_avg)
  
  
  eps_an <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["eps_an"]))
  theta_bar <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["theta_bar"]))
  Sd <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_d"]))
  b <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["b"]))
  pb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["p_b"]))
  Sb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_b"]))
  alpha_dfe <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * as.numeric(overall[[i]]$alpha[1])))
  
  (print(b))
  prime_values <- data.frame(species, plotstr, t(disc_dfe_model_avg[,1]), b, Sd, pb, Sb, alpha_dfe, eps_an)
  colnames(prime_values) <- c('Species', 'name', "<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1", 'b', 'S_d', 'p_b', 'S_b', 'alpha', 'eps_an')
  model_average_statistics <- rbind(prime_values, model_average_statistics)
  
  
  
  dfe1 <- data.frame(t(dfe1))
  dfe1 <- cbind(dfe1, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  dfe1 <- cbind(dfe1, c(species), c('mod1'))
  
  dfe2 <- data.frame(t(dfe2))
  dfe2 <- cbind(dfe2, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  dfe2 <- cbind(dfe2, c(species), c('mod2'))
  
  dfe3 <- data.frame(t(dfe3))
  dfe3 <- cbind(dfe3, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  dfe3 <- cbind(dfe3, c(species), c('mod3'))
  
  dfe4 <- data.frame(t(dfe4))
  dfe4 <- cbind(dfe4, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  dfe4 <- cbind(dfe4, c(species), c('mod4'))
  
  colnames(dfe1) <- c('fraction', 'categories', 'Species', 'Model')
  dfe1$categories <- factor(dfe1$categories, levels = dfe1$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe1)
  
  colnames(dfe2) <- c('fraction', 'categories', 'Species', 'Model')
  dfe2$categories <- factor(dfe2$categories, levels = dfe2$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe2)
  
  colnames(dfe3) <- c('fraction', 'categories', 'Species', 'Model')
  dfe3$categories <- factor(dfe3$categories, levels = dfe3$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe3)
  
  colnames(dfe4) <- c('fraction', 'categories', 'Species', 'Model')
  dfe4$categories <- factor(dfe2$categories, levels = dfe4$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe4)
  
  num = num + 1
}

discretised_plot_list_model_avg_GC$fraction <- as.numeric(discretised_plot_list_model_avg_GC$fraction)
discretised_plot_list_model_avg_GC$categories <- factor(discretised_plot_list_model_avg_GC$categories, levels = unique(discretised_plot_list_model_avg_GC$categories))

discretised_plot_list_model_avg_GC$Mutations <- c('GC-conservative')
discretised_plot_list_model_avg$Mutations <- c('All')
discretised_plot_list_model_avg_GC$Species <- sub("_", " ", discretised_plot_list_model_avg_GC$Species)
discretised_plot_list$Species <- sub("_", " ", discretised_plot_list$Species)
discretised_plot_list$Species <- factor(discretised_plot_list$Species, levels = c('Fagus sylvatica','Quercus petraea','Betula pendula','Populus nigra','Picea abies','Pinus pinaster'))
discretised_plot_list_model_avg_GC$Species <- factor(discretised_plot_list_model_avg_GC$Species, levels = c('Fagus sylvatica','Quercus petraea','Betula pendula','Populus nigra','Picea abies','Pinus pinaster'))

discretised_plot_list_merge <- rbind(discretised_plot_list_model_avg_GC, discretised_plot_list_model_avg)


ggplot(discretised_plot_list, aes(fill = Species, colour = Model, x=categories, y=fraction)) +
  scale_fill_manual(values = colourlist) +
  #  scale_color_manual(values = c('black', 'gold')) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() + xlab(expression(paste("Fitness effects (", italic("N"['e']), italic("s"), ")" ) )) + ylab('Fraction of new mutations, GC conservative only') + 
  theme(legend.text = element_text(face = c("italic")))+
  coord_cartesian(ylim = c(0, 0.8))

ggplot(discretised_plot_list_merge, aes(fill = Species, colour = Mutations, x=categories, y=fraction)) +
  scale_fill_manual(values = colourlist) +
  #  scale_color_manual(values = c('black', 'gold')) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() + xlab(expression(paste("Fitness effects (", italic("N"['e']), italic("s"), ")" ) )) + ylab('Fraction of new mutations') + 
  theme(legend.text = element_text(face = c("italic")))+
  coord_cartesian(ylim = c(0, 0.8))




################################################ Suppl. Fig. 4: Discretised DFEs, comparing all SNPs to GC-conservative ################################################################

path = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/'

species_list = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster')
sp_num = 0

plot_list <- list()

for (species in species_list){
  sp_num = sp_num+1
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  
  print(sp_num)
  sp_path = paste(path, sp_short, sep = '')
  
  indep <- paste(sp_path, 'PolyDFE_all40GC_conservative.indep.out', sep = '/')
  share <- paste(sp_path, 'PolyDFE_all40GC_conservative.share.out', sep = '/')
  
  est = c(parseOutput(indep),
          parseOutput(share))
  grad = sapply(est, function(e) e$criteria)
  print(species)
  print(grad)
  print(compareModels(est[1], est[2]))
  
  colorRampPalette(c(viridis(8)[3], "white"))(4)[1]
  
  dfe_indep <- getDiscretizedDFE(est[[1]], c(-100, -10, -1, 0, 1))
  dfe_share <- getDiscretizedDFE(est[[2]], c(-100, -10, -1, 0, 1))
  
  ###The last two rows of indep_share are for models where parameters are shared across datasets, and so will be identical. I'm therefore clipping one row off here
  indep_share <- rbind(dfe_indep, dfe_share)[1:3,]
  
  i_S <- data.frame(t(indep_share))
  i_S <- cbind(c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"), i_S)
  i_S_1 <- data.frame(c(i_S[1], i_S[2]))
  i_S_1['Model'] <- "All mutations, independent"
  colnames(i_S_1) <- c('categories', 'fraction', 'Model')
  i_S_2 <- data.frame(c(i_S[1], i_S[3]))
  i_S_2['Model'] <- "GC-conservative, independent"
  colnames(i_S_2) <- c('categories', 'fraction', 'Model')
  i_S_3 <- data.frame(c(i_S[1], i_S[4]))
  i_S_3['Model'] <- "Shared"
  colnames(i_S_3) <- c('categories', 'fraction', 'Model')
  
  indep_share<- rbind(i_S_1, i_S_2, i_S_3)
  
  indep_share$categories <- factor(indep_share$categories, levels = c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  
#  print(ggplot(indep_share, aes(fill = Model, x=categories, y=fraction)) +
#          scale_fill_manual(values = c(colourlist[sp_num], colorRampPalette(c(colourlist[sp_num], "white"))(4)[3], colorRampPalette(c(colourlist[sp_num], "white"))(4)[2])) +
#          geom_bar(position="dodge", stat="identity") +
#          theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()))
  
  plot_list[[sp_num]] <-
    ggplot(indep_share, aes(fill = Model, x=categories, y=fraction)) +
    scale_fill_manual(values = c(colourlist[sp_num], colorRampPalette(c(colourlist[sp_num], "white"))(4)[3], colorRampPalette(c(colourlist[sp_num], "white"))(4)[2])) +
    geom_bar(position="dodge", stat="identity") +
    theme_classic() + theme(axis.title.x = element_blank()) + ylab('Fraction of new mutations')
  
}

a <- plot_list[[1]]
b <- plot_list[[2]]
c <- plot_list[[3]]
d <- plot_list[[4]]
e <- plot_list[[5]]
f <- plot_list[[6]]

grid.arrange(
  plot_grid(a + coord_cartesian(ylim = c(0, 0.8)) + theme(axis.text.x = element_text(angle = 90))+ theme(legend.position = c(0.7, 0.8)), #+ annotate(geom = "text", label = "*", x =4, y = 0.8, size = 12), 
            b + coord_cartesian(ylim = c(0, 0.8)) + theme(axis.text.x = element_text(angle = 90))+theme(legend.position = c(0.7, 0.8)), #+ annotate(geom = "text", label = "*", x =4, y = 0.8, size = 12), 
            c + coord_cartesian(ylim = c(0, 0.8)) + theme(axis.text.x = element_text(angle = 90))+theme(legend.position = c(0.7, 0.8)), #+ annotate(geom = "text", label = "*", x =4, y = 0.8, size = 12), 
            d + coord_cartesian(ylim = c(0, 0.8)) + theme(axis.text.x = element_text(angle = 90))+theme(legend.position = c(0.7, 0.8)),
            e + coord_cartesian(ylim = c(0, 0.8)) + theme(axis.text.x = element_text(angle = 90))+theme(legend.position = c(0.7, 0.8)), #+ annotate(geom = "text", label = "*", x =4, y = 0.8, size = 12), 
            f + coord_cartesian(ylim = c(0, 0.8)) + theme(axis.text.x = element_text(angle = 90))+theme(legend.position = c(0.7, 0.8)), 
            ncol = 3), legend, 
  layout_matrix = rbind(c(1,1,1,1,2),
                        c(1,1,1,1,2))
)




################################################ Suppl. Fig. 5: Model averaged population DFEs, deleterious only ################################################################

pop_statistics <- data.frame()
species_list = c('Betula_pendula', 'Fagus_sylvatica', 'Picea_abies', 'Pinus_pinaster', 'Populus_nigra', 'Quercus_petraea')

for (species in species_list){
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  print(species)
  path = paste('/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/', sp_short, sep = '')
  popfiles <- list.files(path = path, pattern = '_40.in.indep.out')
  pops <- lapply(popfiles, function(x) str_sub(x, 1, -11))
  
  for (pop in pops){
    plotstr <-sub("_*_40.in", "", pop)
    plotstr <-sub("PolyDFE_*", "", plotstr)
    
    Sp_pop <- paste(species, plotstr, sep = '')
    
    est_mod1 <- paste(path, '/', pop, '.4model3.pop.out', sep = '')
    est_mod2 <- paste(path, '/', pop, '.4model4.pop.out', sep = '')
    overall <- c(parseOutput(est_mod1), parseOutput(est_mod2))
    aic = getAICweights(overall)
    disc_dfe <- sapply(1:length(overall), function(i) getDiscretizedDFE(overall[[i]], c(-100, -10, -1)))
    
    disc_dfe_model_avg = apply(disc_dfe, 1, function(x) 
      sum(sapply(1:length(est), function(i) aic[i, "weight"] * x[i])) )
    
    
    disc_dfe_model_avg <- data.frame(t(disc_dfe_model_avg))
    
    eps_an <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["eps_an"]))
    theta_bar <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["theta_bar"]))
    Sd <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_d"]))
    b <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["b"]))

    prime_values <- data.frame(species, Sp_pop, plotstr, disc_dfe_model_avg, b, Sd, eps_an)
    colnames(prime_values) <- c('species', 'Sp_pop', "pop", "<-100","-100 to -10","-10 to -1","-1 to 0", 'b', 'S_d', 'eps_an')
    pop_statistics <- rbind(prime_values, pop_statistics)
    
  }
}

pop_statistics$pop[pop_statistics$Sp_pop == "Fagus_sylvaticaFRY"] <- "FR-S"
pop_statistics$pop[pop_statistics$Sp_pop == "Fagus_sylvaticaFRB"] <- "FR-C"
pop_statistics$pop[pop_statistics$Sp_pop == "Pinus_pinasterESG"] <- 	"ES-S"
pop_statistics$pop[pop_statistics$Sp_pop == "Pinus_pinasterESP"] <- 	"ES-N"
pop_statistics$pop[pop_statistics$Sp_pop == "Pinus_pinasterFRP"] <- 	"FR-W"
pop_statistics$pop[pop_statistics$Sp_pop == "Pinus_pinasterFRB"] <- 	"FR-E"
pop_statistics$pop[pop_statistics$Sp_pop == "Pinus_pinasterFRY"] <- 	"FR-C"
pop_statistics$pop[pop_statistics$Sp_pop == "Pinus_pinasterITY"] <- 	"IT-S"
pop_statistics$pop[pop_statistics$Sp_pop == "Pinus_pinasterITB"] <- 	"IT-N"
pop_statistics$pop[pop_statistics$Sp_pop == "Populus_nigraITG"] <- 	"IT-N"
pop_statistics$pop[pop_statistics$Sp_pop == "Populus_nigraITP"] <- 	"IT-S"
pop_statistics$pop[pop_statistics$Sp_pop == "Quercus_petraeaITY"] <- 	"IT-N"
pop_statistics$pop[pop_statistics$Sp_pop == "Quercus_petraeaITP"] <- 	"IT-S"

pop_statistics$species <- factor(pop_statistics$species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster'))


is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}


dat <- pop_statistics %>% tibble::rownames_to_column(var="outlier") %>% group_by(species) %>% mutate(is_outlier=ifelse(is_outlier(b), b, ''))
dat$is_outlier <- ifelse(dat$is_outlier != '', dat$pop, NA)
b_pop<- ggplot(dat, aes(x = as.factor(species), y = b)) +
  geom_boxplot(aes(x = as.factor(species), y = b, fill = species)) + 
  scale_fill_manual(values = colourlist) +
  geom_text_repel(aes(label = is_outlier), na.rm = TRUE, nudge_y = 0.5)+
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90)) 

dat <- pop_statistics %>% tibble::rownames_to_column(var="outlier") %>% group_by(species) %>% mutate(is_outlier=ifelse(is_outlier(S_d), S_d, ''))
dat$is_outlier <- ifelse(dat$is_outlier != '', dat$pop, NA)
Sd_pop<- ggplot(dat, aes(x = as.factor(species), y = S_d)) +
  geom_boxplot(aes(x = as.factor(species), y = S_d, fill = species)) + 
  scale_fill_manual(values = colourlist) +
  geom_text_repel(aes(label = is_outlier), na.rm = TRUE, nudge_y = -10)+
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90)) 

dat <- pop_statistics %>% tibble::rownames_to_column(var="outlier") %>% group_by(species) %>% mutate(is_outlier=ifelse(is_outlier(`-1 to 0`), `-1 to 0`, ''))
dat$is_outlier <- ifelse(dat$is_outlier != '', dat$pop, NA)
fraction_pop<- ggplot(dat, aes(x = as.factor(species), y = `-1 to 0`)) +
  geom_boxplot(aes(x = as.factor(species), y = `-1 to 0`, fill = species)) + 
  scale_fill_manual(values = colourlist) +
  geom_text_repel(aes(label = is_outlier), na.rm = TRUE)+
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90)) 


b_pop <- b_pop + theme(legend.position="none") + scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp"))+ theme(axis.title.y=element_text(face="italic")) 
Sd_pop <- Sd_pop + theme(legend.position="none")+ scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp")) +labs(y= "Sd") + theme(axis.title.y=element_text(face="italic"))

fraction_pop <- fraction_pop + labs(y= 'Nearly neutral fraction of mutations')+ scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp"))+ theme(legend.position="none")

grid.arrange(plot_grid(b_pop, Sd_pop, fraction_pop, ncol = 2, align = "hv", labels = c("A", "B", "C", "D"),
                       label_x = 0.05, label_y = 0.95, label_fontface = "plain"),
             legend,
             layout_matrix = rbind(c(1,1,1,2))
)




################################################ Suppl. Fig. 6: Discretised DFEs, all populations ################################################################


samplenumber = 40

pop_statistics <- data.frame()
###Pulling out the independent runs from population plots
species_list = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster')

discretised_plot_list = data.frame()

for (species in species_list){
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  print(species)
  path = paste('/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/', sp_short, sep = '')
  popfiles <- list.files(path = path, pattern = '40.in.indep.out')
  per_sp_statistics <- data.frame()
  
  for (file in popfiles) {
    print(file)
    est <-parseOutput(paste(path, file, sep = '/'))  
    print(est[[1]]$values[[2]])
    # print(est[[1]]$values[[1]])
    
    dfe <- data.frame(getDiscretizedDFE(est[[1]], c(-100, -10, -1, 0, 1)))[2,]
    values <- est[[1]]$values[[2]]
    b_val <- values["b"]
    S_d_val <- values["S_d"]
    p_b_val <- values["p_b"]
    S_b_val <- values["S_b"]
    criteria_val <- est[[1]]$criteria
    eps_an_val <- values["eps_an"]
    
    alpha_val<- est[[1]]$alpha[[2]]
    
    plotstr <-sub("_*_40.in.indep.out", "", file)
    plotstr <-sub("_*_40.in.all_indep_del.out", "", plotstr)
    plotstr <-sub("_*_40.in.all_del_indep.out", "", plotstr)
    pop <-sub("PolyDFE_*", "", plotstr)
    Sp_pop <- paste(species, pop, sep = '')
    
    prime_values <- data.frame(species, Sp_pop, file, pop, dfe, b_val, S_d_val, p_b_val, S_b_val, alpha_val, eps_an_val, criteria_val)
    colnames(prime_values) <- c('species', 'Sp_pop', 'name', 'pop', "<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1", 'b', 'S_d', 'p_b', 'S_b', 'alpha', 'eps_an', 'criteria')
    pop_statistics <- rbind(prime_values, pop_statistics)
    per_sp_statistics <- rbind(prime_values, per_sp_statistics)
  }
  
  
  
  disc_dfe <- per_sp_statistics[2:10]

  df1<- cbind(disc_dfe$Sp_pop, disc_dfe$`<-100`, 'categories' = "<-100")
  df2<- cbind(disc_dfe$Sp_pop, disc_dfe$`-100 to -10`, 'categories' = "-100 to -10")
  df3<- cbind(disc_dfe$Sp_pop, disc_dfe$`-10 to -1`, 'categories' = "-10 to -1")
  df4<- cbind(disc_dfe$Sp_pop, disc_dfe$`-1 to 0`, 'categories' = "-1 to 0")
  df5<- cbind(disc_dfe$Sp_pop, disc_dfe$`0 to 1`, 'categories' = "0 to 1")
  df6<- cbind(disc_dfe$Sp_pop, disc_dfe$`>1`, 'categories' = ">1")
  
  
  disc_dfe <- rbind(df1, df2, df3, df4, df5, df6)
  colnames(disc_dfe) <- c('Sp_pop', 'fraction', 'categories')
  
  disc_dfe <- data.frame(disc_dfe)
  disc_dfe[2] <- sapply(disc_dfe[,2], as.numeric) 
  disc_dfe <- as.data.frame(lapply(disc_dfe, unlist))
  
  disc_dfe$categories <- factor(disc_dfe$categories, levels = c("<-100", "-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  
  print(
    ggplot(disc_dfe, aes(fill = substring(gsub("[^::A-Z::]","",Sp_pop), 2), x=categories, y=fraction)) +
      scale_fill_manual(values = inferno(length(unique(disc_dfe$Sp_pop))+1)) +
      geom_bar(position="dodge", stat="identity") +
      theme_classic() + xlab(expression(paste("Fitness effects (", italic("N"['e']), italic("s"), ")" ) )) + ylab('Fraction of new mutations') +
      labs(fill = sub("_", " ", species)) +
      theme(legend.title=element_text(face="italic")) +
      theme(axis.text.x = element_text(angle = 90)) +
      coord_cartesian(ylim = c(0,1)) 
  )
  
}



################################################ Suppl. Fig. 7: All pairwise population comparisons ################################################################

Betula_pendula_Fst <- data.frame(
  c('Betula_pendula','LT-RU','LT',0.014581254412199667,0.015705156922296303),
  c('Betula_pendula','LT-SE','LT',0.02238603724244902,0.025195164770244655),
  c('Betula_pendula','LT-NO','LT',0.01658052129321898,0.017370271417814824),
  c('Betula_pendula','LT-GB','LT',0.024464175493462175,0.024399576953225737),
  c('Betula_pendula','LT-FI','LT',0.017630415650389215,0.019538834756669947),
  c('Betula_pendula','LT-FR','LT',0.02725149864793306,0.026285814693106584),
  c('Betula_pendula','LT-ES','LT',0.056395023362984964,0.05611379255014874),
  c('Betula_pendula','LT-DE','LT',0.019001707288271705,0.019166012863892634),
  c('Betula_pendula','LT-CH','LT',0.02737577165361005,0.028332353029920862),
  c('Betula_pendula','LT-IT','LT',0.05062483607548205,0.048275024101277327),
  c('Betula_pendula','RU-LT','RU',0.014581254412199665,0.0157051569222963),
  c('Betula_pendula','RU-SE','RU',0.023271776936900614,0.02468508245434476),
  c('Betula_pendula','RU-NO','RU',0.016077744550771725,0.01619520934393265),
  c('Betula_pendula','RU-GB','RU',0.023351327008634116,0.02489904842721507),
  c('Betula_pendula','RU-FI','RU',0.018484747428473516,0.019787838918940252),
  c('Betula_pendula','RU-FR','RU',0.027415354022741593,0.026669661528730394),
  c('Betula_pendula','RU-ES','RU',0.05567654033632089,0.0559045511371216),
  c('Betula_pendula','RU-DE','RU',0.019341312855573094,0.019435204367346495),
  c('Betula_pendula','RU-CH','RU',0.026251328367185977,0.027926192146055028),
  c('Betula_pendula','RU-IT','RU',0.04689863387377062,0.04588636192767228),
  c('Betula_pendula','SE-LT','SE',0.02238603724244902,0.025195164770244655),
  c('Betula_pendula','SE-RU','SE',0.023271776936900614,0.02468508245434476),
  c('Betula_pendula','SE-NO','SE',0.02034884010610283,0.021218057859175854),
  c('Betula_pendula','SE-GB','SE',0.027777924189733134,0.030576533982146006),
  c('Betula_pendula','SE-FI','SE',0.016832371482216987,0.01705175172372625),
  c('Betula_pendula','SE-FR','SE',0.028513446311998307,0.02926764151367554),
  c('Betula_pendula','SE-ES','SE',0.06002990686116735,0.0631161268840971),
  c('Betula_pendula','SE-DE','SE',0.025622649507799293,0.027500036639014815),
  c('Betula_pendula','SE-CH','SE',0.030262208321791325,0.0315819119103829),
  c('Betula_pendula','SE-IT','SE',0.049536852842160015,0.049227295444420456),
  c('Betula_pendula','NO-LT','NO',0.01658052129321898,0.017370271417814824),
  c('Betula_pendula','NO-RU','NO',0.01607774455077172,0.016195209343932645),
  c('Betula_pendula','NO-SE','NO',0.020348840106102835,0.021218057859175854),
  c('Betula_pendula','NO-GB','NO',0.019934194721350215,0.019468545958158195),
  c('Betula_pendula','NO-FI','NO',0.01713920543382517,0.019378069331372576),
  c('Betula_pendula','NO-FR','NO',0.023074241999908766,0.02152513743329376),
  c('Betula_pendula','NO-ES','NO',0.05473246813705776,0.05408330539706747),
  c('Betula_pendula','NO-DE','NO',0.018289090579337426,0.017822300125531706),
  c('Betula_pendula','NO-CH','NO',0.023480708916325612,0.02443382174169387),
  c('Betula_pendula','NO-IT','NO',0.04661798222898791,0.04707602243275029),
  c('Betula_pendula','GB-LT','GB',0.02446417549346218,0.024399576953225737),
  c('Betula_pendula','GB-RU','GB',0.023351327008634116,0.024899048427215065),
  c('Betula_pendula','GB-SE','GB',0.027777924189733134,0.030576533982146006),
  c('Betula_pendula','GB-NO','GB',0.019934194721350218,0.019468545958158195),
  c('Betula_pendula','GB-FI','GB',0.026977420757835806,0.0285380849424238),
  c('Betula_pendula','GB-FR','GB',0.01973674178665366,0.019721603803654723),
  c('Betula_pendula','GB-ES','GB',0.047525875174667026,0.046637840622334906),
  c('Betula_pendula','GB-DE','GB',0.02159957260198976,0.021186269451083862),
  c('Betula_pendula','GB-CH','GB',0.021641289781490464,0.02266926093653387),
  c('Betula_pendula','GB-IT','GB',0.05144416460321414,0.04964130478469782),
  c('Betula_pendula','FI-LT','FI',0.017630415650389215,0.01953883475666995),
  c('Betula_pendula','FI-RU','FI',0.018484747428473516,0.019787838918940252),
  c('Betula_pendula','FI-SE','FI',0.016832371482216987,0.017051751723726246),
  c('Betula_pendula','FI-NO','FI',0.017139205433825175,0.019378069331372576),
  c('Betula_pendula','FI-GB','FI',0.0269774207578358,0.0285380849424238),
  c('Betula_pendula','FI-FR','FI',0.027511561828161356,0.027973801022988094),
  c('Betula_pendula','FI-ES','FI',0.0587219587617724,0.060958405886595554),
  c('Betula_pendula','FI-DE','FI',0.022868232862222403,0.02496138637849051),
  c('Betula_pendula','FI-CH','FI',0.028994792606470168,0.030859121170448707),
  c('Betula_pendula','FI-IT','FI',0.04994002638533169,0.05013288731444898),
  c('Betula_pendula','FR-LT','FR',0.02725149864793306,0.02628581469310658),
  c('Betula_pendula','FR-RU','FR',0.02741535402274159,0.026669661528730394),
  c('Betula_pendula','FR-SE','FR',0.02851344631199831,0.02926764151367554),
  c('Betula_pendula','FR-NO','FR',0.023074241999908762,0.021525137433293765),
  c('Betula_pendula','FR-GB','FR',0.01973674178665366,0.019721603803654723),
  c('Betula_pendula','FR-FI','FR',0.027511561828161356,0.02797380102298809),
  c('Betula_pendula','FR-ES','FR',0.04381895994300434,0.04245138489958724),
  c('Betula_pendula','FR-DE','FR',0.02463900601454987,0.023846190829232724),
  c('Betula_pendula','FR-CH','FR',0.01894428524172519,0.01888503453280237),
  c('Betula_pendula','FR-IT','FR',0.049123890468593695,0.04567730301249879),
  c('Betula_pendula','ES-LT','ES',0.05639502336298497,0.05611379255014874),
  c('Betula_pendula','ES-RU','ES',0.055676540336320884,0.0559045511371216),
  c('Betula_pendula','ES-SE','ES',0.06002990686116735,0.0631161268840971),
  c('Betula_pendula','ES-NO','ES',0.05473246813705776,0.054083305397067474),
  c('Betula_pendula','ES-GB','ES',0.04752587517466703,0.0466378406223349),
  c('Betula_pendula','ES-FI','ES',0.05872195876177241,0.060958405886595554),
  c('Betula_pendula','ES-FR','ES',0.043818959943004346,0.042451384899587254),
  c('Betula_pendula','ES-DE','ES',0.05641761911393999,0.054337708945286924),
  c('Betula_pendula','ES-CH','ES',0.05899473681537678,0.05822712513847533),
  c('Betula_pendula','ES-IT','ES',0.0815390182456732,0.07633124255288008),
  c('Betula_pendula','DE-LT','DE',0.019001707288271705,0.01916601286389263),
  c('Betula_pendula','DE-RU','DE',0.019341312855573094,0.0194352043673465),
  c('Betula_pendula','DE-SE','DE',0.02562264950779929,0.027500036639014812),
  c('Betula_pendula','DE-NO','DE',0.018289090579337426,0.0178223001255317),
  c('Betula_pendula','DE-GB','DE',0.021599572601989765,0.021186269451083862),
  c('Betula_pendula','DE-FI','DE',0.022868232862222396,0.02496138637849051),
  c('Betula_pendula','DE-FR','DE',0.02463900601454987,0.023846190829232724),
  c('Betula_pendula','DE-ES','DE',0.056417619113940005,0.054337708945286924),
  c('Betula_pendula','DE-CH','DE',0.024419386973917177,0.023893977185242866),
  c('Betula_pendula','DE-IT','DE',0.05239350441447925,0.0501282989041598),
  c('Betula_pendula','CH-LT','CH',0.02737577165361005,0.028332353029920862),
  c('Betula_pendula','CH-RU','CH',0.026251328367185973,0.02792619214605502),
  c('Betula_pendula','CH-SE','CH',0.030262208321791325,0.031581911910382905),
  c('Betula_pendula','CH-NO','CH',0.023480708916325612,0.024433821741693874),
  c('Betula_pendula','CH-GB','CH',0.021641289781490464,0.022669260936533875),
  c('Betula_pendula','CH-FI','CH',0.028994792606470168,0.030859121170448707),
  c('Betula_pendula','CH-FR','CH',0.01894428524172519,0.01888503453280237),
  c('Betula_pendula','CH-ES','CH',0.05899473681537678,0.05822712513847533),
  c('Betula_pendula','CH-DE','CH',0.024419386973917177,0.023893977185242866),
  c('Betula_pendula','CH-IT','CH',0.04974154676458146,0.04771567750629942),
  c('Betula_pendula','IT-LT','IT',0.050624836075482044,0.048275024101277327),
  c('Betula_pendula','IT-RU','IT',0.04689863387377063,0.04588636192767229),
  c('Betula_pendula','IT-SE','IT',0.049536852842160015,0.04922729544442045),
  c('Betula_pendula','IT-NO','IT',0.0466179822289879,0.04707602243275029),
  c('Betula_pendula','IT-GB','IT',0.05144416460321413,0.04964130478469781),
  c('Betula_pendula','IT-FI','IT',0.04994002638533169,0.05013288731444898),
  c('Betula_pendula','IT-FR','IT',0.04912389046859369,0.04567730301249879),
  c('Betula_pendula','IT-ES','IT',0.0815390182456732,0.0763312425528801),
  c('Betula_pendula','IT-DE','IT',0.05239350441447924,0.0501282989041598),
  c('Betula_pendula','IT-CH','IT',0.049741546764581465,0.04771567750629942)
)


Fagus_sylvatica_Fst <- data.frame(
  c('Fagus_sylvatica','SE-NO','SE',0.07354518806454623,0.0725862724361456),
  c('Fagus_sylvatica','SE-GB','SE',0.058155839292920945,0.059502815600203936),
  c('Fagus_sylvatica','SE-FRY','SE',0.061454249068996956,0.06344979844021044),
  c('Fagus_sylvatica','SE-FRB','SE',0.08732557078075849,0.0832664024548793),
  c('Fagus_sylvatica','SE-ES','SE',0.09803386622883313,0.09870637241069871),
  c('Fagus_sylvatica','SE-AT','SE',0.04169673425557095,0.0405636386436924),
  c('Fagus_sylvatica','SE-DE','SE',0.032615071463570386,0.033698210878487594),
  c('Fagus_sylvatica','SE-CH','SE',0.04368654457093126,0.0444686405385181),
  c('Fagus_sylvatica','SE-SI','SE',0.03937054484021066,0.03912555969010162),
  c('Fagus_sylvatica','SE-IT','SE',0.08625171885480809,0.08477079709464014),
  c('Fagus_sylvatica','NO-SE','NO',0.07354518806454621,0.0725862724361456),
  c('Fagus_sylvatica','NO-GB','NO',0.10226092385459584,0.0971330018763892),
  c('Fagus_sylvatica','NO-FRY','NO',0.10983087013867465,0.10137610001230986),
  c('Fagus_sylvatica','NO-FRB','NO',0.1363895820878565,0.12457231388992032),
  c('Fagus_sylvatica','NO-ES','NO',0.14313057109091765,0.13155597805455227),
  c('Fagus_sylvatica','NO-AT','NO',0.09194714787708602,0.08595703333418436),
  c('Fagus_sylvatica','NO-DE','NO',0.07999436464808753,0.07630952454291917),
  c('Fagus_sylvatica','NO-CH','NO',0.09050854775506442,0.0857859671386802),
  c('Fagus_sylvatica','NO-SI','NO',0.08790658282559377,0.08274151265594265),
  c('Fagus_sylvatica','NO-IT','NO',0.13520045840848693,0.1261634299988073),
  c('Fagus_sylvatica','GB-SE','GB',0.058155839292920945,0.05950281560020392),
  c('Fagus_sylvatica','GB-NO','GB',0.10226092385459583,0.0971330018763892),
  c('Fagus_sylvatica','GB-FRY','GB',0.051311024284474274,0.05239958547664105),
  c('Fagus_sylvatica','GB-FRB','GB',0.09658643742306723,0.09676289441355182),
  c('Fagus_sylvatica','GB-ES','GB',0.08260196712747545,0.08325212389898211),
  c('Fagus_sylvatica','GB-AT','GB',0.05654295649254167,0.0547364041183244),
  c('Fagus_sylvatica','GB-DE','GB',0.0415902100904066,0.04333774790618554),
  c('Fagus_sylvatica','GB-CH','GB',0.03702648193554919,0.03762995726220434),
  c('Fagus_sylvatica','GB-SI','GB',0.053350335831040006,0.05313293899129703),
  c('Fagus_sylvatica','GB-IT','GB',0.09249900249179926,0.09430518928368357),
  c('Fagus_sylvatica','FRY-SE','FRY',0.061454249068996956,0.06344979844021047),
  c('Fagus_sylvatica','FRY-NO','FRY',0.10983087013867462,0.10137610001230986),
  c('Fagus_sylvatica','FRY-GB','FRY',0.05131102428447428,0.05239958547664105),
  c('Fagus_sylvatica','FRY-FRB','FRY',0.07992428902333555,0.07903351100560489),
  c('Fagus_sylvatica','FRY-ES','FRY',0.04995228172077263,0.04744145885868371),
  c('Fagus_sylvatica','FRY-AT','FRY',0.05725301050519581,0.054243600600443725),
  c('Fagus_sylvatica','FRY-DE','FRY',0.046708208524112636,0.04678020221341998),
  c('Fagus_sylvatica','FRY-CH','FRY',0.027946583322225178,0.02880648003050291),
  c('Fagus_sylvatica','FRY-SI','FRY',0.054411120066405236,0.052053197566639946),
  c('Fagus_sylvatica','FRY-IT','FRY',0.07107014243134113,0.07267704435403032),
  c('Fagus_sylvatica','FRB-SE','FRB',0.0873255707807585,0.0832664024548793),
  c('Fagus_sylvatica','FRB-NO','FRB',0.13638958208785656,0.12457231388992035),
  c('Fagus_sylvatica','FRB-GB','FRB',0.09658643742306723,0.09676289441355186),
  c('Fagus_sylvatica','FRB-FRY','FRB',0.07992428902333555,0.0790335110056049),
  c('Fagus_sylvatica','FRB-ES','FRB',0.11347959792536304,0.1118925020090781),
  c('Fagus_sylvatica','FRB-AT','FRB',0.06830491981540888,0.0648053873579324),
  c('Fagus_sylvatica','FRB-DE','FRB',0.06875932984235099,0.06572817741196554),
  c('Fagus_sylvatica','FRB-CH','FRB',0.0814821943453329,0.07668628510033557),
  c('Fagus_sylvatica','FRB-SI','FRB',0.06769934086769329,0.06457471173717268),
  c('Fagus_sylvatica','FRB-IT','FRB',0.04901292395817855,0.052096394621094735),
  c('Fagus_sylvatica','ES-SE','ES',0.09803386622883313,0.09870637241069868),
  c('Fagus_sylvatica','ES-NO','ES',0.14313057109091767,0.13155597805455227),
  c('Fagus_sylvatica','ES-GB','ES',0.08260196712747543,0.08325212389898211),
  c('Fagus_sylvatica','ES-FRY','ES',0.04995228172077263,0.04744145885868371),
  c('Fagus_sylvatica','ES-FRB','ES',0.11347959792536304,0.11189250200907808),
  c('Fagus_sylvatica','ES-AT','ES',0.09118769835651154,0.08651685789027709),
  c('Fagus_sylvatica','ES-DE','ES',0.08039196672586876,0.07878499569525071),
  c('Fagus_sylvatica','ES-CH','ES',0.059953037911809906,0.05890608782790217),
  c('Fagus_sylvatica','ES-SI','ES',0.09010298522228173,0.08565452036668111),
  c('Fagus_sylvatica','ES-IT','ES',0.10022995977019244,0.10155539862767433),
  c('Fagus_sylvatica','AT-SE','AT',0.041696734255570966,0.0405636386436924),
  c('Fagus_sylvatica','AT-NO','AT',0.091947147877086,0.08595703333418435),
  c('Fagus_sylvatica','AT-GB','AT',0.056542956492541666,0.0547364041183244),
  c('Fagus_sylvatica','AT-FRY','AT',0.05725301050519582,0.054243600600443725),
  c('Fagus_sylvatica','AT-FRB','AT',0.06830491981540888,0.0648053873579324),
  c('Fagus_sylvatica','AT-ES','AT',0.09118769835651155,0.08651685789027709),
  c('Fagus_sylvatica','AT-DE','AT',0.02194648427088831,0.020958327370636656),
  c('Fagus_sylvatica','AT-CH','AT',0.03920046466027391,0.035488847435520235),
  c('Fagus_sylvatica','AT-SI','AT',0.0171521307522232,0.016623691611505943),
  c('Fagus_sylvatica','AT-IT','AT',0.06516531439319811,0.06269203056059991),
  c('Fagus_sylvatica','DE-SE','DE',0.032615071463570386,0.033698210878487594),
  c('Fagus_sylvatica','DE-NO','DE',0.07999436464808753,0.07630952454291919),
  c('Fagus_sylvatica','DE-GB','DE',0.041590210090406615,0.04333774790618554),
  c('Fagus_sylvatica','DE-FRY','DE',0.04670820852411263,0.04678020221341998),
  c('Fagus_sylvatica','DE-FRB','DE',0.06875932984235097,0.06572817741196553),
  c('Fagus_sylvatica','DE-ES','DE',0.08039196672586876,0.07878499569525071),
  c('Fagus_sylvatica','DE-AT','DE',0.02194648427088831,0.020958327370636656),
  c('Fagus_sylvatica','DE-CH','DE',0.0277195444186459,0.02775556149309832),
  c('Fagus_sylvatica','DE-SI','DE',0.020228824594831295,0.020040917986144834),
  c('Fagus_sylvatica','DE-IT','DE',0.06546107502641042,0.0641790013862767),
  c('Fagus_sylvatica','CH-SE','CH',0.04368654457093125,0.044468640538518116),
  c('Fagus_sylvatica','CH-NO','CH',0.09050854775506442,0.0857859671386802),
  c('Fagus_sylvatica','CH-GB','CH',0.03702648193554919,0.03762995726220433),
  c('Fagus_sylvatica','CH-FRY','CH',0.027946583322225178,0.0288064800305029),
  c('Fagus_sylvatica','CH-FRB','CH',0.0814821943453329,0.07668628510033558),
  c('Fagus_sylvatica','CH-ES','CH',0.0599530379118099,0.05890608782790217),
  c('Fagus_sylvatica','CH-AT','CH',0.03920046466027391,0.035488847435520235),
  c('Fagus_sylvatica','CH-DE','CH',0.0277195444186459,0.02775556149309832),
  c('Fagus_sylvatica','CH-SI','CH',0.03652878928111662,0.03433614220367398),
  c('Fagus_sylvatica','CH-IT','CH',0.07488822804027642,0.07268012201005755),
  c('Fagus_sylvatica','SI-SE','SI',0.03937054484021066,0.03912555969010162),
  c('Fagus_sylvatica','SI-NO','SI',0.08790658282559377,0.08274151265594266),
  c('Fagus_sylvatica','SI-GB','SI',0.05335033583104002,0.05313293899129703),
  c('Fagus_sylvatica','SI-FRY','SI',0.054411120066405236,0.05205319756663995),
  c('Fagus_sylvatica','SI-FRB','SI',0.06769934086769329,0.06457471173717266),
  c('Fagus_sylvatica','SI-ES','SI',0.09010298522228172,0.08565452036668111),
  c('Fagus_sylvatica','SI-AT','SI',0.017152130752223202,0.016623691611505943),
  c('Fagus_sylvatica','SI-DE','SI',0.020228824594831295,0.020040917986144834),
  c('Fagus_sylvatica','SI-CH','SI',0.036528789281116615,0.03433614220367398),
  c('Fagus_sylvatica','SI-IT','SI',0.0659396439786132,0.06333523326742581),
  c('Fagus_sylvatica','IT-SE','IT',0.08625171885480809,0.08477079709464014),
  c('Fagus_sylvatica','IT-NO','IT',0.13520045840848693,0.1261634299988073),
  c('Fagus_sylvatica','IT-GB','IT',0.09249900249179928,0.09430518928368357),
  c('Fagus_sylvatica','IT-FRY','IT',0.07107014243134113,0.07267704435403032),
  c('Fagus_sylvatica','IT-FRB','IT',0.04901292395817855,0.052096394621094735),
  c('Fagus_sylvatica','IT-ES','IT',0.10022995977019244,0.10155539862767433),
  c('Fagus_sylvatica','IT-AT','IT',0.06516531439319811,0.06269203056059994),
  c('Fagus_sylvatica','IT-DE','IT',0.06546107502641042,0.0641790013862767),
  c('Fagus_sylvatica','IT-CH','IT',0.07488822804027642,0.07268012201005755),
  c('Fagus_sylvatica','IT-SI','IT',0.0659396439786132,0.06333523326742582)
)

Picea_abies_Fst <- data.frame(
  c('Picea_abies','LT-SE','LT',0.0496827124360583,0.04873507604765536),
  c('Picea_abies','LT-NO','LT',0.03368969756232694,0.03327121547440294),
  c('Picea_abies','LT-FI','LT',0.04940045355078621,0.0498727108753926),
  c('Picea_abies','LT-FR','LT',0.10073594374598821,0.10967672952806543),
  c('Picea_abies','LT-AT','LT',0.04898755988481637,0.05410205306428937),
  c('Picea_abies','LT-DE','LT',0.0482056561463583,0.05339556853579877),
  c('Picea_abies','LT-CH','LT',0.06673606599922807,0.07280803331404795),
  c('Picea_abies','LT-GR','LT',0.04645330695606679,0.046145336233545366),
  c('Picea_abies','LT-IT','LT',0.05249435022811834,0.057718903936838574),
  c('Picea_abies','SE-LT','SE',0.0496827124360583,0.048735076047655354),
  c('Picea_abies','SE-NO','SE',0.02065875981089464,0.021382378099152367),
  c('Picea_abies','SE-FI','SE',0.023886713870908682,0.023530725887228102),
  c('Picea_abies','SE-FR','SE',0.154303971555427,0.16274746293708095),
  c('Picea_abies','SE-AT','SE',0.105225772026837,0.10978652932941609),
  c('Picea_abies','SE-DE','SE',0.10425798456018075,0.10887908510859022),
  c('Picea_abies','SE-CH','SE',0.12166630048659013,0.1269507441340443),
  c('Picea_abies','SE-GR','SE',0.09792508193634764,0.09613589315160494),
  c('Picea_abies','SE-IT','SE',0.10819310921753067,0.11332010814539126),
  c('Picea_abies','NO-LT','NO',0.03368969756232694,0.03327121547440294),
  c('Picea_abies','NO-SE','NO',0.02065875981089464,0.021382378099152367),
  c('Picea_abies','NO-FI','NO',0.026392775509968526,0.026778103850046916),
  c('Picea_abies','NO-FR','NO',0.13429109043885845,0.1432427628223937),
  c('Picea_abies','NO-AT','NO',0.08375185676164994,0.0885043625758955),
  c('Picea_abies','NO-DE','NO',0.08296395072619224,0.08784773832591058),
  c('Picea_abies','NO-CH','NO',0.10033785932106788,0.10619558354571237),
  c('Picea_abies','NO-GR','NO',0.07717487790485773,0.076862496772467),
  c('Picea_abies','NO-IT','NO',0.08761885838495552,0.09178369491871279),
  c('Picea_abies','FI-LT','FI',0.04940045355078621,0.0498727108753926),
  c('Picea_abies','FI-SE','FI',0.02388671387090868,0.023530725887228106),
  c('Picea_abies','FI-NO','FI',0.02639277550996853,0.02677810385004692),
  c('Picea_abies','FI-FR','FI',0.15878995603237125,0.16774246178543237),
  c('Picea_abies','FI-AT','FI',0.10956147708195703,0.11535159116178481),
  c('Picea_abies','FI-DE','FI',0.1087994088907524,0.11461547995932225),
  c('Picea_abies','FI-CH','FI',0.12713009548119011,0.13247419011061234),
  c('Picea_abies','FI-GR','FI',0.10211742024648618,0.10262961016827743),
  c('Picea_abies','FI-IT','FI',0.11301437890437833,0.11922574958193029),
  c('Picea_abies','FR-LT','FR',0.10073594374598822,0.10967672952806543),
  c('Picea_abies','FR-SE','FR',0.154303971555427,0.16274746293708098),
  c('Picea_abies','FR-NO','FR',0.13429109043885845,0.1432427628223937),
  c('Picea_abies','FR-FI','FR',0.15878995603237128,0.16774246178543237),
  c('Picea_abies','FR-AT','FR',0.07144550475882107,0.0724917418689798),
  c('Picea_abies','FR-DE','FR',0.071306653037199,0.07452113856595453),
  c('Picea_abies','FR-CH','FR',0.07795089271687015,0.08103517111439483),
  c('Picea_abies','FR-GR','FR',0.09359124544951843,0.09706567285710675),
  c('Picea_abies','FR-IT','FR',0.06895149941729636,0.0695559557284274),
  c('Picea_abies','AT-LT','AT',0.048987559884816394,0.05410205306428936),
  c('Picea_abies','AT-SE','AT',0.10522577202683699,0.10978652932941609),
  c('Picea_abies','AT-NO','AT',0.08375185676164994,0.0885043625758955),
  c('Picea_abies','AT-FI','AT',0.10956147708195706,0.11535159116178481),
  c('Picea_abies','AT-FR','AT',0.07144550475882105,0.0724917418689798),
  c('Picea_abies','AT-DE','AT',0.014341816990602229,0.015040014614520052),
  c('Picea_abies','AT-CH','AT',0.03494722345315881,0.03356175297518092),
  c('Picea_abies','AT-GR','AT',0.04381189066890921,0.04593698060325608),
  c('Picea_abies','AT-IT','AT',0.020870429865389648,0.020100392718416123),
  c('Picea_abies','DE-LT','DE',0.0482056561463583,0.053395568535798764),
  c('Picea_abies','DE-SE','DE',0.10425798456018075,0.10887908510859022),
  c('Picea_abies','DE-NO','DE',0.08296395072619221,0.08784773832591058),
  c('Picea_abies','DE-FI','DE',0.10879940889075243,0.11461547995932227),
  c('Picea_abies','DE-FR','DE',0.07130665303719899,0.07452113856595453),
  c('Picea_abies','DE-AT','DE',0.014341816990602227,0.015040014614520049),
  c('Picea_abies','DE-CH','DE',0.03480340954181754,0.035784829109800825),
  c('Picea_abies','DE-GR','DE',0.04305395008976347,0.0456753702344851),
  c('Picea_abies','DE-IT','DE',0.02152412718437724,0.02197851253560751),
  c('Picea_abies','CH-LT','CH',0.06673606599922809,0.07280803331404799),
  c('Picea_abies','CH-SE','CH',0.12166630048659013,0.12695074413404434),
  c('Picea_abies','CH-NO','CH',0.10033785932106788,0.10619558354571239),
  c('Picea_abies','CH-FI','CH',0.12713009548119011,0.13247419011061237),
  c('Picea_abies','CH-FR','CH',0.07795089271687013,0.08103517111439483),
  c('Picea_abies','CH-AT','CH',0.03494722345315881,0.033561752975180915),
  c('Picea_abies','CH-DE','CH',0.034803409541817536,0.03578482910980083),
  c('Picea_abies','CH-GR','CH',0.060591143964958356,0.06235301046564198),
  c('Picea_abies','CH-IT','CH',0.03175861884895608,0.031329634327021),
  c('Picea_abies','GR-LT','GR',0.04645330695606679,0.04614533623354536),
  c('Picea_abies','GR-SE','GR',0.09792508193634765,0.09613589315160492),
  c('Picea_abies','GR-NO','GR',0.07717487790485773,0.076862496772467),
  c('Picea_abies','GR-FI','GR',0.10211742024648619,0.10262961016827744),
  c('Picea_abies','GR-FR','GR',0.09359124544951841,0.09706567285710675),
  c('Picea_abies','GR-AT','GR',0.04381189066890921,0.04593698060325607),
  c('Picea_abies','GR-DE','GR',0.04305395008976347,0.04567537023448509),
  c('Picea_abies','GR-CH','GR',0.060591143964958356,0.06235301046564196),
  c('Picea_abies','GR-IT','GR',0.04576819662258748,0.0477074083361337),
  c('Picea_abies','IT-LT','IT',0.052494350228118325,0.05771890393683858),
  c('Picea_abies','IT-SE','IT',0.10819310921753068,0.11332010814539127),
  c('Picea_abies','IT-NO','IT',0.08761885838495553,0.0917836949187128),
  c('Picea_abies','IT-FI','IT',0.11301437890437833,0.1192257495819303),
  c('Picea_abies','IT-FR','IT',0.06895149941729638,0.0695559557284274),
  c('Picea_abies','IT-AT','IT',0.02087042986538965,0.020100392718416127),
  c('Picea_abies','IT-DE','IT',0.02152412718437724,0.02197851253560751),
  c('Picea_abies','IT-CH','IT',0.03175861884895609,0.031329634327021),
  c('Picea_abies','IT-GR','IT',0.045768196622587465,0.04770740833613369)
)


Pinus_pinaster_Fst <- data.frame(
  c('Pinus_pinaster','ESP-FRY','ESP',0.1672121481438271,0.15021926561135648),
  c('Pinus_pinaster','ESP-FRB','ESP',0.19137505445744601,0.17144002222002783),
  c('Pinus_pinaster','ESP-MA','ESP',0.2907749271202109,0.2563501367029243),
  c('Pinus_pinaster','ESP-ITY','ESP',0.2138042977604751,0.19298592304736006),
  c('Pinus_pinaster','ESP-ITB','ESP',0.1771421015525075,0.16696349600458552),
  c('Pinus_pinaster','ESP-ESG','ESP',0.0910159688135102,0.07981335333319566),
  c('Pinus_pinaster','ESP-FRP','ESP',0.07863302504326455,0.07574790461422208),
  c('Pinus_pinaster','FRY-ESP','FRY',0.1672121481438271,0.15021926561135646),
  c('Pinus_pinaster','FRY-FRB','FRY',0.10114144432415145,0.09620831398118275),
  c('Pinus_pinaster','FRY-MA','FRY',0.3150454553622783,0.2927905609459647),
  c('Pinus_pinaster','FRY-ITY','FRY',0.07496056499506287,0.06794655753511027),
  c('Pinus_pinaster','FRY-ITB','FRY',0.07698787235467455,0.07303306721745044),
  c('Pinus_pinaster','FRY-ESG','FRY',0.17362189383512514,0.15543229363956254),
  c('Pinus_pinaster','FRY-FRP','FRY',0.23573344714974542,0.220822768900677),
  c('Pinus_pinaster','FRB-ESP','FRB',0.19137505445744601,0.17144002222002788),
  c('Pinus_pinaster','FRB-FRY','FRB',0.10114144432415145,0.09620831398118278),
  c('Pinus_pinaster','FRB-MA','FRB',0.34850527043242335,0.31111453342198997),
  c('Pinus_pinaster','FRB-ITY','FRB',0.1322734375315963,0.12107941330003906),
  c('Pinus_pinaster','FRB-ITB','FRB',0.04782259381046575,0.04472792565124972),
  c('Pinus_pinaster','FRB-ESG','FRB',0.2018227015711804,0.17589916568269048),
  c('Pinus_pinaster','FRB-FRP','FRB',0.2584432093678814,0.24164015901153243),
  c('Pinus_pinaster','MA-ESP','MA',0.2907749271202109,0.2563501367029242),
  c('Pinus_pinaster','MA-FRY','MA',0.3150454553622782,0.29279056094596473),
  c('Pinus_pinaster','MA-FRB','MA',0.34850527043242335,0.31111453342198997),
  c('Pinus_pinaster','MA-ITY','MA',0.35215814408345597,0.3267045759270273),
  c('Pinus_pinaster','MA-ITB','MA',0.3326124597200716,0.30441394144654566),
  c('Pinus_pinaster','MA-ESG','MA',0.1726951517775408,0.159170455666004),
  c('Pinus_pinaster','MA-FRP','MA',0.33836351967880796,0.3043984429061587),
  c('Pinus_pinaster','ITY-ESP','ITY',0.2138042977604751,0.19298592304736012),
  c('Pinus_pinaster','ITY-FRY','ITY',0.07496056499506289,0.06794655753511027),
  c('Pinus_pinaster','ITY-FRB','ITY',0.1322734375315963,0.12107941330003907),
  c('Pinus_pinaster','ITY-MA','ITY',0.352158144083456,0.3267045759270273),
  c('Pinus_pinaster','ITY-ITB','ITY',0.09553550615766218,0.07940947883845048),
  c('Pinus_pinaster','ITY-ESG','ITY',0.21390922966244344,0.1927529101526048),
  c('Pinus_pinaster','ITY-FRP','ITY',0.27983847779125137,0.2583340353917367),
  c('Pinus_pinaster','ITB-ESP','ITB',0.1771421015525075,0.1669634960045855),
  c('Pinus_pinaster','ITB-FRY','ITB',0.07698787235467455,0.07303306721745044),
  c('Pinus_pinaster','ITB-FRB','ITB',0.04782259381046575,0.04472792565124972),
  c('Pinus_pinaster','ITB-MA','ITB',0.3326124597200716,0.3044139414465457),
  c('Pinus_pinaster','ITB-ITY','ITB',0.09553550615766217,0.07940947883845048),
  c('Pinus_pinaster','ITB-ESG','ITB',0.18472137530118118,0.16867237656826894),
  c('Pinus_pinaster','ITB-FRP','ITB',0.2423079251738297,0.22977072746344188),
  c('Pinus_pinaster','ESG-ESP','ESG',0.0910159688135102,0.07981335333319564),
  c('Pinus_pinaster','ESG-FRY','ESG',0.1736218938351252,0.15543229363956254),
  c('Pinus_pinaster','ESG-FRB','ESG',0.2018227015711804,0.1758991656826905),
  c('Pinus_pinaster','ESG-MA','ESG',0.17269515177754086,0.15917045566600393),
  c('Pinus_pinaster','ESG-ITY','ESG',0.21390922966244344,0.19275291015260484),
  c('Pinus_pinaster','ESG-ITB','ESG',0.18472137530118118,0.16867237656826897),
  c('Pinus_pinaster','ESG-FRP','ESG',0.14685424381693557,0.1300339267597454),
  c('Pinus_pinaster','FRP-ESP','FRP',0.07863302504326455,0.07574790461422208),
  c('Pinus_pinaster','FRP-FRY','FRP',0.23573344714974542,0.220822768900677),
  c('Pinus_pinaster','FRP-FRB','FRP',0.2584432093678814,0.24164015901153243),
  c('Pinus_pinaster','FRP-MA','FRP',0.3383635196788079,0.30439844290615864),
  c('Pinus_pinaster','FRP-ITY','FRP',0.27983847779125137,0.2583340353917367),
  c('Pinus_pinaster','FRP-ITB','FRP',0.2423079251738297,0.22977072746344188),
  c('Pinus_pinaster','FRP-ESG','FRP',0.1468542438169356,0.1300339267597454)
)

Populus_nigra_Fst <- data.frame(
  c('Populus_nigra','ITP-GB','ITP',0.21083930071456855,0.21594881620010806),
  c('Populus_nigra','ITP-MA','ITP',0.4881430230619548,0.49498161694275544),
  c('Populus_nigra','ITP-FR','ITP',0.1005581391896736,0.09706676915874281),
  c('Populus_nigra','ITP-DE','ITP',0.13461068070038754,0.13695922434895685),
  c('Populus_nigra','ITP-CH','ITP',0.12032135360437385,0.12013932486921049),
  c('Populus_nigra','ITP-ITG','ITP',0.09955173446855817,0.09854891763592294),
  c('Populus_nigra','GB-ITP','GB',0.21083930071456855,0.21594881620010808),
  c('Populus_nigra','GB-MA','GB',0.5485905905684131,0.5575112374157495),
  c('Populus_nigra','GB-FR','GB',0.15795608934128458,0.166324002589936),
  c('Populus_nigra','GB-DE','GB',0.1785508536877693,0.1839276809664831),
  c('Populus_nigra','GB-CH','GB',0.20668162333314982,0.21113642003880964),
  c('Populus_nigra','GB-ITG','GB',0.20116225784454878,0.20078085291682704),
  c('Populus_nigra','MA-ITP','MA',0.4881430230619548,0.49498161694275544),
  c('Populus_nigra','MA-GB','MA',0.5485905905684132,0.5575112374157494),
  c('Populus_nigra','MA-FR','MA',0.4581064784027005,0.46318674183282404),
  c('Populus_nigra','MA-DE','MA',0.5086630498479411,0.5220025618220947),
  c('Populus_nigra','MA-CH','MA',0.5197345186794073,0.5259908132017024),
  c('Populus_nigra','MA-ITG','MA',0.5118905804884255,0.5215268599678176),
  c('Populus_nigra','FR-ITP','FR',0.10055813918967359,0.09706676915874278),
  c('Populus_nigra','FR-GB','FR',0.15795608934128458,0.166324002589936),
  c('Populus_nigra','FR-MA','FR',0.4581064784027005,0.46318674183282404),
  c('Populus_nigra','FR-DE','FR',0.1106065177054086,0.1168305249841997),
  c('Populus_nigra','FR-CH','FR',0.112888895099491,0.11168631283217381),
  c('Populus_nigra','FR-ITG','FR',0.10700679361476005,0.1046311034412361),
  c('Populus_nigra','DE-ITP','DE',0.13461068070038754,0.13695922434895688),
  c('Populus_nigra','DE-GB','DE',0.1785508536877693,0.1839276809664831),
  c('Populus_nigra','DE-MA','DE',0.5086630498479411,0.5220025618220947),
  c('Populus_nigra','DE-FR','DE',0.11060651770540861,0.1168305249841997),
  c('Populus_nigra','DE-CH','DE',0.0853750646508007,0.08573235368365734),
  c('Populus_nigra','DE-ITG','DE',0.07646577340150132,0.07389864330841515),
  c('Populus_nigra','CH-ITP','CH',0.12032135360437385,0.12013932486921051),
  c('Populus_nigra','CH-GB','CH',0.20668162333314985,0.21113642003880964),
  c('Populus_nigra','CH-MA','CH',0.5197345186794073,0.5259908132017024),
  c('Populus_nigra','CH-FR','CH',0.112888895099491,0.11168631283217383),
  c('Populus_nigra','CH-DE','CH',0.0853750646508007,0.08573235368365735),
  c('Populus_nigra','CH-ITG','CH',0.03967861268707009,0.03841758125060023),
  c('Populus_nigra','ITG-ITP','ITG',0.09955173446855817,0.09854891763592294),
  c('Populus_nigra','ITG-GB','ITG',0.20116225784454878,0.2007808529168271),
  c('Populus_nigra','ITG-MA','ITG',0.5118905804884255,0.5215268599678176),
  c('Populus_nigra','ITG-FR','ITG',0.10700679361476005,0.1046311034412361),
  c('Populus_nigra','ITG-DE','ITG',0.07646577340150133,0.07389864330841515),
  c('Populus_nigra','ITG-CH','ITG',0.03967861268707009,0.03841758125060024)
)

Quercus_petraea_Fst<- data.frame(
  c('Quercus_petraea','ITP-LT','ITP',0.12521238799984194,0.11550553304566111),
  c('Quercus_petraea','ITP-NO','ITP',0.12041859296968853,0.11368895130029838),
  c('Quercus_petraea','ITP-GB','ITP',0.1190959182325433,0.11330424304499435),
  c('Quercus_petraea','ITP-FR','ITP',0.11714937274987447,0.11378090563422102),
  c('Quercus_petraea','ITP-ES','ITP',0.1013244163715,0.09553365183738514),
  c('Quercus_petraea','ITP-ITY','ITP',0.11193583218123014,0.10901207325096672),
  c('Quercus_petraea','ITP-DE','ITP',0.08270335437999281,0.08100646524534087),
  c('Quercus_petraea','ITP-CH','ITP',0.11557320628783255,0.11019762368453553),
  c('Quercus_petraea','LT-ITP','LT',0.1252123879998419,0.11550553304566114),
  c('Quercus_petraea','LT-NO','LT',0.07089984646771424,0.0674812607750665),
  c('Quercus_petraea','LT-GB','LT',0.07264453209990153,0.06984179256317054),
  c('Quercus_petraea','LT-FR','LT',0.0750521038279894,0.07188548834725617),
  c('Quercus_petraea','LT-ES','LT',0.0705800744835302,0.06647210755448188),
  c('Quercus_petraea','LT-ITY','LT',0.07964436266343634,0.07912835855591228),
  c('Quercus_petraea','LT-DE','LT',0.06027905209039644,0.058207183777950824),
  c('Quercus_petraea','LT-CH','LT',0.07173126795398482,0.06947642742355686),
  c('Quercus_petraea','NO-ITP','NO',0.12041859296968854,0.11368895130029838),
  c('Quercus_petraea','NO-LT','NO',0.07089984646771424,0.0674812607750665),
  c('Quercus_petraea','NO-GB','NO',0.023212899979955348,0.023837666934076976),
  c('Quercus_petraea','NO-FR','NO',0.029872387200077174,0.03171849268288767),
  c('Quercus_petraea','NO-ES','NO',0.02714775721602412,0.029662884380597585),
  c('Quercus_petraea','NO-ITY','NO',0.04376987310397582,0.04533878820300899),
  c('Quercus_petraea','NO-DE','NO',0.033993187941809534,0.033958614390641954),
  c('Quercus_petraea','NO-CH','NO',0.027646339103084585,0.02786389925634021),
  c('Quercus_petraea','GB-ITP','GB',0.1190959182325433,0.11330424304499435),
  c('Quercus_petraea','GB-LT','GB',0.07264453209990153,0.06984179256317055),
  c('Quercus_petraea','GB-NO','GB',0.023212899979955348,0.023837666934076972),
  c('Quercus_petraea','GB-FR','GB',0.032227338031250094,0.03229767679730283),
  c('Quercus_petraea','GB-ES','GB',0.026705133215237927,0.02759267692297064),
  c('Quercus_petraea','GB-ITY','GB',0.04792397365933747,0.0479167402799965),
  c('Quercus_petraea','GB-DE','GB',0.03825793362045434,0.03726449192934688),
  c('Quercus_petraea','GB-CH','GB',0.030887641643239934,0.02912692115857181),
  c('Quercus_petraea','FR-ITP','FR',0.11714937274987444,0.11378090563422102),
  c('Quercus_petraea','FR-LT','FR',0.0750521038279894,0.07188548834725617),
  c('Quercus_petraea','FR-NO','FR',0.029872387200077177,0.03171849268288767),
  c('Quercus_petraea','FR-GB','FR',0.03222733803125009,0.03229767679730283),
  c('Quercus_petraea','FR-ES','FR',0.021468142365495116,0.02194426159619893),
  c('Quercus_petraea','FR-ITY','FR',0.026618790978035933,0.027208924173810354),
  c('Quercus_petraea','FR-DE','FR',0.02838559943058567,0.028951914926040653),
  c('Quercus_petraea','FR-CH','FR',0.018464393891237017,0.018663880782903874),
  c('Quercus_petraea','ES-ITP','ES',0.1013244163715,0.09553365183738514),
  c('Quercus_petraea','ES-LT','ES',0.07058007448353022,0.06647210755448188),
  c('Quercus_petraea','ES-NO','ES',0.027147757216024117,0.02966288438059759),
  c('Quercus_petraea','ES-GB','ES',0.026705133215237924,0.027592676922970634),
  c('Quercus_petraea','ES-FR','ES',0.021468142365495112,0.02194426159619893),
  c('Quercus_petraea','ES-ITY','ES',0.030571530693649655,0.031949412562886303),
  c('Quercus_petraea','ES-DE','ES',0.027788400428661007,0.02807150895164682),
  c('Quercus_petraea','ES-CH','ES',0.022061224867299594,0.023401782907187432),
  c('Quercus_petraea','ITY-ITP','ITY',0.11193583218123013,0.10901207325096672),
  c('Quercus_petraea','ITY-LT','ITY',0.07964436266343636,0.07912835855591228),
  c('Quercus_petraea','ITY-NO','ITY',0.04376987310397582,0.045338788203008996),
  c('Quercus_petraea','ITY-GB','ITY',0.04792397365933747,0.0479167402799965),
  c('Quercus_petraea','ITY-FR','ITY',0.026618790978035933,0.027208924173810347),
  c('Quercus_petraea','ITY-ES','ITY',0.03057153069364965,0.031949412562886303),
  c('Quercus_petraea','ITY-DE','ITY',0.03198169617373433,0.03319979389530078),
  c('Quercus_petraea','ITY-CH','ITY',0.029282895528400212,0.03127737760333512),
  c('Quercus_petraea','DE-ITP','DE',0.0827033543799928,0.08100646524534086),
  c('Quercus_petraea','DE-LT','DE',0.060279052090396454,0.058207183777950824),
  c('Quercus_petraea','DE-NO','DE',0.033993187941809534,0.03395861439064197),
  c('Quercus_petraea','DE-GB','DE',0.038257933620454335,0.03726449192934689),
  c('Quercus_petraea','DE-FR','DE',0.028385599430585662,0.028951914926040653),
  c('Quercus_petraea','DE-ES','DE',0.027788400428661007,0.02807150895164682),
  c('Quercus_petraea','DE-ITY','DE',0.031981696173734336,0.03319979389530078),
  c('Quercus_petraea','DE-CH','DE',0.02660709856916925,0.026970450898949654),
  c('Quercus_petraea','CH-ITP','CH',0.1155732062878325,0.11019762368453553),
  c('Quercus_petraea','CH-LT','CH',0.07173126795398482,0.06947642742355686),
  c('Quercus_petraea','CH-NO','CH',0.027646339103084585,0.027863899256340214),
  c('Quercus_petraea','CH-GB','CH',0.03088764164323993,0.02912692115857181),
  c('Quercus_petraea','CH-FR','CH',0.018464393891237017,0.01866388078290387),
  c('Quercus_petraea','CH-ES','CH',0.022061224867299587,0.023401782907187432),
  c('Quercus_petraea','CH-ITY','CH',0.02928289552840021,0.03127737760333512),
  c('Quercus_petraea','CH-DE','CH',0.026607098569169242,0.026970450898949654)
)

species_list = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster')
sp_num = 0
plot_list = list()

for (species in species_list){
  
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  sp_num = sp_num+1
  Fst <- mget(paste(species,'_Fst', sep = ''))
  
  Fst <- data.frame(t(Fst[[1]]))
  colnames(Fst) <- c('Species','pop', 'focal_pop', 'Fst', 'Fst0')
  Fst[4:5] <- sapply(Fst[,4:5], as.numeric) 
  
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  sp_path = paste('/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/',sp_short, sep = '')
  
  resultsfiles <- list.files(path = sp_path, pattern = '_40.in.indep.out')

  differences <- data.frame()
  
  for (focus_file in resultsfiles){
    focpop <-sub("*_40.in.indep.out", "", focus_file)
    
    focpop <-sub("PolyDFE_*", "", focpop)
    complist <- resultsfiles[which(resultsfiles != focus_file)]
    for (file in complist){
      est <- parseOutput(paste(sp_path, file, sep = '/'))
      comp <- est[[1]]$values[[2]]
      
      est <- parseOutput(paste(sp_path, focus_file, sep = '/'))
      focus <- est[[1]]$values[[2]]
      
      focus_values <- c(focus["b"], focus["S_d"], focus["p_b"], focus["S_b"], est[[1]]$alpha[[1]], focus["eps_an"], est[[1]]$criteria)  
      comp_values <-  c(comp["b"], comp["S_d"], comp["p_b"], comp["S_b"], est[[1]]$alpha[[2]], comp["eps_an"], 0)    
      diffs <- focus_values-comp_values
      pop <-sub("*_40.in.indep.out", "", file)
      pop <-sub("PolyDFE_*", "", pop)
      pop <- paste(focpop, pop, sep = '-')
      diffs <- c(pop, diffs)
      differences <- rbind(differences, diffs)
      
    }
  }
  
  colnames <- c('pop', 'b', 'S_d', 'p_b', 'S_b', 'alpha', 'eps_an', 'criteria')
  colnames(differences) <- colnames
  
  diff_Fst <- merge(differences, Fst, by = "pop")

  
  ####Flag- remember there is a mistake in alpha, cannot compare across until fixed
  ####And good to do this with nearly neutral fraction too? But simpler to use pi and pin/pis
  
  ### because values are compared through subtraction, removing negative values removes duplicates
  
  Sd_plot <- ggplot(data =  diff_Fst[which(diff_Fst$S_d >= 0),])+
    geom_point(aes(x = as.numeric(Fst), y = as.numeric(S_d), fill = focal_pop), shape=23, size=4, fill = colourlist[sp_num])+
    #    scale_fill_manual(values = inferno(length(unique(diff_Fst[which(diff_Fst$S_d >= 0),]$focal_pop))))+
    geom_smooth(aes(x = as.numeric(Fst), y = as.numeric(S_d)), method = "lm") +
    xlab(expression("F"['ST'])) + ylab("Difference in Sd") + theme_bw()
  
  b_plot <- ggplot(data =  diff_Fst[which(diff_Fst$b >= 0),])+
    geom_point(aes(x = as.numeric(Fst), y = as.numeric(b), fill = focal_pop), shape=23, size=4, fill = colourlist[sp_num])+
    #    scale_fill_manual(values = inferno(length(unique(diff_Fst[which(diff_Fst$b >= 0),]$focal_pop))))+
    geom_smooth(aes(x = as.numeric(Fst), y = as.numeric(b)), method = "lm")+
    xlab(expression("F"['ST'])) + ylab("Difference in b") + theme_bw()

  ### Do the same with pi values - pin/pis
  print(species)
  pi <- read.table('/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PiCalculate.txt', header = TRUE)
  Sp_pi <- pi[which(pi$Species == species),]
  rownames(Sp_pi) <- c(Sp_pi$Sp_pop) 
  
  focal_pop_list <- Sp_pi$Sp_pop
  
  pi_diff = data.frame()
  
  for (pop in focal_pop_list){
    focal_pop <- sub(species,"",pop)
    print(focal_pop)
    diff <- c(Sp_pi$pi0fold_to_pi4fold - Sp_pi[pop, 'pi0fold_to_pi4fold'])
    comp_pop <- c(sub(species,"",Sp_pi$Sp_pop))
    pop <- paste(focal_pop, comp_pop, sep =  '-')
    pi_ratio_diff <- cbind(diff, pop)
    pi_diff <- rbind(pi_diff, pi_ratio_diff)
  }
  
  pi_diff <- pi_diff[which(pi_diff$diff != 0),]
  
  ### merge with Fst
  pi_diff_fst <- merge(Fst, pi_diff, by = "pop")
  ### because pi is compared, removing negative values removes duplicates.
  pi_diff_fst <- pi_diff_fst[which(pi_diff_fst$diff >= 0),] 
  
  pi_plot <- ggplot(data = pi_diff_fst)+
    geom_point(aes(x = as.numeric(Fst), y = abs(as.numeric(diff)), fill = focal_pop), shape=23, size=4, fill = colourlist[sp_num])+
    #    scale_fill_manual(values = inferno(length(unique(pi_diff_fst$focal_pop))))+
    geom_smooth(aes(x = as.numeric(Fst), y = abs(as.numeric(diff))), method = "lm") +
    xlab(expression("F"['ST'])) + ylab(expression(paste("Difference in ", pi[0]," /", pi[4])) ) + theme_bw()
  
  
  plot_list[[sp_num]] <- grid.arrange(b_plot + theme(legend.position="none"), 
                                      Sd_plot + theme(legend.position="none"), 
                                      pi_plot+ theme(legend.position="none"),
                                      layout_matrix = rbind(c(1,2,3))
  )
  
  
  print(summary(lm(diff_Fst[which(diff_Fst$b >= 0),]$b ~ diff_Fst[which(diff_Fst$b >= 0),]$Fst)))
  print(summary(lm(diff_Fst[which(diff_Fst$S_d >= 0),]$S_d ~ diff_Fst[which(diff_Fst$S_d >= 0),]$Fst)))
  print(summary(lm(as.numeric(pi_diff_fst$diff) ~ pi_diff_fst$Fst)))
  
  
}

a <- plot_list[[1]]
b <- plot_list[[2]]
c <- plot_list[[3]]
d <- plot_list[[4]]
e <- plot_list[[5]]
f <- plot_list[[6]]


plot_grid(a, b, c, d, e, f, legend, ncol = 2)

grid.arrange(plot_grid(a, b, c, d, e, f, ncol = 2),
             legend,
             layout_matrix = rbind(c(1,1,1,1,2))
)

