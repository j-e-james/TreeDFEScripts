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

################################################ Figure 1: Pi, Fst, Rxy ################################################################

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


###Fst calculated using python script, FstCalculate.py
Fst <- data.frame(
  c('Betula_pendula','Betula_pendulaCH','DE','CH',0.024419386973917177,0.023893977185242866),
  c('Betula_pendula','Betula_pendulaES','DE','ES',0.05641761911393999,0.054337708945286924),
  c('Betula_pendula','Betula_pendulaFI','DE','FI',0.022868232862222403,0.02496138637849051),
  c('Betula_pendula','Betula_pendulaFR','DE','FR',0.02463900601454987,0.023846190829232724),
  c('Betula_pendula','Betula_pendulaGB','DE','GB',0.02159957260198976,0.021186269451083862),
  c('Betula_pendula','Betula_pendulaIT','DE','IT',0.05239350441447924,0.0501282989041598),
  c('Betula_pendula','Betula_pendulaLT','DE','LT',0.019001707288271705,0.019166012863892634),
  c('Betula_pendula','Betula_pendulaNO','DE','NO',0.018289090579337426,0.017822300125531706),
  c('Betula_pendula','Betula_pendulaRU','DE','RU',0.019341312855573094,0.019435204367346495),
  c('Betula_pendula','Betula_pendulaSE','DE','SE',0.025622649507799293,0.027500036639014815),
  c('Fagus_sylvatica','Fagus_sylvaticaAT','CH','AT',0.03920046466027391,0.035488847435520235),
  c('Fagus_sylvatica','Fagus_sylvaticaDE','CH','DE',0.0277195444186459,0.02775556149309832),
  c('Fagus_sylvatica','Fagus_sylvaticaES','CH','ES',0.059953037911809906,0.05890608782790217),
  c('Fagus_sylvatica','Fagus_sylvaticaFRY','CH','FRY',0.027946583322225178,0.02880648003050291),
  c('Fagus_sylvatica','Fagus_sylvaticaFRB','CH','FRB',0.0814821943453329,0.07668628510033557),
  c('Fagus_sylvatica','Fagus_sylvaticaGB','CH','GB',0.03702648193554919,0.03762995726220434),
  c('Fagus_sylvatica','Fagus_sylvaticaIT','CH','IT',0.07488822804027642,0.07268012201005755),
  c('Fagus_sylvatica','Fagus_sylvaticaNO','CH','NO',0.09050854775506442,0.0857859671386802),
  c('Fagus_sylvatica','Fagus_sylvaticaSE','CH','SE',0.04368654457093126,0.0444686405385181),
  c('Fagus_sylvatica','Fagus_sylvaticaSI','CH','SI',0.036528789281116615,0.03433614220367398),
  c('Picea_abies','Picea_abiesAT','LT','AT',0.048987559884816394,0.05410205306428936),
  c('Picea_abies','Picea_abiesCH','LT','CH',0.06673606599922809,0.07280803331404799),
  c('Picea_abies','Picea_abiesDE','LT','DE',0.0482056561463583,0.053395568535798764),
  c('Picea_abies','Picea_abiesFI','LT','FI',0.04940045355078621,0.0498727108753926),
  c('Picea_abies','Picea_abiesFR','LT','FR',0.10073594374598822,0.10967672952806543),
  c('Picea_abies','Picea_abiesGR','LT','GR',0.04645330695606679,0.04614533623354536),
  c('Picea_abies','Picea_abiesIT','LT','IT',0.052494350228118325,0.05771890393683858),
  c('Picea_abies','Picea_abiesNO','LT','NO',0.03368969756232694,0.03327121547440294),
  c('Picea_abies','Picea_abiesSE','LT','SE',0.0496827124360583,0.048735076047655354),
  c('Pinus_pinaster','Pinus_pinasterESG','FRY','ESG',0.1736218938351252,0.15543229363956254),
  c('Pinus_pinaster','Pinus_pinasterESP','FRY','ESP',0.1672121481438271,0.15021926561135648),
  c('Pinus_pinaster','Pinus_pinasterFRP','FRY','FRP',0.23573344714974542,0.220822768900677),
  c('Pinus_pinaster','Pinus_pinasterFRB','FRY','FRB',0.10114144432415145,0.09620831398118278),
  c('Pinus_pinaster','Pinus_pinasterITY','FRY','ITY',0.07496056499506289,0.06794655753511027),
  c('Pinus_pinaster','Pinus_pinasterITB','FRY','ITB',0.07698787235467455,0.07303306721745044),
  c('Pinus_pinaster','Pinus_pinasterMA','FRY','MA',0.3150454553622782,0.29279056094596473),
  c('Pinus_sylvestris','Pinus_sylvestrisCH','DE','CH',0.02670962913985506,0.028680820607451824),
  c('Pinus_sylvestris','Pinus_sylvestrisES','DE','ES',0.04069141252362084,0.039904425486197895),
  c('Pinus_sylvestris','Pinus_sylvestrisFR','DE','FR',0.023700855382000816,0.027199592888275566),
  c('Pinus_sylvestris','Pinus_sylvestrisGB','DE','GB',0.03105738802823849,0.024911843095293432),
  c('Pinus_sylvestris','Pinus_sylvestrisGR','DE','GR',0.026474274672395916,0.02353985712357884),
  c('Pinus_sylvestris','Pinus_sylvestrisIT','DE','IT',0.07387766464877403,0.06680448325291945),
  c('Pinus_sylvestris','Pinus_sylvestrisLT','DE','LT',0.015699059004022284,0.013476315087318879),
  c('Pinus_sylvestris','Pinus_sylvestrisNO','DE','NO',0.015687100823942,0.014578670935053797),
  c('Pinus_sylvestris','Pinus_sylvestrisRU','DE','RU',0.030955396625602125,0.023324411232285173),
  c('Populus_nigra','Populus_nigraCH','ITG','CH',0.03967861268707009,0.03841758125060023),
  c('Populus_nigra','Populus_nigraDE','ITG','DE',0.07646577340150132,0.07389864330841515),
  c('Populus_nigra','Populus_nigraFR','ITG','FR',0.10700679361476005,0.1046311034412361),
  c('Populus_nigra','Populus_nigraGB','ITG','GB',0.20116225784454878,0.20078085291682704),
  c('Populus_nigra','Populus_nigraITP','ITG','ITP',0.09955173446855817,0.09854891763592294),
  c('Populus_nigra','Populus_nigraMA','ITG','MA',0.5118905804884255,0.5215268599678176),
  c('Quercus_petraea','Quercus_petraeaDE','CH','DE',0.02660709856916925,0.026970450898949654),
  c('Quercus_petraea','Quercus_petraeaES','CH','ES',0.022061224867299594,0.023401782907187432),
  c('Quercus_petraea','Quercus_petraeaFR','CH','FR',0.018464393891237017,0.018663880782903874),
  c('Quercus_petraea','Quercus_petraeaGB','CH','GB',0.030887641643239934,0.02912692115857181),
  c('Quercus_petraea','Quercus_petraeaITY','CH','ITY',0.029282895528400212,0.03127737760333512),
  c('Quercus_petraea','Quercus_petraeaITP','CH','ITP',0.11557320628783255,0.11019762368453553),
  c('Quercus_petraea','Quercus_petraeaLT','CH','LT',0.07173126795398482,0.06947642742355686),
  c('Quercus_petraea','Quercus_petraeaNO','CH','NO',0.027646339103084585,0.02786389925634021)
)  

Fst <- data.frame(t(Fst))
colnames(Fst) <- c('Species', 'Sp_pop','Ref','Focal', 'Fst', 'Fst0')
Fst[5:6] <- sapply(Fst[,5:6], as.numeric) 

Lat <- read.table('/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/Latitude_pi_average.txt', header = TRUE)

pi <- read.table('/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PiCalculate.txt', header = TRUE)

pi_lat <- merge(pi, Lat)
pi_lat <- pi_lat[which(pi_lat$Species != 'Pinus_sylvestris'),]

pi_lat$Species <- sub("_", " ", pi_lat$Species)
pi_lat$Species <- factor(pi_lat$Species, levels = c('Fagus sylvatica','Quercus petraea','Betula pendula','Populus nigra','Picea abies','Pinus pinaster'))

pi_ratio_plot <- ggplot(data = pi_lat, aes(x = Latitude , y = pi0fold_to_pi4fold, fill = Species))+
  geom_smooth (method='lm', aes(col = Species), alpha = 0.05, linetype = "dotted") +
  scale_color_manual(values = colourlist) +
  geom_point(size = 3, shape=21) + 
  scale_fill_manual(values = colourlist) +
  theme_classic() + xlab('Latitude') + ylab(expression(paste(pi[0]," /", pi[4]))) +
  theme(legend.text = element_text(face = c("italic")))+
  theme(legend.title= element_blank())

pi_lat <- merge(pi, Lat)
pi_lat <- merge(pi_lat, Fst)
pi_lat <- pi_lat[which(pi_lat$Species != 'Pinus_sylvestris'),]

pi_lat$Species <- sub("_", " ", pi_lat$Species)
pi_lat$Species <- factor(pi_lat$Species, levels = c('Fagus sylvatica','Quercus petraea','Betula pendula','Populus nigra','Picea abies','Pinus pinaster'))


Fst_plot <- ggplot(data = pi_lat, aes(x = Latitude , y = Fst, fill = Species))+
  geom_smooth (method='lm', aes(col = Species), alpha = 0.05, linetype = "dotted") +
  scale_color_manual(values = colourlist) +
  geom_point(size = 3, shape=21) + 
  scale_fill_manual(values = colourlist) +
  theme(legend.text = element_text(face = c("italic"))) +
  theme_classic() + xlab('Latitude (degrees)') + ylab(expression(paste("F"["ST"]))) + coord_cartesian(ylim = c(0, 0.5))


###calculating SE based on 100 jackknife samples, used formula that is quite standard, of form: 
###https://blogs.sas.com/content/iml/2017/06/21/jackknife-estimate-standard-error-sas.html#:~:text=The%20jackknife%20method%20estimates%20the,jackknife%20samples%20from%20the%20data.
path = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/Rxy_to_centralpop_norm/'
sp_list = c('Betula_pendula', 'Fagus_sylvatica', 'Picea_abies', 'Pinus_pinaster', 'Populus_nigra', 'Quercus_petraea')
species_list = c('Bpendula', 'Fsylvatica', 'Pabies', 'Ppinaster', 'Pnigra', 'Qpetraea')

Rxy <- data.frame()


### I have changed the Rxy file name format, so this code must change - will delineate progress
num = 1
for (species in species_list){
  files <- c(list.files(path = path, pattern = species))
  print(files)
  for (file in files){
    popfile = paste(path, file, sep = '')
    
    popfile <- read.table(popfile, skip = 1)
    popvalue <- read.table(paste(path, file, sep = ''), nrows = 1)
    colnames(popvalue) <- c('Rxy', 'Ryx', 'Rxy0', 'Ryx0', 'Rxy0norm', 'Ryx0norm')
    
    
    colnames(popfile) <- c('JKNum', 'Rx', 'Ry', 'Rxy', 'Rx0', 'Ry0', 'Rxy0', 'Rxy0norm', 'Ryx0norm')
    Sp <- sub('.+_(.+)', '\\1', file)
    Sp <- sub('*.tsv', "", Sp)
    
    i_function <- function(JKRxy){(JKRxy - Tavg)**2 }
    
    Tavg <- mean(popfile$Rxy) 
    i_sum <- sapply(popfile$Rxy, FUN = i_function)
    Tavg <- mean(popfile$Rxy0) 
    i_sum0 <- sapply(popfile$Rxy0, FUN = i_function)
    Tavg <- mean(popfile$Rxy0/popfile$Rxy) 
    i_sum0over4 <- sapply(popfile$Rxy0/popfile$Rxy, FUN = i_function)
    Tavg <- mean(popfile$Rxy0norm) 
    i_sum0norm <- sapply(popfile$Rxy0norm, FUN = i_function)     
    
    
    SEjack <- sqrt( (100-1)/100 * sum(i_sum) )
    SEjack0 <- sqrt( (100-1)/100 * sum(i_sum0) )
    SEjack0over4 <- sqrt( (100-1)/100 * sum(i_sum0over4) )
    SEjack0norm <- sqrt( (100-1)/100 * sum(i_sum0norm) )
    
    popfile$Species <- sp_list[num] 
    popfile$Sp_pop <- paste(sp_list[num], Sp, sep = '')
    popfile$Sp <- Sp
    popfile$Rxy_value <- popvalue$Rxy
    popfile$Rxy0_value <- popvalue$Rxy0
    popfile$Rxy0norm_value <- popvalue$Rxy0norm
    
    
    popfile$SEjack <- SEjack
    popfile$SEjack0 <- SEjack0
    popfile$SEjack0over4 <- SEjack0over4
    popfile$SEjack0norm <- SEjack0norm
    
    Rxy <- rbind(popfile, Rxy)
    
  }
  num = num+1
}


Rxy <- merge(pi_lat, Rxy, by = 'Sp_pop')
Rxy$Species.x <- Rxy$Species

Rxy$Species <- factor(Rxy$Species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster'))

Rxy_order <- Rxy[order(Rxy$Species, as.factor(Rxy$Latitude)),]

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


Rxy_order[Rxy_order == "Fagus_sylvaticaFRY"] <- "Fagus_sylvaticaFR-South"
Rxy_order[Rxy_order == "Fagus_sylvaticaFRB"] <- "Fagus_sylvaticaFR-Corsica"
Rxy_order[Rxy_order == "Pinus_pinasterESG"] <- 	"Pinus_pinasterES-South"
Rxy_order[Rxy_order == "Pinus_pinasterESP"] <- 	"Pinus_pinasterES-North"
Rxy_order[Rxy_order == "Pinus_pinasterFRP"] <- 	"Pinus_pinasterFR-West"
Rxy_order[Rxy_order == "Pinus_pinasterFRB"] <- 	"Pinus_pinasterFR-East"
Rxy_order[Rxy_order == "Pinus_pinasterFRY"] <- 	"Pinus_pinasterFR-Corsica"
Rxy_order[Rxy_order == "Pinus_pinasterITY"] <- 	"Pinus_pinasterIT-South"
Rxy_order[Rxy_order == "Pinus_pinasterITB"] <- 	"Pinus_pinasterIT-North"
Rxy_order[Rxy_order == "Populus_nigraITG"] <- 	"Populus_nigraIT-North"
Rxy_order[Rxy_order == "Populus_nigraITP"] <- 	"Populus_nigraIT-South"
Rxy_order[Rxy_order == "Quercus_petraeaITY"] <- 	"Quercus_petraeaIT-North"
Rxy_order[Rxy_order == "Quercus_petraeaITP"] <- 	"Quercus_petraeaIT-South"


Rxy_order$Sp_pop <- factor(Rxy_order$Sp_pop, levels = unique(Rxy_order$Sp_pop))
pop_labels <- as.list(levels(Rxy_order$Sp_pop))


Rxy0_plot <- ggplot(data = Rxy_order) +
  geom_errorbar(aes(x = as.factor(Sp_pop), ymin = Rxy0_value - 2*SEjack0, ymax = Rxy0_value+ 2*SEjack0, color = Species), size = 0.5) +
  scale_color_manual(values =colourlist) +
  #  geom_violin(aes(x = as.factor(Sp_pop), y = Rxy0, fill = Species), scale = "width") + 
  geom_point(aes(x = as.factor(Sp_pop), y = Rxy0_value), shape=23, size=2, colour = "white", fill = "black")+
  scale_fill_manual(values = colourlist) + 
#  scale_x_discrete(labels = substring(gsub("[^::A-Z::]","",pop_labels), 2)) +
  scale_x_discrete(labels =str_extract( pop_labels, '[A-Z][A-Z].?.?')) +
  geom_hline(yintercept = 1, linetype = 'dashed')  +
  labs( y = expression(paste("R"['XY']," 0-fold degenerate")) ) +
  theme_bw() + theme(axis.title.x = element_blank() ) + theme(axis.text.x = element_text(angle = 90)) + coord_cartesian(ylim = c(0.2, 1.8))

Rxy_norm_plot <- ggplot(data = Rxy_order) +
  geom_errorbar(aes(x = as.factor(Sp_pop), ymin = Rxy0norm_value - 2*SEjack0norm, ymax = Rxy0norm_value + 2*SEjack0norm, color = Species), size = 0.5) +
  scale_color_manual(values = colourlist) +
  #  geom_violin(aes(x = as.factor(Sp_pop), y = Rxy, fill = Species), scale = "width") + 
  geom_point(aes(x = as.factor(Sp_pop), y = Rxy0norm_value), shape=23, size=2, colour = "white", fill = "black")+
  scale_fill_manual(values = colourlist) + 
  scale_x_discrete(labels =str_extract( pop_labels, '[A-Z][A-Z].?.?')) +
  geom_hline(yintercept = 1, linetype = 'dashed')  +
  labs( y = expression(paste("R'"['XY']," 0-fold degenerate")) ) +
  theme_bw() + theme(axis.title.x = element_blank() ) + theme(axis.text.x = element_text(angle = 90)) + coord_cartesian(ylim = c(0.2, 1.8))


legend <- get_legend(pi_ratio_plot)
grid.arrange(pi_ratio_plot + theme(legend.position="none") + annotate(geom = "text", label = "A", x = 65, y = 0.38, size = 7), 
             Fst_plot + theme(legend.position="none") + annotate(geom = "text", label = "B", x = 65, y = 0.495, size = 7),  
             legend,
             Rxy0_plot + theme(legend.position="none")+ annotate(geom = "text", label = "C", x = 48, y = 1.6, size = 7), 
             Rxy_norm_plot + theme(legend.position="none")+ annotate(geom = "text", label = "D",  x = 48, y = 1.6, size = 7), 
             layout_matrix = rbind(c(1,1,1, 2, 2, 2, 3),
                                   c(1,1,1, 2, 2, 2, 3),
                                   c(4,4, 4, 4,4,4, 4 ),
                                   c(5,5,5,5,5,5,5)))





################################################ Figure 2: Model averaged, deleterious-only DFE ################################################################


boot_statistics = data.frame()

species_list = c('Betula_pendula', 'Fagus_sylvatica', 'Picea_abies', 'Pinus_pinaster', 'Populus_nigra', 'Quercus_petraea')
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
  
  boot_filepath = paste(path, 'Bootstrap', sep = '/')
  dfedisc <- list()
  dflist <- list()
  S_d <- list()
  b <- list()
  S_b <- list()
  p_b <- list()
  alpha <- list()
  criteria <- list()
  eps_an <- list()
  
  for (val in 1:100){
    
    est_1_boot = parseOutput(paste(boot_filepath, '/Boot_', val,'_4model3.out', sep = ''))
    est_2_boot = parseOutput(paste(boot_filepath, '/Boot_', val,'_4model4.out', sep = ''))
    
    boot_dfe1 = data.frame(getDiscretizedDFE(est_1_boot[[1]], c(-100, -10, -1)))
    boot_dfe2 = data.frame(getDiscretizedDFE(est_2_boot[[1]], c(-100, -10, -1)))

    overall_boot <- c(est_1_boot, est_2_boot) #, est_3_boot, est_4_boot)
    aic_boot = getAICweights(overall_boot)
    
    disc_dfe <- sapply(1:length(overall_boot), function(i) getDiscretizedDFE(overall_boot[[i]], c(-100, -10, -1)))
    disc_dfe_model_avg = sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * disc_dfe[1,][i]))
    
    disc_dfe_model_avg = apply(disc_dfe, 1, function(x) 
      sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * x[i])) )
    
    dfe_disc <- t(data.frame(disc_dfe_model_avg))
    colnames(dfe_disc) <- c("val.100","val100.10","val10.1","val1.0")
    dfedisc <- rbind(dfedisc, data.frame(dfe_disc))
    
    
    eps_an <- c(eps_an, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["eps_an"])))
    S_d <- c(S_d, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["S_d"])))
    b <- c(b, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["b"])))
    p_b <- c(p_b, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["p_b"])))
    S_b <- c(S_b, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["S_b"])))

  }
  
  ###these steps occur per species
  statistics = data.frame(species, dfedisc, unlist(S_d), unlist(b), unlist(S_b), unlist(p_b), unlist(eps_an))
  colnames(statistics) <- c('Species', "<-100","-100 to -10","-10 to -1","-1 to 0", "S_d","b","S_b","p_b", "eps_an")
  boot_statistics = rbind(boot_statistics, statistics)
  
  upper<- list()
  lower<- list()
  for(i in 1:ncol(dfedisc)) {       # for-loop over columns
    upper <- c(upper, sort(dfedisc[ , i] ,partial=100-2)[100-2])
    lower <- c(lower, sort(dfedisc[ , i] ,partial=3)[3])
  }
  
  est_mod1 = paste(path, 'PolyDFE_4model3_all40.out', sep = '/')
  est_mod2 = paste(path, 'PolyDFE_4model4_all40.out', sep = '/')
  
  est1 <- parseOutput(est_mod1)
  dfe1 <- data.frame(getDiscretizedDFE(est1[[1]], c(-100, -10, -1)))
  est2 <- parseOutput(est_mod2)
  dfe2 <- data.frame(getDiscretizedDFE(est2[[1]], c(-100, -10, -1)))

  overall <- c(parseOutput(est_mod1), parseOutput(est_mod2))#, parseOutput(est_mod3), parseOutput(est_mod4))
  
  aic = getAICweights(overall)
  disc_dfe <- sapply(1:length(overall), function(i) getDiscretizedDFE(overall[[i]], c(-100, -10, -1)))
  disc_dfe_model_avg = sum(sapply(1:length(overall), function(i) aic[i, "weight"] * disc_dfe[1,][i]))
  
  disc_dfe_model_avg = apply(disc_dfe, 1, function(x) 
    sum(sapply(1:length(overall), function(i) aic[i, "weight"] * x[i])) )
  
  
  disc_dfe_model_avg <- data.frame(t(disc_dfe_model_avg))
  
  dfe <- rbind(disc_dfe_model_avg, upper, lower)
  
  dfe <- data.frame(t(dfe))
  dfe <- cbind(dfe, c("<-100","-100 to -10","-10 to -1","-1 to 0"))
  dfe <- cbind(dfe, c(species))
  colnames(dfe) <- c('fraction', 'upper', 'lower', 'categories', 'species')
  dfe$categories <- factor(dfe$categories, levels = dfe$categories)
  
  discretised_plot_list_model_avg = rbind(discretised_plot_list_model_avg, dfe)
  
  
  eps_an <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["eps_an"]))
  theta_bar <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["theta_bar"]))
  Sd <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_d"]))
  b <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["b"]))
  pb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["p_b"]))
  Sb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_b"]))

  prime_values <- data.frame(species, plotstr, disc_dfe_model_avg, b, Sd, pb, Sb, eps_an)
  colnames(prime_values) <- c('Species', 'name', "<-100","-100 to -10","-10 to -1","-1 to 0", 'b', 'S_d', 'p_b', 'S_b', 'eps_an')
  model_average_statistics <- rbind(prime_values, model_average_statistics)
  
  
  ###we might want to visually compare the discretised DFEs from the different models, which is where this comes in
  dfe1 <- data.frame(t(dfe1))
  dfe1 <- cbind(dfe1, c("<-100","-100 to -10","-10 to -1","-1 to 0"))
  dfe1 <- cbind(dfe1, c(species), c('mod1'))
  
  dfe2 <- data.frame(t(dfe2))
  dfe2 <- cbind(dfe2, c("<-100","-100 to -10","-10 to -1","-1 to 0"))
  dfe2 <- cbind(dfe2, c(species), c('mod2'))

  colnames(dfe1) <- c('fraction', 'categories', 'Species', 'Model')
  dfe1$categories <- factor(dfe1$categories, levels = dfe1$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe1)
  
  colnames(dfe2) <- c('fraction', 'categories', 'Species', 'Model')
  dfe2$categories <- factor(dfe2$categories, levels = dfe2$categories)
  discretised_plot_list = rbind(discretised_plot_list, dfe2)
  
  num = num + 1
}


discretised_plot_list_model_avg$species <- factor(discretised_plot_list_model_avg$species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster','Pinus_sylvestris'))


best_models_disc_DFE<- ggplot(discretised_plot_list_model_avg, aes(fill = species, x=categories, y=fraction)) +
  scale_fill_manual(values = colourlist) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2, position= position_dodge(0.9)) +
  theme_classic() + xlab(expression(paste("Fitness effects (", italic("N"['e']), italic("s"), ")" ) )) + ylab('Fraction of new mutations')

boot_statistics$Species <- factor(boot_statistics$Species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster','Pinus_sylvestris'))
model_average_statistics$Species <- factor(model_average_statistics$Species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster','Pinus_sylvestris'))
best_boots <- boot_statistics
best_primes <- model_average_statistics


###remove the values above and below 95% confidence intervals
tops <- best_boots %>% group_by(Species) %>% slice_min(b, n = 98)
tops <- tops %>% group_by(Species) %>% slice_max(b, n = 98)
b_plot <- ggplot()+
  geom_violin(data = tops, aes(factor(Species), b, fill = Species), scale = "width") + scale_fill_manual(values = colourlist) +
  geom_point(data = best_primes, aes(factor(Species), b), shape=23, size=4, colour = "white", fill = "black") + 
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90))

tops <- best_boots %>% group_by(Species) %>% slice_min(S_d, n = 98)
tops <- tops %>% group_by(Species) %>% slice_max(S_d, n = 98)
Sd_plot <- ggplot()+
  geom_violin(data = tops, aes(factor(Species), S_d, fill = Species), scale = "width") + scale_fill_manual(values = colourlist) +
  geom_point(data = best_primes, aes(factor(Species), S_d), shape=23, size=4, colour = "white", fill = "black") + 
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90))

legend <- get_legend(pi_ratio_plot)

b <- b_plot + theme(legend.position="none") + scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp")) + labs(y= "b")+theme(axis.title.y=element_text(face="italic")) 
Sd <- Sd_plot + theme(legend.position="none")+ scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp")) +labs(y= "Sd") + theme(axis.title.y=element_text(face="italic"))


grid.arrange(best_models_disc_DFE  + theme(legend.position="none") + annotate(geom = "text", label = "A", x = 4.5, y = 0.8, size = 7),
             legend,
             ###plot_grid is from cowplot, and allows aligning
             plot_grid(b, Sd, ncol = 2, align = "hv", labels = c("B", "C"), label_x = 0.06, label_y = 0.95, label_fontface = "plain"),
             layout_matrix = rbind(c(1,1,1,1,2),
                                   c(1,1,1,1,2),
                                   c(3,3,3,3,2)))



################################################ Figure 3: All DFE model results, model averaged ################################################################


boot_statistics = data.frame()

species_list = c('Betula_pendula', 'Fagus_sylvatica', 'Picea_abies', 'Pinus_pinaster', 'Populus_nigra', 'Quercus_petraea')
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
  
  boot_filepath = paste(path, 'Bootstrap', sep = '/')
  dfedisc <- list()
  dflist <- list()
  S_d <- list()
  b <- list()
  S_b <- list()
  p_b <- list()
  alpha <- list()
  criteria <- list()
  eps_an <- list()
  
  for (val in 1:100){
    
    est_1_boot = parseOutput(paste(boot_filepath, '/Boot_', val,'_4model1.out', sep = ''))
    est_2_boot = parseOutput(paste(boot_filepath, '/Boot_', val,'_4model2.out', sep = ''))
    est_3_boot = parseOutput(paste(boot_filepath, '/Boot_', val,'_4model3.out', sep = ''))
    est_4_boot = parseOutput(paste(boot_filepath, '/Boot_', val,'_4model4.out', sep = ''))
    
    boot_dfe1 = data.frame(getDiscretizedDFE(est_1_boot[[1]], c(-100, -10, -1, 0, 1)))
    boot_dfe2 = data.frame(getDiscretizedDFE(est_2_boot[[1]], c(-100, -10, -1, 0, 1)))
    boot_dfe3 = data.frame(getDiscretizedDFE(est_3_boot[[1]], c(-100, -10, -1, 0, 1)))
    boot_dfe4 = data.frame(getDiscretizedDFE(est_4_boot[[1]], c(-100, -10, -1, 0, 1)))
    
    overall_boot <- c(est_1_boot, est_2_boot, est_3_boot, est_4_boot)
    aic_boot = getAICweights(overall_boot)
    
    disc_dfe <- sapply(1:length(overall_boot), function(i) getDiscretizedDFE(overall_boot[[i]], c(-100, -10, -1, 0, 1)))
    disc_dfe_model_avg = sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * disc_dfe[1,][i]))
    
    disc_dfe_model_avg = apply(disc_dfe, 1, function(x) 
      sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * x[i])) )
    
    dfe_disc <- t(data.frame(disc_dfe_model_avg))
    colnames(dfe_disc) <- c("val.100","val100.10","val10.1","val1.0","val0.1","val.1")
    dfedisc <- rbind(dfedisc, data.frame(dfe_disc))
    
    eps_an <- c(eps_an, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["eps_an"])))
    S_d <- c(S_d, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["S_d"])))
    b <- c(b, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["b"])))
    p_b <- c(p_b, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["p_b"])))
    S_b <- c(S_b, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * overall_boot[[i]]$values[[1]]["S_b"])))
    
    alpha<- c(alpha, sum(sapply(1:length(overall_boot), function(i) aic_boot[i, "weight"] * as.numeric(overall_boot[[i]]$alpha[1]))))
    
  }
  
  ###these steps occur per species
  statistics = data.frame(species, dfedisc, unlist(S_d), unlist(b), unlist(S_b), unlist(p_b), unlist(alpha), unlist(eps_an))
  colnames(statistics) <- c('Species', "<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1", "S_d","b","S_b","p_b", "alpha", "eps_an")
  boot_statistics = rbind(boot_statistics, statistics)
  
  upper<- list()
  lower<- list()
  for(i in 1:ncol(dfedisc)) {       # for-loop over columns
    upper <- c(upper, sort(dfedisc[ , i] ,partial=100-2)[100-2])
    lower <- c(lower, sort(dfedisc[ , i] ,partial=3)[3])
  }
  
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
  
  overall <- c(parseOutput(est_mod1), parseOutput(est_mod2), parseOutput(est_mod3), parseOutput(est_mod4))
  
  aic = getAICweights(overall)
  disc_dfe <- sapply(1:length(overall), function(i) getDiscretizedDFE(overall[[i]], c(-100, -10, -1, 0, 1)))
  disc_dfe_model_avg = sum(sapply(1:length(overall), function(i) aic[i, "weight"] * disc_dfe[1,][i]))
  
  disc_dfe_model_avg = apply(disc_dfe, 1, function(x) 
    sum(sapply(1:length(overall), function(i) aic[i, "weight"] * x[i])) )
  
  
  disc_dfe_model_avg <- data.frame(t(disc_dfe_model_avg))
  
  dfe <- rbind(disc_dfe_model_avg, upper, lower)
  
  dfe <- data.frame(t(dfe))
  dfe <- cbind(dfe, c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  dfe <- cbind(dfe, c(species))
  colnames(dfe) <- c('fraction', 'upper', 'lower', 'categories', 'species')
  dfe$categories <- factor(dfe$categories, levels = dfe$categories)
  
  discretised_plot_list_model_avg = rbind(discretised_plot_list_model_avg, dfe)
  
  
  eps_an <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["eps_an"]))
  theta_bar <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["theta_bar"]))
  Sd <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_d"]))
  b <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["b"]))
  pb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["p_b"]))
  Sb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_b"]))
  alpha_dfe <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * as.numeric(overall[[i]]$alpha[1])))
  
  prime_values <- data.frame(species, plotstr, disc_dfe_model_avg, b, Sd, pb, Sb, alpha_dfe, eps_an)
  colnames(prime_values) <- c('Species', 'name', "<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1", 'b', 'S_d', 'p_b', 'S_b', 'alpha', 'eps_an')
  model_average_statistics <- rbind(prime_values, model_average_statistics)
  
  
  ###we might want to visually compare teh discretised DFEs from the different models, which is where this comes in
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


discretised_plot_list_model_avg$species <- factor(discretised_plot_list_model_avg$species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster'))

### incorporate pi0/pi4 estimates per species
pi_estimates <- read.table('/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PiCalculateSpecies.txt', header = TRUE)
pi_estimates <- pi_estimates[which(pi_estimates$species != 'Pinus_sylvestris'),]

discretised_plot_list_mod1 <- discretised_plot_list[which(discretised_plot_list$Model == 'mod1'),]
discretised_plot_list_mod3 <- discretised_plot_list[which(discretised_plot_list$Model == 'mod3'),]
discretised_plot_list_mod1$species <- factor(discretised_plot_list_mod1$species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster'))
discretised_plot_list_mod3$species <- factor(discretised_plot_list_mod3$species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster'))

slightly_del_fraction <- merge(discretised_plot_list_mod1[which(discretised_plot_list_mod1$categories == "-1 to 0"),], pi_estimates, by.x = 'Species', by.y = 'species')
slightly_ben_fraction <- merge(discretised_plot_list_mod1[which(discretised_plot_list_mod1$categories == "0 to 1"),], pi_estimates, by.x = 'Species', by.y = 'species')

slightly_del_fraction['ben_fraction'] <- slightly_ben_fraction$fraction
slightly_del_fraction$species <- factor(slightly_del_fraction$Species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster','Pinus_sylvestris'))
slightly_del_fraction_mod1 <- slightly_del_fraction

slightly_del_fraction <- merge(discretised_plot_list_mod3[which(discretised_plot_list_mod3$categories == "-1 to 0"),], pi_estimates, by.x = 'Species', by.y = 'species')
slightly_del_fraction$species <- factor(slightly_del_fraction$Species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster','Pinus_sylvestris'))
slightly_del_fraction_mod3 <- slightly_del_fraction

shapes <- c("Full model" = 18, "Deleterious-ony model" = 19)

pi_ratio_slightlydel_fraction <- ggplot()+
  geom_point(data = slightly_del_fraction_mod1, aes(fraction+ben_fraction, pi0fold_to_pi4fold, color = species, shape = "Full model"), size = 4) +
  geom_point(data = slightly_del_fraction_mod3, aes(fraction, pi0fold_to_pi4fold, color = species, shape = "Deleterious-ony model"), size = 3) +
  scale_color_manual(values = colourlist, guide = "none") +
  labs(x = "Nearly neutral fraction of mutations",
       y = expression(paste(pi[0]," /", pi[4])),
       shape = NULL) +
  scale_shape_manual(values = shapes) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  coord_cartesian(xlim = c(0, 0.4), ylim = c(0.1, 0.5)) +
  theme_classic() 



best_models_disc_DFE<- ggplot(discretised_plot_list_model_avg, aes(fill = species, x=categories, y=fraction)) +
  scale_fill_manual(values = colourlist) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2, position= position_dodge(0.9)) +
  theme_classic() + xlab(expression(paste("Fitness effects (", italic("N"['e']), italic("s"), ")" ) )) + ylab('Fraction of new mutations')

boot_statistics$Species <- factor(boot_statistics$Species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster','Pinus_sylvestris'))
model_average_statistics$Species <- factor(model_average_statistics$Species, levels = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster','Pinus_sylvestris'))
best_boots <- boot_statistics
best_primes <- model_average_statistics


###remove the values above and below 95% confidence intervals
tops <- best_boots %>% group_by(Species) %>% slice_min(b, n = 98)
tops <- tops %>% group_by(Species) %>% slice_max(b, n = 98)
b_plot <- ggplot()+
  geom_violin(data = tops, aes(factor(Species), b, fill = Species), scale = "width") + scale_fill_manual(values = colourlist) +
  geom_point(data = best_primes, aes(factor(Species), b), shape=23, size=4, colour = "white", fill = "black") + 
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90))

tops <- best_boots %>% group_by(Species) %>% slice_min(S_d, n = 98)
tops <- tops %>% group_by(Species) %>% slice_max(S_d, n = 98)
Sd_plot <- ggplot()+
  geom_violin(data = tops, aes(factor(Species), S_d, fill = Species), scale = "width") + scale_fill_manual(values = colourlist) +
  geom_point(data = best_primes, aes(factor(Species), S_d), shape=23, size=4, colour = "white", fill = "black") + 
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90))

tops <- best_boots %>% group_by(Species) %>% slice_min(alpha, n = 98)
tops <- tops %>% group_by(Species) %>% slice_max(alpha, n = 98)
alpha_plot<- ggplot()+
  geom_point(data = best_primes, aes(factor(Species), alpha), shape=23, size=4, colour = "white", fill = "black") + 
  geom_violin(data = tops, aes(factor(Species), alpha, fill = Species), scale = "width") + scale_fill_manual(values = colourlist) +
  geom_point(data = best_primes, aes(factor(Species), alpha), shape=23, size=4, colour = "white", fill = "black") + 
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90))

legend <- get_legend(pi_ratio_plot)

b <- b_plot + theme(legend.position="none") + scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp")) +labs(y= "b")+ theme(axis.title.y=element_text(face="italic")) 
Sd <- Sd_plot + theme(legend.position="none") + scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp")) +labs(y= "Sd") + theme(axis.title.y=element_text(face="italic"))
alpha <- alpha_plot + theme(legend.position="none")+ scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp")) + labs(y= expression(paste(alpha['DFE'])) ) + theme(axis.title.y=element_text(face="italic")) 
pi_ratio <- pi_ratio_slightlydel_fraction + theme(legend.position = c(1.1,0.3)) + theme(legend.text = element_text(size=8)) + theme(legend.key.size = unit(0.01, 'cm'))

grid.arrange(best_models_disc_DFE + coord_cartesian(ylim = c(0,1))  + theme(legend.position="none") + annotate(geom = "text", label = "A", x = 6, y = 0.8, size = 7),
             legend,
             ###plot_grid is from cowplot, and allows aligning
             plot_grid(b, Sd, alpha, pi_ratio, ncol = 2, align = "hv", labels = c("B", "C", "D", "E"), label_x = 0.06, label_y = 0.95, label_fontface = "plain"),
             layout_matrix = rbind(c(1,1,1,1,2),
                                   c(3,3,3,3,2)))




################################################ Figure 4: Discretised DFEs, comparing all genes to all-species orthologs ################################################################

path = '/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/'

species_list = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster')
sp_num = 0

plot_list <- list()

for (species in species_list){
  sp_num = sp_num+1
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  
  print(sp_num)
  sp_path = paste(path, sp_short, sep = '')
  
  indep <- paste(sp_path, 'PolyDFE_all40orthogroups_allspecies.indep.out', sep = '/')
  share <- paste(sp_path, 'PolyDFE_all40orthogroups_allspecies.share.out', sep = '/')
  
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
  i_S_1['Model'] <- "Full genome, independent fit"
  colnames(i_S_1) <- c('categories', 'fraction', 'Model')
  i_S_2 <- data.frame(c(i_S[1], i_S[3]))
  i_S_2['Model'] <- "Orthologs, independent fit"
  colnames(i_S_2) <- c('categories', 'fraction', 'Model')
  i_S_3 <- data.frame(c(i_S[1], i_S[4]))
  i_S_3['Model'] <- "Shared"
  colnames(i_S_3) <- c('categories', 'fraction', 'Model')
  
  indep_share<- rbind(i_S_1, i_S_2, i_S_3)
  
  indep_share$categories <- factor(indep_share$categories, levels = c("<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1"))
  
  print(ggplot(indep_share, aes(fill = Model, x=categories, y=fraction)) +
          scale_fill_manual(values = c(colourlist[sp_num], colorRampPalette(c(colourlist[sp_num], "white"))(4)[3], colorRampPalette(c(colourlist[sp_num], "white"))(4)[2])) +
          geom_bar(position="dodge", stat="identity") +
          theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()))
  
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




################################################ Analysis: Compare species independent and shared models for genome partitions ################################################################

### for comparisons of independent-shared models fitted to GC conservative  as opposed to all mutations
species_list = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster')

for (species in species_list){
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  indep <- paste("/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/",sp_short,"/PolyDFE_all40GC_conservative.indep.out", sep = '')
  share <- paste("/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/",sp_short,"/PolyDFE_all40GC_conservative.share.out", sep = '')
  est = c(parseOutput(indep),
          parseOutput(share))
  grad = sapply(est, function(e) e$criteria)
  print(species)
  print(grad)
  print(compareModels(est[1], est[2]))
}

### for comparisons of independent-shared models fitted to all species orthologs  as opposed to all genes
species_list = c('Fagus_sylvatica','Quercus_petraea','Betula_pendula','Populus_nigra','Picea_abies','Pinus_pinaster')

for (species in species_list){
  sp_short = paste(substr(species,1,1), str_split(species, '_')[[1]][2], sep = '')
  indep <- paste("/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/",sp_short,"/PolyDFE_all40orthogroups_allspecies.indep.out", sep = '')
  share <- paste("/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/",sp_short,"/PolyDFE_all40orthogroups_allspecies.share.out", sep = '')
  est = c(parseOutput(indep),
          parseOutput(share))
  grad = sapply(est, function(e) e$criteria)
  print(species)
  print(grad)
  print(compareModels(est[1], est[2]))
}



################################################ Figure 5: Model averaged population DFEs ################################################################

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
    
    est_mod1 <- paste(path, '/', pop, '.4model1.pop.out', sep = '')
    est_mod2 <- paste(path, '/', pop, '.4model2.pop.out', sep = '')
    est_mod3 <- paste(path, '/', pop, '.4model3.pop.out', sep = '')
    est_mod4 <- paste(path, '/', pop, '.4model4.pop.out', sep = '')
    overall <- c(parseOutput(est_mod1), parseOutput(est_mod2), parseOutput(est_mod3), parseOutput(est_mod4))
    aic = getAICweights(overall)
    disc_dfe <- sapply(1:length(overall), function(i) getDiscretizedDFE(overall[[i]], c(-100, -10, -1, 0, 1)))
    
    disc_dfe_model_avg = apply(disc_dfe, 1, function(x) 
      sum(sapply(1:length(est), function(i) aic[i, "weight"] * x[i])) )
    
    
    disc_dfe_model_avg <- data.frame(t(disc_dfe_model_avg))
    
    eps_an <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["eps_an"]))
    theta_bar <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["theta_bar"]))
    Sd <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_d"]))
    b <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["b"]))
    pb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["p_b"]))
    Sb <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * overall[[i]]$values[[1]]["S_b"]))
    alpha_dfe <- sum(sapply(1:length(overall), function(i) aic[i, "weight"] * as.numeric(overall[[i]]$alpha[1])))
    
    prime_values <- data.frame(species, Sp_pop, plotstr, disc_dfe_model_avg, b, Sd, pb, Sb, alpha_dfe, eps_an)
    colnames(prime_values) <- c('species', 'Sp_pop', "pop", "<-100","-100 to -10","-10 to -1","-1 to 0","0 to 1",">1", 'b', 'S_d', 'p_b', 'S_b', 'alpha', 'eps_an')
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

dat <- pop_statistics %>% tibble::rownames_to_column(var="outlier") %>% group_by(species) %>% mutate(is_outlier=ifelse(is_outlier(alpha), alpha, ''))
dat$is_outlier <- ifelse(dat$is_outlier != '', dat$pop, NA)
alpha_pop<- ggplot(dat, aes(x = as.factor(species), y = alpha)) +
  geom_boxplot(aes(x = as.factor(species), y = alpha, fill = species)) + 
  scale_fill_manual(values = colourlist) +
  geom_text_repel(aes(label = is_outlier), na.rm = TRUE, nudge_y = -0.05)+
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90)) 

dat <- pop_statistics %>% tibble::rownames_to_column(var="outlier") %>% group_by(species) %>% mutate(is_outlier=ifelse(is_outlier(`0 to 1` + `-1 to 0`), `0 to 1` + `-1 to 0`, ''))
dat$is_outlier <- ifelse(dat$is_outlier != '', dat$pop, NA)
fraction_pop<- ggplot(dat, aes(x = as.factor(species), y = `0 to 1` + `-1 to 0`)) +
  geom_boxplot(aes(x = as.factor(species), y = `0 to 1` + `-1 to 0`, fill = species)) + 
  scale_fill_manual(values = colourlist) +
  geom_text_repel(aes(label = is_outlier), na.rm = TRUE)+
  theme_classic() + xlab("Species") + theme(axis.text.x = element_text(angle = 90)) 


b_pop <- b_pop + theme(legend.position="none") + scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp"))+ theme(axis.title.y=element_text(face="italic")) 
Sd_pop <- Sd_pop + theme(legend.position="none")+ scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp")) +labs(y= "Sd") + theme(axis.title.y=element_text(face="italic"))
alpha_pop <- alpha_pop + theme(legend.position="none")+ scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp")) + labs(y= expression(paste(alpha['DFE'])) ) +theme(axis.title.y=element_text(face="italic")) 

fraction_pop <- fraction_pop + labs(y= 'Nearly neutral fraction of mutations')+ scale_x_discrete(labels=c("Fs", "Qp", "Bp", "Pn", "Pa", "Pp"))+ theme(legend.position="none")

grid.arrange(plot_grid(b_pop, Sd_pop, alpha_pop, fraction_pop, ncol = 2, align = "hv", labels = c("A", "B", "C", "D"),
                       label_x = 0.05, label_y = 0.95, label_fontface = "plain"),
             legend,
             layout_matrix = rbind(c(1,1,1,2))
)




################################################ Final discussion analysis: compare the DFEs of more closely related species ################################################################

compareModels("/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/FsylQpet_indep.out", "/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/FsylQpet_share.out")
compareModels("/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/FsylQpet_indep.out", "/Users/jennyjames/Dropbox/LASCOUX/GenTree/assign_ancestral/PolyDFE_GenTree_Clean/FsylQpet_shareb.out")


