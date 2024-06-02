# Progetto sulle Foreste Casentinesi
# Elaborerò alcune immagini, di questo Parco Nazionale, recuperate dal satellite Sentinel 2.
# L'obiettivo è indagare l'impatto del cambiamento climatico e delle attività umane 
# sulla salute vegetativa e la copertura del suolo nelle Foreste Casentinesi (dal 2018 al 2023).

# Installazione pacchetti necessari:

# install.packages("raster")
# install.packages("viridis")
# install.packages("ggplot2")
# install.packages("tidyr")         
# install.packages("cleaR")

# Caricamento dei pacchetti:

library(raster)
library(viridis)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(cleaR)

# per pulire tutto l'environment 
rm(list = ls())

#oppure questa funzione del pacchetto cleaR
clear() 


# Setting della working directory:
setwd("C:/lab/exam_teleriv") 

#Importazione e visualizzazione IMMAGINI Sentinel2:  ----
# Le immagini satellitari utilizzate sono state prodotte da Sentinel 2, con risoluzione 60 m x 60 m.

# Layers <- 1 = BLUE, 2 = GREEN, 3 = RED, 4 = SWIR I, 5 = SWIR II, 6 = NIR

## Immagine 2018 ----
# Creazione di una lista cercando nella woking directory elementi in comune
rlist_2018 <- list.files(pattern = "T32TQP_20180926T101021_B")
rlist_2018

# Applico la funzione raster() all'intera lista
import_2018 <- lapply(rlist_2018, raster)
import_2018

# Unione di tutte le bande presenti nella lista in un solo oggetto
img_2018 <- stack(import_2018)

# Visualizzo le informazioni
img_2018

# Plot img_2018
plot(img_2018)  #abbiamo tutte le immagini in un singolo elemento

# Ritaglio l'area di interesse
ext <- c(712500,740000,4844000,4872500) #coordinate per il ritaglio 
#ext <- c(712500, 712800, 4844000, 4844300) #prova
Forest_2018 <- crop(img_2018, ext)


# Plot dell'immagine ritagliata + esportazione in .pdf
pdf("Forest_2018_prova.pdf") 
par(mfrow = c(1,2))  # grafici sono disposti in 1 riga e 2 colonne  
plotRGB(Forest_2018,3,2,1, stretch = "lin") #colori reali
plotRGB(Forest_2018,6,3,2, stretch = "lin") #NIR
dev.off()


## Altre immagini: 2019-2023 ----
# Creazione di una lista cercando nella woking directory elementi in comune
rlist_2019 <- list.files(pattern = "T32TQP_20190911T101021_B")
rlist_2020 <- list.files(pattern = "T32TQP_20200905T101031_B")
rlist_2021 <- list.files(pattern = "T32TQP_20210811T101031_B")
rlist_2022 <- list.files(pattern = "T32TQP_20220806T100611_B")
rlist_2023 <- list.files(pattern = "T32TQP_20230910T100601_B")

# Applico la funzione raster() all'intera lista
import_2019 <- lapply(rlist_2019, raster)
import_2020 <- lapply(rlist_2020, raster)
import_2021 <- lapply(rlist_2021, raster)
import_2022 <- lapply(rlist_2022, raster)
import_2023 <- lapply(rlist_2023, raster)

# Unione di tutte le bande presenti nella lista in un solo oggetto
img_2019 <- stack(import_2019)
img_2020 <- stack(import_2020)
img_2021 <- stack(import_2021)
img_2022 <- stack(import_2022)
img_2023 <- stack(import_2023)

# Ritaglio l'area di interesse
Forest_2019 <- crop(img_2019, ext)
Forest_2020 <- crop(img_2020, ext)
Forest_2021 <- crop(img_2021, ext)
Forest_2022 <- crop(img_2022, ext)
Forest_2023 <- crop(img_2023, ext)

#rm(img_2018,img_2019,img_2020,img_2021,img_2022,img_2023)


# Plot con RGB l'immagine ritagliata

#Forest_2019
pdf("Forest_2019.pdf")
par(mfrow = c(1,2))    
plotRGB(Forest_2019,3,2,1, stretch = "lin") # colori reali 
plotRGB(Forest_2019,6,3,2, stretch = "lin") # NIR
dev.off()

#Forest_2020
pdf("Forest_2020.pdf")
par(mfrow = c(1,2))
plotRGB(Forest_2020,3,2,1, stretch = "lin") # colori reali 
plotRGB(Forest_2020,6,3,2, stretch = "lin") # NIR
dev.off()

#Forest_2021
pdf("Forest_2021.pdf")
par(mfrow = c(1,2))
plotRGB(Forest_2021,3,2,1, stretch = "lin") # colori reali 
plotRGB(Forest_2021,6,3,2, stretch = "lin") # NIR
dev.off()

#Forest_2022
pdf("Forest_2022.pdf")
par(mfrow = c(1,2))
plotRGB(Forest_2022,3,2,1, stretch = "lin") # colori reali 
plotRGB(Forest_2022,6,3,2, stretch = "lin") # NIR
dev.off()

#Forest_2023
pdf("Forest_2023.pdf")
par(mfrow = c(1,2))
plotRGB(Forest_2023,3,2,1, stretch = "lin") # colori reali 
plotRGB(Forest_2023,6,3,2, stretch = "lin") # NIR
dev.off()


## Plot dellle 6 immagini + esportazione in .pdf ----

# Colori reali
pdf("Forest_TC.pdf")
par(mfrow = c(2,3))
plotRGB(Forest_2018,3,2,1, stretch = "lin", main = "2018_TC")
plotRGB(Forest_2019,3,2,1, stretch = "lin", main = "2019_TC")
plotRGB(Forest_2020,3,2,1, stretch = "lin", main = "2020_TC")
plotRGB(Forest_2021,3,2,1, stretch = "lin", main = "2021_TC")
plotRGB(Forest_2022,3,2,1, stretch = "lin", main = "2022_TC")
plotRGB(Forest_2023,3,2,1, stretch = "lin", main = "2023_TC")
dev.off()

# Colori NIR
pdf("Forest_NIR.pdf")
par(mfrow = c(2,3))
plotRGB(Forest_2018,6,3,2, stretch = "lin", main = "2018_NIR")
plotRGB(Forest_2019,6,3,2, stretch = "lin", main = "2019_NIR")
plotRGB(Forest_2020,6,3,2, stretch = "lin", main = "2020_NIR")
plotRGB(Forest_2021,6,3,2, stretch = "lin", main = "2021_NIR")
plotRGB(Forest_2022,6,3,2, stretch = "lin", main = "2022_NIR")
plotRGB(Forest_2023,6,3,2, stretch = "lin", main = "2023_NIR")
dev.off()

#INDICI SPETTRALI ----

## NDVI (NORMALIZED DIFFERENCE VEGETATION INDEX) -----
# NDVI è un indice applicato per quantificare la salute e la densità della vegetazione
# range da -1 a +1

# DVI = NIR - rosso
DVI_2018 <- Forest_2018[[6]] - Forest_2018[[3]]
DVI_2019 <- Forest_2019[[6]] - Forest_2019[[3]]
DVI_2020 <- Forest_2020[[6]] - Forest_2020[[3]]
DVI_2021 <- Forest_2021[[6]] - Forest_2021[[3]]
DVI_2022 <- Forest_2022[[6]] - Forest_2022[[3]]
DVI_2023 <- Forest_2023[[6]] - Forest_2023[[3]]

# NDVI = (NIR - rosso) / (NIR + rosso) = DVI / (NIR + rosso)
NDVI_2018 <- DVI_2018 / (Forest_2018[[6]] + Forest_2018[[3]])
NDVI_2019 <- DVI_2019 / (Forest_2019[[6]] + Forest_2019[[3]])
NDVI_2020 <- DVI_2020 / (Forest_2020[[6]] + Forest_2020[[3]])
NDVI_2021 <- DVI_2021 / (Forest_2021[[6]] + Forest_2021[[3]])
NDVI_2022 <- DVI_2022 / (Forest_2022[[6]] + Forest_2022[[3]])
NDVI_2023 <- DVI_2023 / (Forest_2023[[6]] + Forest_2023[[3]])

# Plot the NDVI
pdf("NDVI.pdf")
par(mfrow = c(2, 3), oma = c(2, 1, 2, 1))    #oma (inf, sx, sup, dx))
cl1 <- colorRampPalette(c("darkblue", "lightgrey", "darkorange"))(100)

plot(NDVI_2018, col = cl1, main = "NDVI_2018")
plot(NDVI_2019, col = cl1, main = "NDVI_2019")
plot(NDVI_2020, col = cl1, main = "NDVI_2020")
plot(NDVI_2021, col = cl1, main = "NDVI_2021")
plot(NDVI_2022, col = cl1, main = "NDVI_2022")
plot(NDVI_2023, col = cl1, main = "NDVI_2023")

dev.off()


# Differenza temporale 2018-2023
NDVI_def <- NDVI_2018 - NDVI_2023   #osservo un  vegetation loss
#NDVI_def_prova <- NDVI_2023 - NDVI_2018    #osservo un vegetation gain 

# Plot NDVI_def che è il risultato dell'analisi temporale 
cl2 <- magma(100)
pdf("NDVI_def.pdf")
plot(NDVI_def, col = cl2, main = "NDVI difference 2018-2023")

# più è alta la differenza, maggiore è la perdita di vegetazione
# se la differenza è negativa significa che c'è un guadagno nella vegetazione
# ci interessa soprattutto la perdita e l'aumento della copertura forestale

# Sembra che in generale l'NDVI non sia diminuito troppo,all'interno dell'area protetta, dal 2018 al 2023.

dev.off()


## NDMI (NORMALIZED DIFFERENCE MOISTURE INDEX) ----
# Indice spettrale utilizzato per rilevare il contenuto di umidità nella vegetazione. 
# Utile per monitorare lo stress idrico nelle piante e può aiutare a valutare la salute 
# generale della vegetazione in una determinata area.

# NDMI = (NIR - SWIR) / (NIR + SWIR)
# SWIR in questo calcolo è la riflettanza nella banda del medio infrarosso
NDMI_2018 <- (Forest_2018[[6]] - Forest_2018[[4]]) / (Forest_2018[[6]] + Forest_2018[[4]])
NDMI_2019 <- (Forest_2019[[6]] - Forest_2019[[4]]) / (Forest_2019[[6]] + Forest_2019[[4]])
NDMI_2020 <- (Forest_2020[[6]] - Forest_2020[[4]]) / (Forest_2020[[6]] + Forest_2020[[4]])
NDMI_2021 <- (Forest_2021[[6]] - Forest_2021[[4]]) / (Forest_2021[[6]] + Forest_2021[[4]])
NDMI_2022 <- (Forest_2022[[6]] - Forest_2022[[4]]) / (Forest_2022[[6]] + Forest_2022[[4]])
NDMI_2023 <- (Forest_2023[[6]] - Forest_2023[[4]]) / (Forest_2023[[6]] + Forest_2023[[4]])

#plot NDMI
pdf("NDMI.pdf")
par(mfrow = c(2,3), oma = c(2, 1, 2, 1))    #oma (inf, sx, sup, dx)
plot(NDMI_2018, col = cl1, main = "NDMI_2018")
plot(NDMI_2019, col = cl1, main = "NDMI_2019")
plot(NDMI_2020, col = cl1, main = "NDMI_2020")
plot(NDMI_2021, col = cl1, main = "NDMI_2021")
plot(NDMI_2022, col = cl1, main = "NDMI_2022")
plot(NDMI_2023, col = cl1, main = "NDMI_2023")
dev.off()

# Calcolo la differenza fra NDMI_2018 e NDMI_2023
NDMI_def <- NDMI_2018 - NDMI_2023
#cl2 <- magma(100)

# Plot di NDMI_def + esportazione in .pdf
pdf("NDMI_def.pdf")
plot(NDMI_def, col = cl2, main = "NDMI difference 2018-2023")
dev.off()

## NBR (NORMALIZED BURN RATIO) ----
# Indice spettrale applicato principalmente per identificare aree bruciate e 
# monitorare la rigenerazione della vegetazione post-incendio.

# NBR = (NIR - SWIR) / (NIR + SWIR)
# SWIR in questo calcolo è la riflettanza nella banda del corto infrarosso
NBR_2018 <- (Forest_2018[[6]] - Forest_2018[[5]]) / (Forest_2018[[6]] + Forest_2018[[5]])
NBR_2019 <- (Forest_2019[[6]] - Forest_2019[[5]]) / (Forest_2019[[6]] + Forest_2019[[5]])
NBR_2020 <- (Forest_2020[[6]] - Forest_2020[[5]]) / (Forest_2020[[6]] + Forest_2020[[5]])
NBR_2021 <- (Forest_2021[[6]] - Forest_2021[[5]]) / (Forest_2021[[6]] + Forest_2021[[5]])
NBR_2022 <- (Forest_2022[[6]] - Forest_2022[[5]]) / (Forest_2022[[6]] + Forest_2022[[5]])
NBR_2023 <- (Forest_2023[[6]] - Forest_2023[[5]]) / (Forest_2023[[6]] + Forest_2023[[5]])

# plot NBR + esportazione in .pdf
pdf("NBR.pdf")
par(mfrow = c(2,3), oma = c(2, 1, 2, 1))
plot(NBR_2018, col = cl1, main = "NBR_2018")
plot(NBR_2019, col = cl1, main = "NBR_2019")
plot(NBR_2020, col = cl1, main = "NBR_2020")
plot(NBR_2021, col = cl1, main = "NBR_2021")
plot(NBR_2022, col = cl1, main = "NBR_2022")
plot(NBR_2023, col = cl1, main = "NBR_2023")
dev.off()

# Calcolo differenza fra NBR_2018 e NBR_2023
NBR_def <- NBR_2018 - NBR_2023

# Plot di NBR_def + esportazione in .pdf
pdf("NBR_def.pdf")
plot(NBR_def, col = cl2, main = "NBR difference 2018-2023")
dev.off()

## MSI (MOISTURE STRESS INDEX) ----
# Indice spettrale che misura lo stress idrico nella vegetazione
# Valuta il contenuto di umidità nelle piante, è utile per monitare siccità e 
# stress idrico nelle aree forestali.

# MSI = SWIR / NIR
# SWIR è la riflettanza nel medio infrarosso
MSI_2018 <- Forest_2018[[4]] / Forest_2018[[6]]
MSI_2019 <- Forest_2019[[4]] / Forest_2019[[6]]
MSI_2020 <- Forest_2020[[4]] / Forest_2020[[6]]
MSI_2021 <- Forest_2021[[4]] / Forest_2021[[6]]
MSI_2022 <- Forest_2022[[4]] / Forest_2022[[6]]
MSI_2023 <- Forest_2023[[4]] / Forest_2023[[6]]

# plot MSI + esportazione in .pdf
pdf("MSI.pdf")
par(mfrow = c(2,3), oma = c(2, 1, 2, 1))
plot(MSI_2018, col = cl1, main = "MSI_2018")
plot(MSI_2019, col = cl1, main = "MSI_2019")
plot(MSI_2020, col = cl1, main = "MSI_2020")
plot(MSI_2021, col = cl1, main = "MSI_2021")
plot(MSI_2022, col = cl1, main = "MSI_2022")
plot(MSI_2023, col = cl1, main = "MSI_2023")
dev.off()

# Calcolo differenza fra MSI_2018 e MSI_2023
MSI_def <- MSI_2018 - MSI_2023

# Plot di MSI_def + esportazione in .pdf
pdf("MSI_def.pdf")
plot(MSI_def, col = cl2, main = "MSI difference 2018-2023")
dev.off()

## Confronto visivo fra i risultati delle analisi temporali dei vari indici
pdf("Confronto_fra_indici.pdf")
par(mfrow = c(2,2))
plot(NDVI_def, col = cl2, main = "NDVI difference 2018-2023")
plot(NDMI_def, col = cl2, main = "NDMI difference 2018-2023")
plot(NBR_def, col = cl2, main = "NBR difference 2018-2023")
plot(MSI_def, col = cl2, main = "MSI difference 2018-2023")

dev.off()


# MULTIVARIATE ANALYSIS -----

# Unisco gli indici NDVI_def e MSAVI_def in un unico oggetto
box <- stack(NDVI_def, NDMI_def, NBR_def, MSI_def)

# Plot
plot(box, col = cl2, main = "Differenze fra NDVI_def, NDMI_def, NBR_def, MSI_def", xaxt = "n", yaxt = "n")

# Effettuo un campionamento casuale di 10000 pixel da box
sr <- sampleRandom(box, 10000)

# Effettuo la PCA (Principal Component Analysis)
PCA <- prcomp(sr)

# Visualizzazione delle informazioni relative alla PCA
summary(PCA)

# Plot della varianza spiegata da ciascuna delle componenti
plot(PCA)

# Proiezione dell'oggetto box nello spazio creato precedentemente usando le CP
PCI <- predict(box, PCA, index = 1:4)
PCI

# Plot della PC1
plot(PCI[[1]], col = cl2,)

# Conversione di PC1 in un dataframe
PC_fin <- as.data.frame(PCI[[1]], xy = T)

# Plot + # Esportazione di  in .pdf
ggplot() + 
  geom_raster(PC_fin, mapping = aes(x = x, y = y, fill = PC1)) + 
  scale_fill_viridis(option = "inferno") +
  theme_bw() #+
#  labs(title = "PC1")

dev.off()

# Il grafico ggplot aiuta a comprendere come la variazione, catturata dalla prima componente principale  
# è distribuita nell'area di studio, basandosi sulle differenze negli indici di input.

# Sembra che la mia PCA indichi che, dove i valori scendono molto, la situazione in termini di copertura vegetale è realmente cambiata. 
# Man mano che i valori aumentano, si prevede che la situazione cambi meno, minor variabilità.



# LAND COVER ----

# Estraggo i valori dalle immagini 
single_nr_2018 <- getValues(Forest_2018)
single_nr_2019 <- getValues(Forest_2019)
single_nr_2020 <- getValues(Forest_2020)
single_nr_2021 <- getValues(Forest_2021)
single_nr_2022 <- getValues(Forest_2022)
single_nr_2023 <- getValues(Forest_2023)

# Classificazione
kcluster_2018 <- kmeans(single_nr_2018, centers=3)
kcluster_2019 <- kmeans(single_nr_2019, centers=3)
kcluster_2020 <- kmeans(single_nr_2020, centers=3)
kcluster_2021 <- kmeans(single_nr_2021, centers=3)
kcluster_2022 <- kmeans(single_nr_2022, centers=3)
kcluster_2023 <- kmeans(single_nr_2023, centers=3)

# Set dei valori
Forest2018_class <- setValues(Forest_2018[[1]], kcluster_2018$cluster)
Forest2019_class <- setValues(Forest_2019[[1]], kcluster_2019$cluster)
Forest2020_class <- setValues(Forest_2020[[1]], kcluster_2020$cluster)
Forest2021_class <- setValues(Forest_2021[[1]], kcluster_2021$cluster)
Forest2022_class <- setValues(Forest_2022[[1]], kcluster_2022$cluster)
Forest2023_class <- setValues(Forest_2023[[1]], kcluster_2023$cluster)

# Plots
pdf("classes.pdf")
cl_class <- colorRampPalette(c('blue','yellow','red'))(3)
par(mfrow = c(2,3),oma = c(0.5, 1, 0.5, 1))
plot(Forest2018_class, col = cl_class, main = "2018")
plot(Forest2019_class, col = cl_class, main = "2019")
plot(Forest2020_class, col = cl_class, main = "2020")
plot(Forest2021_class, col = cl_class, main = "2021")
plot(Forest2022_class, col = cl_class, main = "2022")
plot(Forest2023_class, col = cl_class, main = "2023")

# class 1: suolo nudo, aree urbane, assenza di vegetazione
# class 2: bassa vegetation, aree coltivate con vegetazione
# class 3: foresta, terreno coperto da alberi

dev.off()

# Percentuali 

process_percentage <- function(Forest_class) {
  
  # Calcolo delle frequenze
  freq_values <- freq(Forest_class)
  freq_values
  
  # Calcolo dei pixel totali
  tot <- ncell(Forest_class)
  
  # Calcolo delle percentuali
  perc_values <- round((freq_values * 100) / tot, digit = 5)
  perc_values
  
  return(data.frame(
    RasterName = names(Forest_class),
    Class1_Frequency = freq_values[1,2],
    Class2_Frequency = freq_values[2,2],
    Class3_Frequency = freq_values[3,2],
    Class1_Percentage = perc_values[1,2],
    Class2_Percentage = perc_values[2,2],
    Class3_Percentage = perc_values[3,2]
  ))
}

# Creazione della lista di raster per gli anni desiderati
classes <- list(Forest2018_class, Forest2019_class, Forest2020_class,
                Forest2021_class, Forest2022_class, Forest2023_class)

# Inizializzazione del dataframe per i risultati
results_df <- data.frame(
  RasterName = character(0),
  Class1_Frequency = numeric(0),
  Class2_Frequency = numeric(0),
  Class3_Frequency = numeric(0),
  Class1_Percentage = numeric(0),
  Class2_Percentage = numeric(0),
  Class3_Percentage = numeric(0)
)

# Itera attraverso gli anni
for (i in 1:length(classes)) {
  # Chiama la funzione e ottieni i risultati
  result_values <- process_percentage(classes[[i]])
  
  # Aggiungi i risultati al dataframe complessivo
  results_df <- rbind(results_df, result_values)
}

results_df

# I colori delle classi possono variare da un grafico all'altro.
# Dunque osservare le immagini per assegnare i valori alle classi giuste.
# Quindi, salvarli in un dataframe per eseguire alcuni grafici.


# Creazione del dataframe con i risultati per ogni classe
df_classes <- data.frame(
  year = c(2018, 2019, 2020, 2021, 2022, 2023),
  bare_soil = c(16.05746,18.38842,17.91588,17.35003,18.23397,18.61871),
  low_vegetation = c(34.75661,41.27051, 43.75684,47.50448,49.46035,35.36750),
  forest = c(49.18593,40.34107,38.32728,35.14548,32.30568,46.01379))


# RISULTATI ----
classes_long <- pivot_longer(df_classes, -year, names_to = "class", values_to = "percentage")
classes_long$class <- factor(classes_long$class, levels = c("bare_soil", "low_vegetation", "forest"))

# Grafico a barre
classes_long %>% 
  ggplot(aes(x = factor(year), y = percentage, fill = class))+
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Land cover in Foreste Casentinesi by Year",
       x = "Year",
       y = "Percentage") +
  scale_fill_manual(values = c("lightgoldenrod2", "chartreuse3", "darkgreen"),
                    name = "Class",
                    labels = c("bare soil", "low vegetation", "forest")) +
  theme_bw()+
  theme(legend.position = "bottom")
#ggsave("plot_bar.pdf")

# Grafico a linee 
classes_long %>% 
  ggplot(aes(x = year, y = percentage, color = class)) +
  geom_line(linewidth = 2) +
  labs(title = "Land cover in Foreste Casentinesi by Year",
       x = "Year",
       y = "Percentage") +
  scale_color_manual(values = c("lightgoldenrod2", "chartreuse3", "darkgreen"),
                    name = "Class",
                    labels = c("bare soil", "low vegetation", "forest")) +
  theme_bw()+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size = 11))
#ggsave("plot_line.pdf")
