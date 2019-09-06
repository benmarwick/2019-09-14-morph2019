# The application of elliptic Fourier analysis in understanding biface shape and symmetry through the British Acheulean
# Authors: Christian Steven Hoggard, John McNabb and James Nathan Cole
# R Script (last updated: 14/01/2019)

# Abstract: Acheulean biface shape and symmetry has fuelled many discussions on past hominin behaviour in regard to the ‘meaning’ of biface technology. However, few studies have attempted to quantify and investigate their diachronic relationship using a substantial dataset of Acheulean bifaces. Using the British archaeological record as a case study we first perform elliptic Fourier analysis on biface outlines to quantify and better understand the relationship between biface shape and individual interglacial periods. Using the extracted Fourier coefficients we then detail the nature of symmetry throughout this period, before investigating both shape and symmetry in parallel. The importance of size (through biface length) as a factor in biface shape and symmetry is also considered. Results highlight high levels of symmetry from Marine Isotope Stage 13, followed by increasing asymmetry through the British Acheulean. Other observations include a general shift to ‘pointed’ forms throughout the British Acheulean and the importance of size in high biface symmetry levels. This article concludes by discussing the potential importance of secondary deposition and palimpsest sites in skewering the observed relationships throughout the Palaeolithic.

# OSF link: https://osf.io/td92j/

# Rstudio version: v.1.1.456

# 1) Set working directory and load packages ------------------------------

setwd() #Change as applicable
if(!require("Momocs")) install.packages('Momocs', repos='http://cran.us.r-project.org') #momocs v.1.2.9
if(!require("tidyverse")) install.packages('tidyverse', repos='http://cran.us.r-project.org') #tidyverse v.1.2.1
if(!require("cowplot")) install.packages('cowplot', repos='http://cran.us.r-project.org') #cowplot v.0.9.3

# 2) Import and screen all data -----------------------------------------------

database <- read.csv("Hoggard_et_al_2019.csv", header = TRUE, row.names = "ID") #.csv table containing metadata
tpsfile <- import_tps("Hoggard_et_al_2019.tps", curves = TRUE) #.tps file containing the outline data for all bifaces
View(database) #view database (observational screening)
print(tpsfile) #print tpsfile (observational screening)

database$MIS = factor(database$MIS, c("MIS 13", "MIS 11", "MIS 9", "MIS 7", "MIS 4/3")) #reorder according to Marine Isotope Stage (MIS)
database$Context = factor(database$Context, c("Warren Hill", "Boxgrove", "Bowman's Lodge", "Elveden", "Swanscombe", "Broom", "Furze Platt", "Cuxton", "Pontnewydd", "Lynford"))
summary(database$MIS) #count data for the different MISs
summary(database$Context) #count data for the different archaeological contexts


# 3) Convert tps file to 'Out' class (see Momocs guide) -------------------

outlinefile <- Out(tpsfile$coo, fac = database) #creation of an outline file with the database supplying metadata


# 4) Standardise specimens prior EFA --------------------------------------

outlinefile <- coo_close(outlinefile) #close all outlines
outlinefile <- coo_center(outlinefile) #centre outlines to a common centroid (0,0)
outlinefile <- coo_scale(outlinefile) #scale outlines to a common centroid size
stack(outlinefile, title="") #stack for visual examination (see ?pile for further information)
panel(outlinefile, fac = "MIS") #visual examination through MIS factor (chronological order)


# 5) Create the EFA class  ------------------------------------------------

calibrate_harmonicpower_efourier(outlinefile) #confirm how many harmonics equate to 99% harmonic power
calibrate_reconstructions_efourier(outlinefile, range=1:20) #confirm through reconstruction (of a random example)
calibrate_deviations_efourier(outlinefile) #confirm through analysis of centroid deviations
efourierfile <- efourier(outlinefile, nb.h = 13, smooth.it = 0, norm = TRUE, start = FALSE) #creation of EFA class (13 harmonics); normalisation is suitable given previous procedures


# 6) Extraction of symmetry values ----------------------------------------

symmetry_1 <- symmetry(efourierfile) #extraction of symmetry values
symmetry_1<-symmetry_1[match(rownames(database),rownames(symmetry_1)),] #order to match database
database<-cbind(database,symmetry_1) #column bind to database
rm(symmetry_1)
View(database)


# 7) Examination of shape -------------------------------------------------

pca1 <- PCA(efourierfile, scale. = FALSE, center = TRUE, fac = Database) #creation of PCA Class
scores <- data.frame(pca1$x) #creation of scores into a dataframe
scores <- scores[match(rownames(database),rownames(scores)),] #match to row ID on database
database <- cbind(database, scores) #column bind to database
rm(scores)
View(database)

scree(pca1) #produces a tibble of the proportion and cumulative percentage for the inertia
PCcontrib(pca1) #produces the main shape changes (PCs)
pdf("Figure_1.pdf", width = 7, height = 7)
plot(pca1, pca1$MIS, xax = 1, yax = 2, points = FALSE, center.origin = FALSE, zoom = 1, grid = TRUE, pos.shp = "xy", size.shp = 0.4, ellipses = FALSE, chull = FALSE, chull.filled = FALSE, eigen = FALSE, rug = FALSE, title = "Principal Component Analysis (PC1 vs. PC2): XY Warps", labelsgroups=FALSE)
plot(pca1, pca1$MIS, xax = 1, yax = 2, points = FALSE, center.origin = FALSE, zoom = 1, grid = TRUE, morphospace = FALSE, ellipses = TRUE, conf.ellipses = 0.99, chull = FALSE, chull.filled = FALSE, eigen = FALSE, rug = FALSE, title = "Principal Component Analysis (PC1 vs. PC2): Confidence Ellipses (99%)", labelsgroups=TRUE, cex.labelsgroups=1)
dev.off()

figure2a <- ggplot(database, aes(MIS, PC1)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.4) +  coord_flip() + labs(x = "Marine Isotope Stage (MIS)", y = "Principal Component 1 (67.29%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) #PC1 scores categorised by MIS
figure2b <- ggplot(database, aes(MIS, PC2)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.4) +  coord_flip() + labs(x = "Marine Isotope Stage (MIS)", y = "Principal Component 2 (11.46%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) #PC1 scores categorised by context
figure2c <- ggplot(database, aes(Context, PC1)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.5) +  coord_flip() + labs(x = "Context", y = "Principal Component 1 (67.29%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) #PC2 scores categorised by MIS
figure2d <- ggplot(database, aes(Context, PC2)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.5) +  coord_flip() + labs(x = "Context", y = "Principal Component 2 (11.46%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) #PC2 scores categorised by context
figure2 <- plot_grid(figure2a, figure2c, figure2b, figure2d, labels= "AUTO", ncol = 2, align = 'v') #synthesis of the four figures
plot(figure2)
ggsave("Figure_2.tiff", plot = last_plot(), dpi = 400, units = "mm", height = 150, width = 250)

MANOVA(pca1, "MIS", test = "Hotelling", retain = 0.99) #MANOVA against MIS
mismanovapw <- MANOVA_PW(pca1, "MIS", retain = 0.99) #pairwise values 
mismanovapw$stars.tab #pairwise values represented by stars
MANOVA(pca1, "Context", test = "Hotelling", retain = 0.99) #MANOVA against context
contextmanovapw <- MANOVA_PW(pca1, "Context", retain = 0.99) #pairwise values
contextmanovapw$stars.tab #pairwise values represented by stars

lda1 <- LDA(pca1, fac = "MIS") #creation of a discriminant analysis by MIS
lda1 #details of the discriminant analysis
plot(lda1, xax = 1, yax = 2, points = TRUE, pch = 20, cex = 0.6, center.origin = FALSE, zoom = 1.8, grid = TRUE, pos.shp = "circle", size.shp = 0.6, ellipses = TRUE, ellipsesax = FALSE, conf.ellipses = 2/3, chull = FALSE, chull.filled = FALSE, eigen = FALSE, rug = FALSE)

# 8) Examination of symmetry ----------------------------------------------

figure3a <- ggplot(database, aes(sym)) + geom_histogram(bins = 24, binwidth = 0.01, colour = "#E69F00", fill = "#ffd475") + xlim(0.8, 1) + ylim(0, 100) + geom_density(alpha=0.2, colour = "darkgrey") + labs(x = "symmetry (AD harmonic coefficients/amplitude)", y = "count (n=)") + theme_light() + theme(text = element_text(size=10), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) #creation of a histogram focusing on symmetry values
figure3b <- ggplot(database, aes(MIS, sym)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475") + coord_flip() + ylim(0.8, 1) + labs(y = "symmetry (AD harmonic coefficients/amplitude)") + theme(text = element_text(size=10), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) #box and whisker of symmetry values categorised by period
figure3c <- ggplot(database, aes(Context, sym)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475") + coord_flip() + ylim(0.8, 1) + labs(y = "symmetry (AD harmonic coefficients/amplitude)") + theme(text = element_text(size=10), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11))#box and whsisker of symmetry values categorised by context
figure3 <- plot_grid(figure3a, figure3b, figure3c, labels = c('A','B','C'), ncol = 1)
plot(figure3)
ggsave("Figure_3.tiff", plot = last_plot(), dpi = 400, units = "mm", height = 250, width = 150)

database %>%
  group_by(MIS) %>%
  summarise(
    count = n(),
    mean = mean(sym, na.rm = TRUE),
    sd = sd(sym, na.rm = TRUE),
    min = min(sym, na.rm = TRUE),
    max = max(sym, na.rm = TRUE),
    cv = sd(sym, na.rm = TRUE)/mean(sym, na.rm = TRUE)*100)

kruskal.test(sym ~ MIS, data = database)
pairwise.wilcox.test(database$sym, database$MIS, p.adj = "bonf") #pairwise Wilcoxon rank sum test for MIS
kruskal.test(sym ~ Context, data = database)
pairwise.wilcox.test(database$sym, database$Context, p.adj = "bonf") #pairwise Wilcoxon rank sum test

figure5a <- ggplot(database, aes(sym, PC1)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Symmetry") + ylab("Principal Component 1 (67.29%)")
figure5b <- ggplot(database, aes(sym, PC2)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Symmetry") + ylab("Principal Component 2 (11.46%)")
figure5c <- ggplot(database, aes(sym, PC3)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Symmetry") + ylab("Principal Component 3 (7.74%)")
figure5d <- ggplot(database, aes(sym, PC4)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Symmetry") + ylab("Principal Component 4 (2.82%)")
figure5e <- ggplot(database, aes(sym, PC5)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Symmetry") + ylab("Principal Component 5 (2.05%)")
figure5 <- plot_grid(figure5a, figure5b, figure5c, figure5d, figure5e, labels = c('A','B','C', 'D', 'E'), nrow=2, ncol=3)
plot(figure5)
ggsave("Figure_5.tiff", plot = last_plot(), dpi = 400, unit = "mm", width = 200, height = 130)

cor(database$PC1, database$sym) #correlation (PC1 vs symmetry)
cor.test(database$PC1, database$sym) #correlation test

cor(database$PC2, database$sym) #correlation (PC2 vs symmetry)
cor.test(database$PC2, database$sym) #correlation test

cor(database$PC3, database$sym) #correlation (PC3 vs symmetry)
cor.test(database$PC3, database$sym) #correlation test

cor(database$PC4, database$sym) #correlation (PC4 vs symmetry)
cor.test(database$PC4, database$sym) #correlation test

cor(database$PC5, database$sym) #correlation (PC5 vs symmetry)
cor.test(database$PC5, database$sym) #correlation test


# 9) Size as a factor--------------------------------------------------

cor(database$Length, database$sym) #correlation (length vs. symmetry)
cor.test(database$Length, database$sym) #correlation test
cor(database$Length, database$PC1) #Correlation (length vs. symmetry)
cor.test(database$Length, database$PC1) #Correlation test
cor(database$Length, database$PC2) #Correlation (length vs. symmetry)
cor.test(database$Length, database$PC2) #Correlation test
cor(database$Length, database$PC3) #Correlation (length vs. symmetry)
cor.test(database$Length, database$PC3) #Correlation test

pairwise.wilcox.test(database$Length, database$MIS, p.adj = "bonf") #pairwise Wilcoxon rank sum test for MIS

figure6a <- ggplot(database, aes(Length, sym)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + geom_smooth(method = "lm", se = FALSE, colour = "grey") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Length (mm)") + ylab("Symmetry")
figure6b <- ggplot(database, aes(Length, PC1)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + geom_smooth(method = "lm", se = FALSE, colour = "grey") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Length (mm)") + ylab("Principal Component 1 (67.29%)")
figure6c <- ggplot(database, aes(Length, PC2)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + geom_smooth(method = "lm", se = FALSE, colour = "grey") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Length (mm)") + ylab("Principal Component 2 (11.46%)")
figure6d <- ggplot(database, aes(Length, PC3)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + geom_smooth(method = "lm", se = FALSE, colour = "grey") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Length (mm)") + ylab("Principal Component 3 (7.74%)")
figure6 <- plot_grid(figure6a, figure6b, figure6c, figure6d, labels = c('A','B','C','D'), nrow=2, ncol=2)
plot(figure6)
ggsave("Figure_6.tiff", plot = last_plot(), dpi = 400, unit = "mm", width = 150, height = 130)
