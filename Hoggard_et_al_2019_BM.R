# The application of elliptic Fourier analysis in understanding biface shape and symmetry through the British Acheulean
# Authors: Christian Steven Hoggard, John McNabb and James Nathan Cole
# R Script (last updated: 14/01/2019)

# Abstract: Acheulean biface shape and symmetry has fuelled many discussions on past hominin behaviour in regard to the ?meaning? of biface technology. However, few studies have attempted to quantify and investigate their diachronic relationship using a substantial dataset of Acheulean bifaces. Using the British archaeological record as a case study we first perform elliptic Fourier analysis on biface outlines to quantify and better understand the relationship between biface shape and individual interglacial periods. Using the extracted Fourier coefficients we then detail the nature of symmetry throughout this period, before investigating both shape and symmetry in parallel. The importance of size (through biface length) as a factor in biface shape and symmetry is also considered. Results highlight high levels of symmetry from Marine Isotope Stage 13, followed by increasing asymmetry through the British Acheulean. Other observations include a general shift to ?pointed? forms throughout the British Acheulean and the importance of size in high biface symmetry levels. This article concludes by discussing the potential importance of secondary deposition and palimpsest sites in skewering the observed relationships throughout the Palaeolithic.

# OSF link: https://osf.io/td92j/


# ---- Compute-symmetry-and-PCA 
suppressPackageStartupMessages(library(Momocs))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))

MIS_order <-
  c("MIS 13", "MIS 11", "MIS 9", "MIS 7", "MIS 4/3")
Context_order <-
  c(
    "Warren Hill",
    "Boxgrove",
    "Bowman's Lodge",
    "Elveden",
    "Swanscombe",
    "Broom",
    "Furze Platt",
    "Cuxton",
    "Pontnewydd",
    "Lynford"
  )

# reorder according to Marine Isotope Stage (MIS) and Contexts
database_ordered <- 
database %>% 
  mutate(MIS = fct_relevel(MIS, MIS_order)) %>% 
  mutate(Context = fct_relevel(Context, Context_order))

# creation of an outline file with the database supplying metadata
outlinefile <- 
  Out(tpsfile$coo, 
      fac = database_ordered$ID) 

# close all outlines
outlinefile <- coo_close(outlinefile) 
# centre outlines to a common centroid (0,0)
outlinefile <- coo_center(outlinefile) 
# scale outlines to a common centroid size
outlinefile <- coo_scale(outlinefile) 

# creation of EFA class (13 harmonics, we skipped the tests for that); 
# normalisation is suitable given previous procedures
efourierfile <-
  efourier(
    outlinefile,
    nb.h = 13,
    smooth.it = 0,
    norm = TRUE,
    start = FALSE
  )

# extraction of symmetry values
symmetry_df <- 
  symmetry(efourierfile) %>% 
  as.data.frame() %>% 
  rownames_to_column("ID")

# join symmetry values with database
symmetry_df_with_data <- 
  left_join(database_ordered, 
            symmetry_df)

# creation of PCA Class
pca1 <- 
  PCA(efourierfile, 
      scale = FALSE, 
      center = TRUE, 
      fac = database_ordered$ID) 

# creation of scores into a dataframe
scores <- 
  pca1$x %>% 
  as.data.frame() %>% 
  rownames_to_column("ID")
  
# join PCA score values with database
symmetry_df_with_data_and_pca_scores <- 
  left_join(symmetry_df_with_data, 
            scores)

# ---- figure-1

# plot PCA with individual specimens in morphospace, and with confidence ellipses
par(mfrow=c(2,1))
  plot(
    pca1,
    symmetry_df_with_data_and_pca_scores$MIS,
    xax = 1,
    yax = 2,
    points = FALSE,
    center.origin = FALSE,
    zoom = 1,
    grid = TRUE,
    pos.shp = "xy",
    size.shp = 0.4,
    ellipses = FALSE,
    chull = FALSE,
    chull.filled = FALSE,
    eigen = FALSE,
    rug = FALSE,
    title = "Principal Component Analysis (PC1 vs. PC2): XY Warps",
    labelsgroups = FALSE
  )

  plot(
    pca1,
    symmetry_df_with_data_and_pca_scores$MIS,
    xax = 1,
    yax = 2,
    points = FALSE,
    center.origin = FALSE,
    zoom = 1,
    grid = TRUE,
    morphospace = FALSE,
    ellipses = TRUE,
    conf.ellipses = 0.99,
    chull = FALSE,
    chull.filled = FALSE,
    eigen = FALSE,
    rug = FALSE,
    title = "Principal Component Analysis (PC1 vs. PC2): Confidence Ellipses (99%)",
    labelsgroups = TRUE,
    cex.labelsgroups = 1
  )

invisible(dev.off())


# ---- figure-3


# creation of a histogram focusing on symmetry values
figure3a <-
  ggplot(symmetry_df_with_data_and_pca_scores, 
         aes(sym)) + 
  geom_histogram(
    bins = 24,
    binwidth = 0.01,
    colour = "#E69F00",
    fill = "#ffd475"
  ) +
  xlim(0.8, 1) + 
  ylim(0, 100) + 
  geom_density(alpha = 0.2, 
               colour = "darkgrey") + 
  labs(x = "symmetry (AD harmonic coefficients/amplitude)", 
       y = str_glue("count (n={nrow(symmetry_df_with_data_and_pca_scores)})")) + 
  theme_light() + 
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)
  )

# box and whisker of symmetry values categorised by period
figure3b <-
  ggplot(symmetry_df_with_data_and_pca_scores, 
         aes(MIS, sym)) + 
  geom_boxplot(colour = "#E69F00", 
               fill = "#ffd475") + 
  coord_flip() + 
  ylim(0.8, 1) + 
  labs(y = "symmetry (AD harmonic coefficients/amplitude)") + 
  theme_light() +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)
  ) 

# box and whsisker of symmetry values categorised by context
figure3c <-
  ggplot(symmetry_df_with_data_and_pca_scores, 
         aes(Context, sym)) + 
  geom_boxplot(colour = "#E69F00", 
               fill = "#ffd475") + 
  coord_flip() + ylim(0.8, 1) + 
  labs(y = "symmetry (AD harmonic coefficients/amplitude)") + 
  theme_light() +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)
  )


figure3_panel <-
  plot_grid(figure3a,
            figure3b,
            figure3c,
            labels = c('A', 'B', 'C'),
            ncol = 1,
            axis = 'l',
            align = "v")

figure3_panel

# ---- table-3
table_3 <- 
symmetry_df_with_data_and_pca_scores %>%
  group_by(MIS) %>%
  summarise(
    count = n(),
    mean = mean(sym, na.rm = TRUE),
    sd = sd(sym, na.rm = TRUE),
    min = min(sym, na.rm = TRUE),
    max = max(sym, na.rm = TRUE),
    cv = sd(sym, na.rm = TRUE)/mean(sym, na.rm = TRUE)*100)



