setwd("C:/Users/Praca_inzynierska/er")
library(DESeq2)
library(biomaRt)
library(edgeR)
library(gplots)
library(ggplot2)
library(rlang)
library(ballgown)
library(ggrepel)
library(tidyverse)
library(ape)
library(vegan)
library(GGally)
library(arrayQualityMetrics)
library(rgl)
library(adegenet)
library(MASS)
library(data.table)
library(plyr)
library(lmtest)
library(reshape2)
library(Rmisc)
library(lmerTest)
library(ComplexHeatmap)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer) #Przydatne do kolorowania
library(vegan)
library(UpSetR)

count_table <- read.csv2("counts_EV.csv", row.names = 1) #wczytanie danych i ustawienie wiersza pierszego jako nagłówek
pheno_data <- data.frame(colnames(count_table), rep(c("CTR","EV"),c(12,8)), c(0,12,24,0,12,24,0,12,24,0,12,24,12,24,12,24,12,24,12,24)) #utworzenie zmiennej pheno_data (12 wierszy CTR i 8 wierszy EV, rozdzielenie czasu w trzeciej kolumnie)
rownames(pheno_data) <- pheno_data[,1] #numerowanie wierszy według danych z kolumny pierwszej
pheno_data <- pheno_data[,-1] #usunięcie kolumny pierwszej
colnames(pheno_data) <- c("condition", "time") #zmiana nazw pozostałych kolumn na condition i time
pheno_data$condition <- as.factor(pheno_data$condition) #konwersja na faktor
pheno_data$time <- as.factor(pheno_data$time) #konwersja na faktor
pheno_data$treatment <- paste0(pheno_data$condition,"_", pheno_data$time) #dodanie kolumny z kombinacją typu i czasu (czy jest to CTR czy EV i czas)
pheno_data$treatment <- as.factor(pheno_data$treatment) #konwersja na faktor


####Utworzenie zmiennych dla próbek 12 i 24h oraz przypisanie im danych####
pheno_data_12 <- pheno_data[pheno_data$time == 12,]
pheno_data_24 <- pheno_data[pheno_data$time == 24,]

####Obrobienie zawartości pheno_data 12 i 24####
pheno_data_12[,2] <- pheno_data_12[,1] #zastąpienie kolumny drugiej pierwszą
pheno_data_12[,1] <- rownames(pheno_data_12)
pheno_data_12 <- pheno_data_12[,-3] #usunięcie kolumny trzeciej
colnames(pheno_data_12) <- c("sample","treatment") #zmiana nazw kolumn

pheno_data_24[,2] <- pheno_data_24[,1] #zastąpienie kolumny drugiej pierwszą
pheno_data_24[,1] <- rownames(pheno_data_24)
pheno_data_24 <- pheno_data_24[,-3]
colnames(pheno_data_24) <- c("sample","treatment")

#####Utworzenie tabel genów dla 12 i 24 godzin####
count_table_12 <- count_table[,colnames(count_table) %in% rownames(pheno_data_12)] 
count_table_24 <- count_table[,colnames(count_table) %in% rownames(pheno_data_24)]

#####Nawiązanie połączenia z bazą danych ensembl. Wybranie danych dot. człowieka i wskazanie konkretnych atrybutów.####
ensemblHs = useMart(host="ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
allgenes.Ensembl = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "gene_biotype", "entrezgene_id", "description", "entrezgene_accession"),mart=ensemblHs) # tabela z wybranymi atrybutami genów
head(allgenes.Ensembl)


#####Przystosowanie próbek. Wygładzenie oraz ujednolicenie.
pheno_data$grouping <- paste(pheno_data$treatment, pheno_data$time, sep = ".")
Group <- factor(paste(pheno_data$time, pheno_data$treatment, sep = "_"))
y <- DGEList(counts = count_table, remove.zeros = TRUE) #usunięcie zer
y <- calcNormFactors(y) #Normalizacja próbek
df_log <- cpm(y, log = TRUE, prior.count = 2) #cpm count per milion
dds.pcoa = pcoa(vegdist(t(df_log <- cpm(y, log = TRUE, prior.count = 2)),
                        method = "euclidean") / 1000)
scores <- dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
cumulative_percent_variance <- (percent / sum( percent)) * 100



################################################################################
##### ANALIZA PCA #####
################################################################################
#Wykonanie PCA (analiza głównych składowych)
df_log_transposed <- t(df_log) #Próbki w wierszach, geny w kolumnach
pca_results <- prcomp(df_log_transposed, scale. = TRUE) 

#Przygotowanie danych do ggplot2
pca_scores <- as.data.frame(pca_results$x)
pca_scores$Sample <- rownames(pca_scores)

#Połączenie wyników PCA z danymi pheno_data
pca_scores <- merge(pca_scores, pheno_data, by.x = "Sample", by.y = "row.names")

#Obliczenie procentu wariancji wyjaśnianej przez każdą składową
percent_variance <- round((pca_results$sdev^2 / sum(pca_results$sdev^2)) * 100, 1)

#Generowanie wykresu PCA
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, 
  color = condition, 
  shape = time)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_shape_manual(values = c("0" = 16, "12" = 17, "24" = 15)) + #Zdefiniowane kształty: Kółko, Trójkąt, Kwadrat
  scale_color_manual(values = c("CTR" = "chartreuse4", "EV" = "darkturquoise")) +
  #Etykiety osi z procentem wyjaśnionej wariancji
  xlab(paste0("PC1 (", percent_variance[1], "%)")) +
  ylab(paste0("PC2 (", percent_variance[2], "%)")) +
  #Etykiety punktów
  geom_text_repel(aes(label = Sample), size = 3, color = "black", max.overlaps = 20) +
  ggtitle("PCA próbek (Normalizowane logCPM)") +
  coord_fixed(ratio = 1) +
  theme_bw() + 
  theme(legend.title = element_text(face = "bold"))

print(pca_plot)

dev.copy(png, "PCA_Analiza_Eksploracyjna.png", width=700, height=500) #zapis do pliku
dev.off() #zamknięcie okna

################################################################################
##### nMDS #####
################################################################################
set.seed(42) #Ustawienie ziarna dla powtarzalności wyników nMDS
nmds_results <- metaMDS(df_log_transposed, 
  distance = "euclidean", 
  k = 2, #2 wymiary (X i Y)
  trymax = 100) #Zwiększenie prób

print(paste("nMDS Stress Value:", nmds_results$stress))

#Przygotowanie danych do ggplot2
#Wyciągnięcie współrzędnych dla próbek (sites)
nmds_scores <- as.data.frame(scores(nmds_results, display = "sites"))
nmds_scores$Sample <- rownames(nmds_scores)

#Połączenie wyników nMDS z danymi fenotypowymi (pheno_data)
nmds_scores <- merge(nmds_scores, pheno_data, by.x = "Sample", by.y = "row.names")

#Generowanie wykresu nMDS
nmds_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  #Elipsy liczone tylko na podstawie grupy condition
  stat_ellipse(aes(group = condition, color = condition), 
    type = "t", linetype = "dashed", alpha = 0.5) +
  #Punkty z zachowaniem kolorów i kształtów
  geom_point(aes(color = condition, shape = time), size = 4, alpha = 0.8) +
  coord_fixed() +
  ggtitle("nMDS próbek (Normalizowane logCPM)") +
  scale_shape_manual(values = c("0" = 16, "12" = 17, "24" = 15)) + 
  scale_color_manual(values = c("CTR" = "chartreuse4", "EV" = "darkturquoise")) +
  geom_text_repel(aes(label = Sample), size = 3, color = "black", max.overlaps = 20) +
  theme_bw()

nmds_plot

dev.copy(png, "nMDS_Analiza_Eksploracyjna.png", width=700, height=500)
dev.off()


nmds_plot_czas <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = time)) +
  geom_point(size = 4) +
  stat_ellipse(aes(fill = time), geom = "polygon", alpha = 0.1) + #Elipsy dla grup czasowych
  coord_fixed() +
  theme_bw() +
  ggtitle("Grupowanie próbek względem czasu (0h, 12h, 24h)")

print(nmds_plot_czas)
dev.copy(png, "nMDS_grupowanie_probek.png", width=700, height=500)
dev.off()


################################################################################
#####klastrowanie hierarchiczne #####
################################################################################

#Obliczenie macierzy odległości
dist_matrix <- dist(df_log_transposed, method = "euclidean")

#Wykonanie klastrowania
hc_results <- hclust(dist_matrix, method = "ward.D2")

#Wyświetlenie dendrogramu
klaster <- plot(hc_results, 
  main = "Hierarchiczne klastrowanie próbek (logCPM)", 
  sub = "",
  xlab = "",
  labels = pheno_data$treatment,
  ylab = "Wysokość (Dystans euklidesowy)")

dev.copy(png, "klastrowanie_hierarchiczne.png", width=700, height=400)
dev.off()


################################################################################
#####Szablon do Vulcano Plot#####
################################################################################
wyglad_volcano = theme(
  text = element_text(family = "serif"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  axis.text.x = element_text(size = 12, face="plain", color = "black"),     #edycja wartości na osi x
  axis.text.y = element_text(size = 12, face="plain", color = "black"),     #edycja wartości na osi y
  axis.title.x = element_text(size = 12, face="bold", color = "black"),
  axis.title.y = element_text(size = 12, face="bold", color = "black"),
  plot.title = element_text(hjust = 0.5, size = 14, face="bold", family = "serif"),
  #legend.position = c(0.95, 0.65),
  #legend.justification = c("right", "top"),   #Punkt zakotwiczenia legendy
  legend.position = "bottom",
  legend.background = element_rect(fill = "white", color = NA), #Tło
  legend.box.background = element_rect(fill = "white", color = NA)) #Obramowanie


                                                                      ####Analiza 1 i 2####
#Porównanie aktywności genów CTR oraz EV - próbka 12h
#Porównuje EV względem CTR
#Czerwone - ekspresja EV wyższa niż w CTR
#Niebieskie - ekspresja CTR wyższa niż w EV

################################################################################
##### ANALIZA 1 #####
##### CTR VS EV 12H #####
################################################################################

#Utworzenie DESeq2 dla 12h
dds_12 <- DESeqDataSetFromMatrix(countData = count_table_12,
  colData = pheno_data_12,
  design = ~ treatment)

#Filtrowanie i analiza
dds_12$treatment <- relevel(dds_12$treatment, ref = "CTR")
dds_12 <- dds_12[rowSums(counts(dds_12)) > 1, ]
dds_12 <- DESeq(dds_12)

#Wyniki porównania EV vs CTR dla 12h
res_12 <- results(dds_12, contrast = c("treatment", "EV", "CTR"))
res_12 <- as.data.frame(res_12)
res_12 <- res_12[!is.na(res_12$padj), ]

#Oznaczenie zmian
res_12$change <- "Bez zmian"
res_12$change <- ifelse(res_12$log2FoldChange > 1 & res_12$padj < 0.05, "Zwiększona ekspresja", res_12$change)
res_12$change <- ifelse(res_12$log2FoldChange < -1 & res_12$padj < 0.05, "Zmniejszona ekspresja", res_12$change)

#Volcano plot dla 12h

CTR_EV_12H_plot <- ggplot(res_12, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
  geom_point(alpha=0.6, size=2) +
  scale_x_continuous(breaks = seq(-8, 8, 2), limits = c(-10, 10)) +
  scale_color_manual(values=c("Zwiększona ekspresja"="red", "Zmniejszona ekspresja"="blue", "Bez zmian"="grey50")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
  xlab("log2 Fold Change") + ylab("-log10(padj)") +
  ggtitle("12h: EV vs CTR - Zmiany w ekspresji genów") +
  theme_minimal() + wyglad_volcano

CTR_EV_12H_plot
dev.copy(png, "Analiza_1_Volcano_CTR_EV_12H.png", width=600, height=400)
dev.off()

#Heatmapa dla analizy 1 #

#Transformacja danych (Variance Stabilizing Transformation) aby dane były porównywalne na heatmapie
vsd_12 <- vst(dds_12, blind=FALSE)

# 2. Wybieramy geny istotne statystycznie (p.adj < 0.05 i |log2FC| > 1)
# Bierzemy nazwy genów z Twojego obiektu res_12 (przed zamianą na numeric!)
sig_genes_names_12 <- rownames(res_12)[res_12$padj < 0.05 & abs(res_12$log2FoldChange) > 1 & !is.na(res_12$padj)]
print(length(sig_genes_names_12)) 

# Opcjonalnie: wyciągnięcie genów (100) o najmniejszym p-value:
#top_genes <- head(order(res_12$padj), 100)
#sig_genes_names <- rownames(res_12)[top_genes]

# 3. Tworzymy macierz tylko dla tych genów
mat_heatmap_12 <- assay(vsd_12)[sig_genes_names_12, ]

# 4. Skalowanie wierszami (Z-score) - żeby pokazać różnice względne (wysoka/niska ekspresja)
mat_scaled_12 <- t(scale(t(mat_heatmap_12)))
mat_heatmap_12 <- mat_heatmap_12[, order(colnames(mat_heatmap_12))] #sortowanie żeby uzyskać kolejność CTR, EV

# Przygotowanie adnotacji (kolorowy pasek na górze informujący czy to CTR czy EV)
annotation_col_analiza_1 <- as.data.frame(colData(dds_12)[, c("treatment")])
rownames(annotation_col_analiza_1) <- colnames(dds_12)
colnames(annotation_col_analiza_1) <- "Grupa"

Analiza_1_heatmapa <- pheatmap(mat_heatmap_12, 
  scale = "row",                # Skalowanie wierszami (kluczowe!)
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  gaps_col = 0,
  annotation_col = annotation_col_analiza_1, # Dodaje pasek z grupami
  annotation_names_col = FALSE,
  show_rownames = FALSE,        # Ukryj nazwy genów jeśli jest ich dużo (zmień na TRUE jeśli mało)
  show_colnames = TRUE, 
  color = colorRampPalette(c("blue", "white", "red"))(100), # Kolory: niebieski(nisko)-czerwony(wysoko)
  border_color = NA,            # Brak ramek wokół kratek
  annotation_colors = list(
    Grupa = c(CTR = "chartreuse", 
      EV = "darkturquoise"),      
    Czas = c(`12` = "lightblue")),
  main = "Heatmapa genów zmienionych (EV vs CTR - 12h)")

Analiza_1_heatmapa
dev.copy(png, "Analiza_1_Heatmapa_EV_CTR_12h.png", width=700, height=500)
dev.off()


################################################################################
##### ANALIZA 2 #####
##### CTR VS EX 24H #####
################################################################################

#Porównanie aktywności genów CTR oraz EV - próbka 24h

#Przygotowanie danych dla 24h
pheno_data_24_filtered <- pheno_data_24
count_table_24_filtered <- count_table_24

#Utworzenie obiektu DESeq2 dla 24h
dds_24 <- DESeqDataSetFromMatrix(countData = count_table_24_filtered,
  colData = pheno_data_24_filtered,
  design = ~ treatment)

#Filtrowanie i analiza
dds_24$treatment <- relevel(dds_24$treatment, ref = "CTR")
dds_24 <- dds_24[rowSums(counts(dds_24)) > 1, ]
dds_24 <- DESeq(dds_24)

#Wyniki porównania EV vs CTR dla 24h
res_24 <- results(dds_24, contrast = c("treatment", "EV", "CTR"))
res_24 <- as.data.frame(res_24)
res_24 <- res_24[!is.na(res_24$padj), ]

#Oznaczenie zmian
res_24$change <- "Bez znaczącej zmiany"
res_24$change <- ifelse(res_24$log2FoldChange > 1 & res_24$padj < 0.05, "Zwiększona ekspresja", res_24$change)
res_24$change <- ifelse(res_24$log2FoldChange < -1 & res_24$padj < 0.05, "Zmniejszona ekspresja", res_24$change)

#Volcano plot dla 24h
CTR_EV_24H_plot <- ggplot(res_24, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
  geom_point(alpha=0.6, size=1.2) +
  scale_x_continuous(breaks = seq(-8, 8, 2), limits = c(-10, 10)) +
  scale_color_manual(values=c("Zwiększona ekspresja"="red", "Zmniejszona ekspresja"="blue", "Bez znaczącej zmiany"="grey50")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
  xlab("log2 Fold Change") + ylab("-log10(padj)") +
  ggtitle("24h: EV vs CTR - Zmiany w ekspresji genów") +
  theme_minimal() + wyglad_volcano

CTR_EV_24H_plot

dev.copy(png, "Analiza_2_Volcano_CTR_EV_24H.png", width=700, height=400)
dev.off()

# Heatmapa dla analizy 2 #

#Transformacja danych (Variance Stabilizing Transformation) aby dane były porównywalne na heatmapie
vsd_24 <- vst(dds_24, blind=FALSE)

# 2. Wybieramy geny istotne statystycznie (p.adj < 0.05 i |log2FC| > 1)
# Bierzemy nazwy genów z Twojego obiektu res_12 (przed zamianą na numeric!)
sig_genes_names_24 <- rownames(res_24)[res_24$padj < 0.05 & abs(res_24$log2FoldChange) > 1 & !is.na(res_24$padj)]
print(length(sig_genes_names_24)) 

# Opcjonalnie: wyciągnięcie genów (100) o najmniejszym p-value:
#top_genes <- head(order(res_12$padj), 100)
#sig_genes_names <- rownames(res_12)[top_genes]

# 3. Tworzymy macierz tylko dla tych genów
mat_heatmap_24 <- assay(vsd_24)[sig_genes_names_24, ]

# 4. Skalowanie wierszami (Z-score) - żeby pokazać różnice względne (wysoka/niska ekspresja)
mat_scaled_24 <- t(scale(t(mat_heatmap_24)))
mat_heatmap_24 <- mat_heatmap_24[, order(colnames(mat_heatmap_24))] #sortowanie żeby uzyskać kolejność CTR, EV


#Przygotowanie adnotacji
annotation_col_analiza_2 <- as.data.frame(colData(dds_24)[, c("treatment")])
rownames(annotation_col_analiza_2) <- colnames(dds_24)
colnames(annotation_col_analiza_2) <- "Grupa"

Analiza_2_heatmapa <- pheatmap(mat_heatmap_24, 
  scale = "row",#Skalowanie wierszami
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  gaps_col = 0,
  annotation_col = annotation_col_analiza_2, #pasek z grupami
  annotation_names_col = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE, 
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color = NA,
  annotation_colors = list(
    Grupa = c(CTR = "chartreuse", 
      EV = "darkturquoise"),      
    Czas = c(`24` = "cadetblue")),
  main = "Heatmapa genów zmienionych (EV vs CTR - 24h)")

dev.copy(png, "Analiza_2_Heatmapa_EV_CTR_24h.png", width=700, height=500)
dev.off()

################################################################################
##### ANALIZA 3 #####
##### DWUCZYNNIKOWA #####
################################################################################
#Analiza dwuczynnikowa EV vs CTR × czas (12h/24h)
#Przygotowanie danych
pheno_data_analiza_3 <- pheno_data[
  (pheno_data$condition %in% c("EV", "CTR")) & 
    (pheno_data$time %in% c(12, 24)), ]
pheno_data_analiza_3 <- droplevels(pheno_data_analiza_3)

count_table_analiza_3 <- count_table[, rownames(pheno_data_analiza_3)]
count_table_analiza_3 <- as.matrix(count_table_analiza_3)
mode(count_table_analiza_3) <- "numeric"

pheno_data_analiza_3$condition <- as.factor(pheno_data_analiza_3$condition)
pheno_data_analiza_3$condition <- relevel(pheno_data_analiza_3$condition, ref = "CTR")
pheno_data_analiza_3$time <- as.factor(pheno_data_analiza_3$time)
pheno_data_analiza_3$time <- relevel(pheno_data_analiza_3$time, ref = "12")
pheno_data_analiza_3 <- droplevels(pheno_data_analiza_3)

#Obiekt DGEList i normalizacja
dge_analiza_3 <- DGEList(counts = count_table_analiza_3, remove.zeros = TRUE)
dge_analiza_3 <- calcNormFactors(dge_analiza_3)

#Macierz modelu dwuczynnikowego
design_analiza_3 <- model.matrix(~ condition * time, data = pheno_data_analiza_3)
print(colnames(design_analiza_3))

#Estymacja dyspersji i fitting GLM
dge_analiza_3 <- estimateDisp(dge_analiza_3, design_analiza_3, robust = TRUE)
fit_analiza_3 <- glmQLFit(dge_analiza_3, design_analiza_3)

#Test efektu głównego condition (EV vs CTR)
test_treatment_analiza_3 <- glmQLFTest(fit_analiza_3, coef = 2)

#Test efektu głównego time (24h vs 12h)  
test_time_analiza_3 <- glmQLFTest(fit_analiza_3, coef = 3)

#Test interakcji condition × time
test_interaction_analiza_3 <- glmQLFTest(fit_analiza_3, coef = 4)

#Wyciągnięcie top genów z interakcji
top_interaction_analiza_3 <- topTags(test_interaction_analiza_3, n = Inf, adjust.method = "BH", p.value = 0.01)
genes_interaction_analiza_3 <- rownames(top_interaction_analiza_3$table)

full_logCPM_analiza_3 <- cpm(dge_analiza_3, log = TRUE)
logCPM_analiza_3_filtered <- full_logCPM_analiza_3[genes_interaction_analiza_3, ]

#Tworzenie kolumny kombinacji Grupa_Czas
pheno_data_analiza_3$comb_sort <- paste0(pheno_data_analiza_3$condition, "_", pheno_data_analiza_3$time, "h")
#Definicja pożądanej kolejności poziomów
set_order <- c("CTR_12h", "EV_12h", "CTR_24h", "EV_24h")
#Zmiana kolejności poziomów na faktorze
pheno_data_analiza_3$comb_sort <- factor(pheno_data_analiza_3$comb_sort, 
  levels = set_order)
#utalenie kolejnosci probek oraz sortownie
kolejnosc_probek <- order(pheno_data_analiza_3$comb_sort, 
  rownames(pheno_data_analiza_3))
#Wyodrębnienie nazw próbek w kolejności
kolejnosc_probek <- rownames(pheno_data_analiza_3)[kolejnosc_probek]
#Sortowanie macierzy ekspresji genów (kolumn)
logCPM_sorted_analiza_3 <- logCPM_analiza_3_filtered[, kolejnosc_probek]
#Sortowanie macierzy adnotacji (wierszy)
annotation_col <- pheno_data_analiza_3[kolejnosc_probek, c("condition", "time")]
colnames(annotation_col) <- c("Grupa", "Czas")

Analiza_3_heatmapa <- pheatmap(logCPM_sorted_analiza_3, 
  scale = "row",              
  cluster_cols = FALSE,
  cluster_rows = TRUE,        
  gaps_col = c(4, 8, 12), #Przerwy po grupach: CTR_12h, EV_12h, CTR_24h
  annotation_col = annotation_col, 
  annotation_names_col = FALSE, 
  show_rownames = FALSE,      
  show_colnames = TRUE,
  annotation_colors = list(
    Grupa = c(CTR = "chartreuse", EV = "darkturquoise"),      
    Czas = c(`12` = "lightblue",`24` = "cadetblue")),
  color = colorRampPalette(c("blue", "white", "red"))(100), 
  main = "Ekspresja genów - analiza dwuczynnikowa - treatment x time",
  name = "Z-score")           

Analiza_3_heatmapa
dev.copy(png, "Analiza_3_heatmapa_dwuczynnikowa.png", width=800, height=600)
dev.off()

    #diagram venna
#A.Efekt Główny Grupy (EV vs CTR)
top_treatment <- topTags(test_treatment_analiza_3, n = Inf, adjust.method = "BH", p.value = 0.01)
genes_treatment <- rownames(top_treatment$table)

#B.Efekt Główny Czasu (24h vs 12h)
top_time <- topTags(test_time_analiza_3, n = Inf, adjust.method = "BH", p.value = 0.01)
genes_time <- rownames(top_time$table)

#C.Efekt Interakcji (Już mamy z poprzedniej analizy)
genes_interaction <- rownames(top_interaction_analiza_3$table) 

#Tworzenie listy zbiorów
lista_dwuczynnikowa <- list(
  "A: Efekt Grupy (EV vs CTR)" = genes_treatment,
  "B: Efekt Czasu (24h vs 12h)" = genes_time,
  "C: Efekt Interakcji (Grupa x Czas)" = genes_interaction)

#Weryfikacja liczebności
print(lapply(lista_dwuczynnikowa, length))
filename_dwuczynnikowy <- "Analiza_3_Diagram_Venna_Dwuczynnikowy.png"

#Usunięcie poprzedniego pliku
if (file.exists(filename_dwuczynnikowy)) file.remove(filename_dwuczynnikowy)

#Generowanie diagramu
venn.plot.dwuczynnikowy <- venn.diagram(
  x = lista_dwuczynnikowa,
  filename = filename_dwuczynnikowy,
  category.names = names(lista_dwuczynnikowa),
  main = "Analiza Dwuczynnikowa: Geny Efektów Głównych vs Interakcji",
  lwd = 2,
  fill = brewer.pal(3, "Pastel1"),
  alpha = 0.50,
  label.col = "black",
  cex = 1.0, 
  cat.cex = 0.8, 
  cat.fontfamily = "sans",
  #Dostosowanie pozycji etykiet
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.06, 0.06, 0.06))

#Podsumowanie wyników
summary(decideTests(test_interaction_analiza_3))


################################################################################
##### ANALIZA 4 #####
##### Wszystkie próbki EV w porównaniu do wszystkich próbek CTR #####
################################################################################

dge_analiza_4 <- DGEList(counts = count_table, remove.zeros = TRUE)
dge_analiza_4 <- calcNormFactors(dge_analiza_4)

model_analiza_4 <- model.matrix(~ treatment, data = pheno_data)
dge_analiza_4 <- estimateDisp(dge_analiza_4, model_analiza_4, robust = TRUE)
fit_analiza_4 <- glmQLFit(dge_analiza_4, model_analiza_4)

test_EV_vs_CTR <- glmQLFTest(fit_analiza_4, coef = 2)
top_EV_vs_CTR  <- topTags(test_EV_vs_CTR, n = Inf, adjust.method = 'BH', p.value = 0.01)
genes_EV_vs_CTR <- rownames(top_EV_vs_CTR$table)

res_EV_vs_CTR <- test_EV_vs_CTR$table
res_EV_vs_CTR$FDR <- p.adjust(res_EV_vs_CTR$PValue, method="BH")
res_EV_vs_CTR$change <- "Bez znaczącej zmiany"
res_EV_vs_CTR$change <- ifelse(res_EV_vs_CTR$logFC > 1 & res_EV_vs_CTR$FDR < 0.05, "Ekspresja zwiększona w EV", res_EV_vs_CTR$change)
res_EV_vs_CTR$change <- ifelse(res_EV_vs_CTR$logFC < -1 & res_EV_vs_CTR$FDR < 0.05, "Ekspresja zwiększona w CTR", res_EV_vs_CTR$change)


Analiza_4_plot <- ggplot(res_EV_vs_CTR, aes(x=logFC, y=-log10(FDR), color=change)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("Ekspresja zwiększona w EV"="red", "Ekspresja zwiększona w CTR"="blue", "Bez znaczącej zmiany"="grey50")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
  xlab("log2 Fold Change") + ylab("-log10(FDR)") +
  ggtitle("EV vs CTR - Zmiany w ekspresji genów (wszystkie godziny)") +
  theme_minimal() + wyglad_volcano

Analiza_4_plot
dev.copy(png, "Analiza_4_volcano_EV_CTR.png", width=700, height=400)
dev.off()

      #heatmapa dla analizy 4#
#Obliczenie pełnej macierzy logCPM (na podstawie obiektu dge_analiza_4)
full_logCPM_analiza_4 <- cpm(dge_analiza_4, log = TRUE)

#Wyfiltrowanie logCPM dla istotnych genów
#genes_EV_vs_CTR to wynik testu glmQLFTest
logCPM_heatmap_analiza_4 <- full_logCPM_analiza_4[genes_EV_vs_CTR, ]

#Przygotowanie adnotacji kolumn (wykorzystanie pełnego pheno_data)
annotation_col_analiza_4 <- pheno_data[, c("condition", "time")]
colnames(annotation_col_analiza_4) <- c("Grupa", "Czas")
annotation_col_analiza_4$Czas <- as.factor(annotation_col_analiza_4$Czas)

sortowanie_analiza_4 <- order(annotation_col_analiza_4$Grupa, 
  annotation_col_analiza_4$Czas)

sortowanie_analiza_4 <- rownames(annotation_col_analiza_4)[sortowanie_analiza_4]

logCPM_sorted_analiza_4 <- logCPM_heatmap_analiza_4[, sortowanie_analiza_4]
annotation_col_sorted_analiza_4 <- annotation_col_analiza_4[sortowanie_analiza_4, ]

Analiza_4_heatmapa <- pheatmap(logCPM_sorted_analiza_4, 
  scale = "row",              
  cluster_cols = FALSE,
  cluster_rows = TRUE,        
  gaps_col = c(4, 8, 12, 16, 20), #Zakładając 4 repliki na punkt czasowy
  annotation_col = annotation_col_sorted_analiza_4, 
  annotation_colors = list(
    Grupa = c(CTR = "chartreuse", EV = "darkturquoise"),
    Czas = c(`0` = "grey", `12` = "lightblue", `24` = "cadetblue")),
  annotation_names_col = FALSE, 
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100), 
  main = "Heatmapa genów zmienionych (EV vs CTR - wszystkie godziny)",
  name = "Z-score")

Analiza_4_heatmapa
dev.copy(png, "Analiza_4_heatmapa_EV_vs_CTR_wszystkie_godziny.png", width=800, height=1000)
dev.off()

#Liczenie genów Up i Down w Analizie 4
summary_A4 <- table(res_EV_vs_CTR$change)
print(summary_A4)


################################################################################
##### ANALIZA 5 #####
##### PRÓBKI CTR WZGLĘDEM CZASU #####
################################################################################

#Filtrowanie próbek CTR
pheno_data_ctr <- pheno_data[pheno_data$condition == "CTR", ]
pheno_data_ctr <- droplevels(pheno_data_ctr)

count_table_ctr <- count_table[, rownames(pheno_data_ctr)]
count_table_ctr <- as.matrix(count_table_ctr)
mode(count_table_ctr) <- "numeric"

dge_ctr <- DGEList(counts = count_table_ctr, remove.zeros = TRUE)
dge_ctr <- calcNormFactors(dge_ctr)

pheno_data_ctr$time <- as.factor(pheno_data_ctr$time)
pheno_data_ctr$time <- relevel(pheno_data_ctr$time, ref = "0")

#Model tylko na czas
design_ctr <- model.matrix(~ time, data = pheno_data_ctr)

dge_ctr <- estimateDisp(dge_ctr, design_ctr, robust = TRUE)
fit_ctr <- glmQLFit(dge_ctr, design_ctr)

#Test efektu czasu
test_time12_ctr <- glmQLFTest(fit_ctr, coef = 2) #Test 12h vs 0h
test_time24_ctr <- glmQLFTest(fit_ctr, coef = 3) #Test 24h vs 0h

top_time12_ctr <- topTags(test_time12_ctr, n = Inf, adjust.method = "BH", p.value = 0.01)
top_time24_ctr <- topTags(test_time24_ctr, n = Inf, adjust.method = "BH", p.value = 0.01)

pheno_ctr <- pheno_data_ctr
pheno_ctr$time <- as.factor(pheno_ctr$time)

dds_ctr <- DESeqDataSetFromMatrix(countData = count_table_ctr,
  colData = pheno_ctr,
  design = ~ time)
dds_ctr$time <- relevel(dds_ctr$time, ref = "0")
dds_ctr <- dds_ctr[rowSums(counts(dds_ctr)) > 1, ]
dds_ctr <- DESeq(dds_ctr)

res12 <- results(dds_ctr, contrast = c("time", "12", "0"))
res24 <- results(dds_ctr, contrast = c("time", "24", "0"))

#Volcano plot dla 12h vs 0h (edgeR)
res_time12 <- as.data.frame(test_time12_ctr$table)
res_time12$FDR <- p.adjust(res_time12$PValue, method="BH")
res_time12$change <- "Bez znaczącej zmiany"
res_time12$change <- ifelse(res_time12$logFC > 1 & res_time12$FDR < 0.05, "Ekspresja zwiększona w 12h", res_time12$change)
res_time12$change <- ifelse(res_time12$logFC < -1 & res_time12$FDR < 0.05, "Ekspresja zmiejszona w 12h", res_time12$change)

#Wykres dla analizy 5
ggplot(res_time12, aes(x=logFC, y=-log10(FDR), color=change)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("Ekspresja zwiększona w 12h"="red", "Ekspresja zmiejszona w 12h"="blue", "Bez znaczącej zmiany"="grey50")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
  xlab("log2 Fold Change") + ylab("-log10(FDR)") +
  ggtitle("CTR: 12h vs 0h - Zmiany w ekspresji genów") +
  theme_minimal() + wyglad_volcano
dev.copy(png, "Analiza_5_volcano_CTR_12h_vs_0h.png", width=600, height=600)
dev.off()

#Volcano plot dla 24h vs 0h (edgeR)
res_time24 <- as.data.frame(test_time24_ctr$table)
res_time24$FDR <- p.adjust(res_time24$PValue, method="BH")
res_time24$change <- "Bez znaczącej zmiany"
res_time24$change <- ifelse(res_time24$logFC > 1 & res_time24$FDR < 0.05, "Ekspresja zwiększona w 24h", res_time24$change)
res_time24$change <- ifelse(res_time24$logFC < -1 & res_time24$FDR < 0.05, "Ekspresja zmiejszona w 24h", res_time24$change)

#Wykres dla analizy 5
ggplot(res_time24, aes(x=logFC, y=-log10(FDR), color=change)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("Ekspresja zwiększona w 24h"="red", "Ekspresja zmiejszona w 24h"="blue", "Bez znaczącej zmiany"="grey50")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
  xlab("log2 Fold Change") + ylab("-log10(FDR)") +
  ggtitle("CTR: 24h vs 0h - Zmiany w ekspresji genów") +
  theme_minimal() + wyglad_volcano
dev.copy(png, "Analiza_5_volcano_CTR_24h_vs_0h.png", width=600, height=600)
dev.off()

#24h vs 12h
contrast_24vs12 <- c(0, -1, 1) 

test_time24vs12_ctr <- glmQLFTest(fit_ctr, contrast = contrast_24vs12)

top_time24vs12_ctr <- topTags(test_time24vs12_ctr, n = Inf, adjust.method = "BH", p.value = 0.01)

#Volcano plot dla 24h vs 12h (edgeR)
res_time24vs12 <- as.data.frame(test_time24vs12_ctr$table)
res_time24vs12$FDR <- p.adjust(res_time24vs12$PValue, method="BH")
res_time24vs12$change <- "Bez znaczącej zmiany"

#Oznaczenia zmian
#logFC > 1: Ekspresja zwiększona w 24h (względem 12h)
res_time24vs12$change <- ifelse(res_time24vs12$logFC > 1 & res_time24vs12$FDR < 0.05, 
  "Ekspresja zwiększona w 24h", res_time24vs12$change)

#logFC < -1: Ekspresja zmniejszona w 24h (względem 12h)
res_time24vs12$change <- ifelse(res_time24vs12$logFC < -1 & res_time24vs12$FDR < 0.05, 
  "Ekspresja zmiejszona w 24h", res_time24vs12$change)

ggplot(res_time24vs12, aes(x=logFC, y=-log10(FDR), color=change)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("Ekspresja zwiększona w 24h"="red", 
    "Ekspresja zmiejszona w 24h"="blue", 
    "Bez znaczącej zmiany"="grey50")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
  xlab("log2 Fold Change") + ylab("-log10(FDR)") +
  ggtitle("CTR: 24h vs 12h - Zmiany w ekspresji genów") +
  theme_minimal() + wyglad_volcano
dev.copy(png, "Analiza_5_volcano_CTR_24h_vs_12h.png", width=600, height=600)
dev.off()

#Test 24h vs 12h w DESeq2
res24vs12_deseq <- results(dds_ctr, contrast = c("time", "24", "12"))

#heatmapy dla analizy 5:
#obrobienie danych; Macierz logCPM dla wszystkich próbek CTR
logCPM_CTR <- cpm(dge_ctr, log = TRUE)

#Adnotacja kolumn (czas)
annotation_ctr <- pheno_data_ctr[, c("time"), drop = FALSE]
colnames(annotation_ctr) <- "Czas"
annotation_ctr$Czas <- as.factor(annotation_ctr$Czas)

#Kolory dla paska czasu
kolory_czas <- list(
  Czas = c(`0` = "grey", 
    `12` = "lightblue", 
    `24` = "darkblue")
)

#12h vs 0h
genes_12vs0 <- rownames(top_time12_ctr$table)
mat_12vs0 <- logCPM_CTR[genes_12vs0, ]

#Sortowanie próbek CTR (0h, 12h, 24h)
kolejnosc_analiza_5_heatmapy <- order(annotation_ctr$Czas)
kolejnosc_analiza_5_heatmapy <- rownames(annotation_ctr)[kolejnosc_analiza_5_heatmapy]

mat_12vs0_sorted <- mat_12vs0[, kolejnosc_analiza_5_heatmapy]
annotation_12vs0_sorted <- annotation_ctr[kolejnosc_analiza_5_heatmapy, , drop = FALSE]

pheatmap(mat_12vs0_sorted, 
  scale = "row", cluster_cols = FALSE, cluster_rows = TRUE, 
  gaps_col = c(4, 8),
  annotation_col = annotation_12vs0_sorted,
  annotation_colors = kolory_czas,
  annotation_names_col = FALSE, 
  show_rownames = FALSE, show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100), 
  name = "Z-Score",
  main = paste0("CTR: 12h vs 0h (Liczba genów: ", length(genes_12vs0), ")"))
dev.copy(png, "Analiza_5_Heatmap_CTR_12vs0.png", width=800, height=1000)
dev.off()

#24h vs 0h
genes_24vs0 <- rownames(top_time24_ctr$table)
mat_24vs0 <- logCPM_CTR[genes_24vs0, ]
mat_24vs0_sorted <- mat_24vs0[, kolejnosc_analiza_5_heatmapy] # Używamy tej samej kolejności próbek
annotation_24vs0_sorted <- annotation_ctr[kolejnosc_analiza_5_heatmapy, , drop = FALSE]

pheatmap(mat_24vs0_sorted, 
  scale = "row", cluster_cols = FALSE, cluster_rows = TRUE, 
  gaps_col = c(4, 8), 
  annotation_col = annotation_24vs0_sorted,
  annotation_colors = kolory_czas,
  annotation_names_col = FALSE, 
  show_rownames = FALSE, show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100), 
  name = "Z-Score",
  main = paste0("CTR: 24h vs 0h (Liczba genów: ", length(genes_24vs0), ")"))
dev.copy(png, "Analiza_5_Heatmap_CTR_24vs0.png", width=800, height=1000)
dev.off()

#24h vs 12h
genes_24vs12 <- rownames(top_time24vs12_ctr$table)
mat_24vs12 <- logCPM_CTR[genes_24vs12, ]
mat_24vs12_sorted <- mat_24vs12[, kolejnosc_analiza_5_heatmapy] # Używamy tej samej kolejności próbek
annotation_24vs12_sorted <- annotation_ctr[kolejnosc_analiza_5_heatmapy, , drop = FALSE]

pheatmap(mat_24vs12_sorted, 
  scale = "row", cluster_cols = FALSE, cluster_rows = TRUE, 
  gaps_col = c(4, 8), 
  annotation_col = annotation_24vs12_sorted,
  annotation_colors = kolory_czas,
  annotation_names_col = FALSE, 
  show_rownames = FALSE, show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100), 
  name = "Z-Score",
  main = paste0("CTR: 24h vs 12h (Liczba genów: ", length(genes_24vs12), ")"))
dev.copy(png, "Analiza_5_Heatmap_CTR_24vs12.png", width=800, height=1000)
dev.off()

#diagram venna
#Wyciągnięcie list genów
lista_genow <- list(
  "12h vs 0h" = rownames(top_time12_ctr$table),
  "24h vs 0h" = rownames(top_time24_ctr$table),
  "24h vs 12h" = rownames(top_time24vs12_ctr$table))

#Sprawdzenie ilosci genów w danym zbiorze
print(lapply(lista_genow, length))

filename <- "Analiza_5_Diagram_Venna_CTR_Czas.png"
#Usunięcie poprzedniego pliku
if (file.exists(filename)) file.remove(filename)

#Generowanie diagramu
venn.plot <- venn.diagram(
  x = lista_genow,
  filename = filename,
  category.names = names(lista_genow),
  main = "Geny istotnie zmienione w czasie (Grupa CTR)",
  #Wygląd kół
  lwd = 2,
  fill = brewer.pal(3, "Pastel1"),
  alpha = 0.50,
  label.col = "black",
  cex = 1, #Wielkość tekstu
  fontfamily = "sans",
  #Wygląd kategorii
  cat.cex = 1, #Wielkość nazw kategorii
  cat.pos = c(-22, 20, 180),
  cat.dist = c(0.06, 0.06, 0.06))

################################################################################
##### WYKRES TRENDU DLA 19 GENÓW #####
################################################################################
interakcja_data <- df_log[genes_interaction_analiza_3[1:19], ] 
interakcja_melted <- melt(interakcja_data)
colnames(interakcja_melted) <- c("Gene", "Sample", "LogCPM")
#Usunięcie nadmiarowych znaków
interakcja_melted$Gene <- gsub("\\|.*", "", interakcja_melted$Gene)
#złączenie tabel
interakcja_melted <- merge(interakcja_melted, pheno_data, by.x="Sample", by.y="row.names")

#Wykres trendu czasowego dla kluczowych 19 genów
wykres_trendu_19_genow <- ggplot(interakcja_melted, aes(x=time, y=LogCPM, color=treatment, group=treatment)) +
  geom_point() +
  facet_wrap(~Gene, scales="free_y") +
  theme_bw() +
  ggtitle("Trend ekspresji genów interakcji (Zależność EV od czasu)") +
  labs(,
    x = "Czas (h)",
    y = "Poziom ekspresji (LogCPM)")

wykres_trendu_19_genow
dev.copy(png, "Wykres_trendu_19_kluczowych_genow.png", width=700, height=500)
dev.off()

################################################################################
##### WYKRES TRENDU DLA TOP 3 #####
################################################################################

#Wybór genów
wybrane_geny <- c("ENSG00000103522", "ENSG00000133110", "ENSG00000196628")
pelne_nazwy_3 <- genes_interaction_analiza_3[grep(paste(wybrane_geny, collapse="|"), genes_interaction_analiza_3)]

#Przygotowanie tabeli
interakcja_data_3 <- df_log[pelne_nazwy_3, ] 
interakcja_melted_3 <- reshape2::melt(interakcja_data_3)
colnames(interakcja_melted_3) <- c("Gene", "Sample", "LogCPM")
interakcja_melted_3$Gene <- gsub("\\|.*", "", interakcja_melted_3$Gene)

#Łączenie z pheno_data
interakcja_melted_3 <- merge(interakcja_melted_3, pheno_data, by.x="Sample", by.y="row.names")

interakcja_melted_3$condition <- as.factor(interakcja_melted_3$condition)
interakcja_melted_3$time <- as.factor(interakcja_melted_3$time)

#Generowanie wykresu
wykres_trendu_top_3_geny <- ggplot(interakcja_melted_3, 
  aes(x = time, y = LogCPM, color = condition, group = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Gene, scales = "free_y") + 
  theme_bw() +
  scale_color_manual(values = c("CTR" = "chartreuse3", "EV" = "darkorchid1")) +
  labs(
    title = "Dynamika ekspresji wybranych genów",
    subtitle = "Porównanie trendów: EV vs CTR (Top 3)",
    x = "Czas (h)",
    y = "Poziom ekspresji (LogCPM)",
    color = "Grupa")

wykres_trendu_top_3_geny

dev.copy(png, "Wykres_trendu_top_3_geny.png", width=700, height=500)
dev.off()

################################################################################
##### Upset plot dla wszystkich analiz #####
################################################################################
#Przygotowanie listy wszystkich analiz
wszystkie_analizy <- list(
  Time_12h = sig_genes_names_12,
  Time_24h = sig_genes_names_24,
  Multifactor = genes_interaction_analiza_3,
  EV4_vs_CTR_global = genes_EV_vs_CTR,
  CTR_24vs12 = genes_24vs12)

geny_A3_i_A4 <- intersect(wszystkie_analizy$Multifactor, wszystkie_analizy$EV4_vs_CTR_global)
cat("Liczba genów interakcji potwierdzona w analizie globalnej (A4):", length(geny_A3_i_A4), "\n")

#UpSet Plot
upset_data <- fromList(wszystkie_analizy)

porownanie_analiz <- upset(upset_data, 
  sets = c("Time_12h", "Time_24h", "Multifactor", "EV4_vs_CTR_global", "CTR_24vs12"),
  main.bar.color = "darkturquoise", 
  sets.bar.color = "chartreuse4",
  order.by = "freq", 
  text.scale = c(1.3, 1.3, 1, 1, 1.5, 1))

porownanie_analiz
dev.copy(png, "porownanie_analiz.png", width=700, height=500)
dev.off()

################################################################################
##### TABELA PODSUMOWUJĄCA: 19 GENÓW INTERAKCJI VS RESZTA ANALIZ #####
################################################################################
#Przygotowanie czystych ID i Symboli
#Bazowa tabela
tabela_A3 <- data.frame(ensembl_gene_id = clean_genes_19, stringsAsFactors = FALSE)

#Symbole i opisy z bazy biomaRt (allgenes.Ensembl)
tabela_A3 <- merge(tabela_A3, 
  allgenes.Ensembl[, c("ensembl_gene_id", "external_gene_name", "description")], 
  by = "ensembl_gene_id", all.x = TRUE)

#Funkcja sprawdzająca obecność genu w danej liście
check_presence <- function(gene_id, gene_list) {
  #Czyszczenie listy porównawczej
  clean_list <- gsub("\\|.*", "", gene_list)
  return(ifelse(gene_id %in% clean_list, "TAK", "-"))}

#Sprawdzanie obecności w pozostałych analizach
tabela_A3$A1_12h <- sapply(tabela_A3$ensembl_gene_id, check_presence, gene_list = sig_genes_names_12)
tabela_A3$A2_24h <- sapply(tabela_A3$ensembl_gene_id, check_presence, gene_list = sig_genes_names_24)
tabela_A3$A4_Globalna <- sapply(tabela_A3$ensembl_gene_id, check_presence, gene_list = genes_EV_vs_CTR)
tabela_A3$A5_Dynamika_CTR <- sapply(tabela_A3$ensembl_gene_id, check_presence, gene_list = genes_24vs12)

print(tabela_A3)

write.csv2(tabela_A3, "Podsumowanie_19_genow_Interakcji.csv", row.names = FALSE)

################################################################################
##### WYKRES TRENDU DLA UNIKALNYCH 2 GENOW #####
################################################################################

#Wybór genów
wybrane_unikalne <- c("ENSG00000108551", "ENSG00000196628")
wybrane_unikalne_pelne_nazwy_3 <- genes_interaction_analiza_3[grep(paste(wybrane_unikalne, collapse="|"), genes_interaction_analiza_3)]

#Przygotowanie tabeli
wybrane_unikalne_interakcja_data_3 <- df_log[wybrane_unikalne_pelne_nazwy_3, ] 
wybrane_unikalne_interakcja_data_3 <- reshape2::melt(wybrane_unikalne_interakcja_data_3)
colnames(wybrane_unikalne_interakcja_data_3) <- c("Gene", "Sample", "LogCPM")
wybrane_unikalne_interakcja_data_3$Gene <- gsub("\\|.*", "", wybrane_unikalne_interakcja_data_3$Gene)

#Łączenie z pheno_data
wybrane_unikalne_interakcja_data_3 <- merge(wybrane_unikalne_interakcja_data_3, pheno_data, by.x="Sample", by.y="row.names")

wybrane_unikalne_interakcja_data_3$condition <- as.factor(wybrane_unikalne_interakcja_data_3$condition)
wybrane_unikalne_interakcja_data_3$time <- as.factor(wybrane_unikalne_interakcja_data_3$time)

#Generowanie wykresu
wykres_trendu_wybrane_unikalne <- ggplot(wybrane_unikalne_interakcja_data_3, 
  aes(x = time, y = LogCPM, color = condition, group = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Gene, scales = "free_y") + 
  theme_bw() +
  scale_color_manual(values = c("CTR" = "chartreuse3", "EV" = "darkorchid1")) +
  labs(
    title = "Dynamika ekspresji wybranych genów",
    subtitle = "Porównanie trendów: EV vs CTR (Unikalne)",
    x = "Czas (h)",
    y = "Poziom ekspresji (LogCPM)",
    color = "Grupa")

wykres_trendu_wybrane_unikalne

dev.copy(png, "wykres_trendu_wybrane_unikalne.png", width=700, height=500)
dev.off()
################################################################################
##### FUNKCJE BIOLOGICZNE (GO) DLA GENÓW UNIKALNYCH (TCF4 i RASD1) #####
################################################################################

#Pobranie szczegółowych danych z biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

unikalne_GO <- getBM(attributes = c('ensembl_gene_id', 
  'external_gene_name', 
  'description',
  'name_1006', #Terminy GO
  'namespace_1003'), #Kategoria (process/function/component)
  filters = 'ensembl_gene_id', 
  values = c("ENSG00000108551", "ENSG00000196628"), 
  mart = ensembl)

#Wyświetlenie wyników w podziale na geny
#TCF4
cat("--- FUNKCJE DLA TCF4 (ENSG00000196628) ---\n")
print(subset(unikalne_GO, external_gene_name == "TCF4" & namespace_1003 == "biological_process")$name_1006[1:5])

#RASD1
cat("\n--- FUNKCJE DLA RASD1 (ENSG00000108551) ---\n")
print(subset(unikalne_GO, external_gene_name == "RASD1" & namespace_1003 == "biological_process")$name_1006[1:5])


################################################################################
##### EKSPORT GENÓW DO CSV #####
################################################################################

#Pobieranie listy genów, dołączenie opisów i zapis do pliku
export_heatmap_genes <- function(gene_list, fileName) {
  clean_ids <- gsub("\\|.*", "", gene_list)
  
  df_to_export <- data.frame(ensembl_gene_id = clean_ids, stringsAsFactors = FALSE)
  #Dołączenie adnotacji z tabeli allgenes.Ensembl
  df_final <- merge(df_to_export, 
    allgenes.Ensembl[, c("ensembl_gene_id", "external_gene_name", "gene_biotype", "description")], 
    by = "ensembl_gene_id", 
    all.x = TRUE)

  write.csv2(df_final, file = fileName, row.names = FALSE)
  cat("Zapisano:", fileName, "- liczba genów:", nrow(df_final), "\n")}

#Analiza 1 (12h)
export_heatmap_genes(sig_genes_names_12, "Geny_Heatmapa_A1_12h.csv")

#Analiza 2 (24h)
export_heatmap_genes(sig_genes_names_24, "Geny_Heatmapa_A2_24h.csv")

#Analiza 3 (Interakcja - 19 genów)
export_heatmap_genes(genes_interaction_analiza_3, "Geny_Heatmapa_A3_Interakcja.csv")

#Analiza 4 (Globalna EV vs CTR)
export_heatmap_genes(genes_EV_vs_CTR, "Geny_Heatmapa_A4_Globalna.csv")

#Analiza 5 (CTR - 3 heatmapy)
export_heatmap_genes(genes_12vs0, "Geny_Heatmapa_A5_12vs0.csv")
export_heatmap_genes(genes_24vs0, "Geny_Heatmapa_A5_24vs0.csv")
export_heatmap_genes(genes_24vs12, "Geny_Heatmapa_A5_24vs12.csv")

#Dostosowanie tabeli do Analizy 3
geny_analiza_3_zapis <- as.data.frame(top_interaction_analiza_3$table)
geny_analiza_3_zapis$ID <- rownames(geny_analiza_3_zapis)
geny_analiza_3_zapis <- geny_analiza_3_zapis[,c(6,1,2,3,4,5)]
rownames(geny_analiza_3_zapis) <- NULL
geny_analiza_3_zapis$ID <- gsub(".*\\|","",geny_analiza_3_zapis$ID)
colnames(geny_analiza_3_zapis)[colnames(geny_analiza_3_zapis) == "logCPM"] <- "AveExpr"
colnames(geny_analiza_3_zapis)[colnames(geny_analiza_3_zapis) == "F"] <- "t"
colnames(geny_analiza_3_zapis)[colnames(geny_analiza_3_zapis) == "PValue"] <- "P.Value"
colnames(geny_analiza_3_zapis)[colnames(geny_analiza_3_zapis) == "FDR"] <- "adj.P.Val"

write.csv2(geny_analiza_3_zapis, "tabela_19_genow.csv", row.names = FALSE)


