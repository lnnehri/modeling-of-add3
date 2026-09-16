# 1. Dosyayı oku
file_path <- "/Users/lemannur/Downloads/isoform_summary.tsv"

df <- read.delim(
  file_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Boyut ve kolonları kontrol et
dim(df)
colnames(df)

# İlk satırlar
head(df)
# Veri yapısını kontrol et
str(df)

# Özellikle group/dataset benzeri kolonları görelim
colnames(df)

# Kategorik kolonlardaki değerleri hızlıca incele
for (x in colnames(df)) {
  if (is.character(df[[x]]) || is.factor(df[[x]])) {
    cat("\n---", x, "---\n")
    print(head(unique(df[[x]]), 20))
  }
}
library(dplyr)
library(tidyr)

comp <- df %>%
  filter(source %in% c("ccle_cancer", "gtex_ovary")) %>%
  select(gene, transcript, source,
         tpm_median, frac_median, pct_samples_tpm_zero) %>%
  pivot_wider(
    names_from = source,
    values_from = c(tpm_median, frac_median, pct_samples_tpm_zero)
  ) %>%
  mutate(
    # Isoform kullanımındaki değişim
    delta_frac = frac_median_ccle_cancer - frac_median_gtex_ovary,

    # Expression değişimi
    log2FC_TPM = log2(
      (tpm_median_ccle_cancer + 0.1) /
      (tpm_median_gtex_ovary + 0.1)
    )
  )

# En fazla değişen isoformlar
comp %>%
  arrange(desc(abs(delta_frac))) %>%
  select(gene, transcript, delta_frac, log2FC_TPM) %>%
  head(30)

library(dplyr)
library(tidyr)

comp <- df %>%
  dplyr::filter(source %in% c("ccle_cancer", "gtex_ovary")) %>%
  dplyr::select(
    gene, transcript, source,
    tpm_median, frac_median, pct_samples_tpm_zero
  ) %>%
  tidyr::pivot_wider(
    names_from = source,
    values_from = c(tpm_median, frac_median, pct_samples_tpm_zero)
  ) %>%
  dplyr::mutate(
    delta_frac = frac_median_ccle_cancer - frac_median_gtex_ovary,
    log2FC_TPM = log2(
      (tpm_median_ccle_cancer + 0.1) /
      (tpm_median_gtex_ovary + 0.1)
    )
  )

head(comp)

candidates <- comp %>%
  dplyr::filter(
    abs(delta_frac) >= 0.20,
    tpm_median_ccle_cancer >= 1 | tpm_median_gtex_ovary >= 1
  ) %>%
  dplyr::arrange(desc(abs(delta_frac)))

candidates %>%
  dplyr::select(
    gene, transcript,
    delta_frac, log2FC_TPM,
    tpm_median_ccle_cancer,
    tpm_median_gtex_ovary
  )

nrow(candidates)

switch_genes <- candidates %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(
    any(delta_frac > 0) &
    any(delta_frac < 0)
  ) %>%
  dplyr::arrange(gene, desc(delta_frac)) %>%
  dplyr::ungroup()

switch_genes %>%
  dplyr::select(
    gene, transcript,
    delta_frac, log2FC_TPM,
    tpm_median_ccle_cancer,
    tpm_median_gtex_ovary
  )

switch_candidates <- candidates %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(
    any(delta_frac >= 0.10) &
    any(delta_frac <= -0.10)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(gene, desc(delta_frac))

print(
  switch_candidates %>%
    dplyr::select(
      gene, transcript, delta_frac, log2FC_TPM,
      tpm_median_ccle_cancer,
      tpm_median_gtex_ovary
    ),
  n = Inf
)

# Kaç gen kaldı?
length(unique(switch_candidates$gene))

# Hangi genler?
unique(switch_candidates$gene)
gene_scores <- switch_candidates %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    max_increase = max(delta_frac),
    max_decrease = min(delta_frac),
    switch_strength = max(delta_frac) - min(delta_frac),
    n_isoforms = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(switch_strength))

print(gene_scores, n = Inf)

library(ggplot2)

ggplot(
  switch_candidates,
  aes(
    x = reorder(gene, delta_frac),
    y = delta_frac
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    x = "Gene",
    y = "Δ Isoform Fraction",
    title = "Isoform Switching: CCLE Cancer vs GTEx Ovary"
  ) +
  theme_classic()

top_genes <- c("DIAPH1", "VANGL1", "TMPO", "ACTN1", "ADD3")

top_switches <- switch_candidates %>%
  dplyr::filter(gene %in% top_genes) %>%
  dplyr::arrange(gene, desc(delta_frac))

print(
  top_switches %>%
    dplyr::select(
      gene, transcript,
      delta_frac, log2FC_TPM,
      tpm_median_ccle_cancer,
      tpm_median_gtex_ovary
    ),
  n = Inf
)

# Kaydet
write.csv(
  top_switches,
  "/Users/lemannur/Downloads/top_isoform_switches.csv",
  row.names = FALSE
)

all_switches <- switch_candidates %>%
  dplyr::mutate(
    node = paste0(gene, "_", transcript),
    direction = ifelse(delta_frac > 0, "Cancer_up", "Cancer_down")
  ) %>%
  dplyr::arrange(gene, desc(delta_frac))

# Kontrol
all_switches %>%
  dplyr::select(
    node, gene, transcript,
    direction, delta_frac, log2FC_TPM
  ) %>%
  print(n = Inf)

# Kaydet
write.csv(
  all_switches,
  "/Users/lemannur/Downloads/35_candidate_isoforms.csv",
  row.names = FALSE
)

list.files(
  "/Users/lemannur/Downloads",
  pattern = "\\.(tsv|csv|txt)$",
  full.names = FALSE
)
library(dplyr)
library(pheatmap)

model_data <- all_switches %>%
  dplyr::select(
    node,
    delta_frac,
    log2FC_TPM,
    tpm_median_ccle_cancer,
    tpm_median_gtex_ovary
  ) %>%
  as.data.frame()

rownames(model_data) <- model_data$node
model_data$node <- NULL

# Ölçekle
model_scaled <- scale(model_data)

# Heatmap + hierarchical clustering
pheatmap(
  model_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  main = "Isoform Switching Patterns: CCLE Cancer vs GTEx Ovary"
)

# SD bilgilerini wide formata getir
sd_data <- df %>%
  dplyr::filter(source %in% c("ccle_cancer", "gtex_ovary")) %>%
  dplyr::select(
    gene, transcript, source,
    tpm_sd, frac_sd
  ) %>%
  tidyr::pivot_wider(
    names_from = source,
    values_from = c(tpm_sd, frac_sd)
  )

# 35 isoforma ekle
all_switches2 <- all_switches %>%
  dplyr::left_join(
    sd_data,
    by = c("gene", "transcript")
  )

# Model matrisi
model_data2 <- all_switches2 %>%
  dplyr::select(
    node,
    delta_frac,
    log2FC_TPM,
    frac_sd_ccle_cancer,
    frac_sd_gtex_ovary
  ) %>%
  as.data.frame()

rownames(model_data2) <- model_data2$node
model_data2$node <- NULL

model_scaled2 <- scale(model_data2)

pheatmap::pheatmap(
  model_scaled2,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  main = "Isoform Switch Pattern + Variability"
)

scored_isoforms <- all_switches2 %>%
  dplyr::mutate(
    mean_frac_sd = (frac_sd_ccle_cancer + frac_sd_gtex_ovary) / 2,
    
    switch_score = abs(delta_frac) / (mean_frac_sd + 0.01)
  ) %>%
  dplyr::arrange(desc(switch_score))

scored_isoforms %>%
  dplyr::select(
    gene, transcript,
    delta_frac,
    frac_sd_ccle_cancer,
    frac_sd_gtex_ovary,
    switch_score
  ) %>%
  print(n = Inf)

gene_model_score <- scored_isoforms %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    best_up = max(switch_score[delta_frac > 0], na.rm = TRUE),
    best_down = max(switch_score[delta_frac < 0], na.rm = TRUE),
    combined_switch_score = best_up + best_down,
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(combined_switch_score))

print(gene_model_score, n = Inf)

model_isoforms <- scored_isoforms %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_max(
    order_by = switch_score,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

# Ama yön başına bir tane istiyoruz:
model_isoforms <- scored_isoforms %>%
  dplyr::mutate(
    direction = ifelse(delta_frac > 0, "UP", "DOWN")
  ) %>%
  dplyr::group_by(gene, direction) %>%
  dplyr::slice_max(
    order_by = switch_score,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(gene, direction)

model_isoforms %>%
  dplyr::select(
    gene, transcript, direction,
    delta_frac, switch_score
  ) %>%
  print(n = Inf)

gene_profiles <- model_isoforms %>%
  dplyr::select(
    gene, direction,
    delta_frac, switch_score
  ) %>%
  tidyr::pivot_wider(
    names_from = direction,
    values_from = c(delta_frac, switch_score)
  )

rownames(gene_profiles) <- gene_profiles$gene
gene_profiles$gene <- NULL

# Standardize
gene_scaled <- scale(gene_profiles)

# Genler arasındaki benzerlik
gene_cor <- cor(t(gene_scaled), method = "spearman")

pheatmap::pheatmap(
  gene_cor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Similarity of Isoform-Switch Profiles"
)

pheatmap::pheatmap(
  gene_cor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  main = "Similarity of Isoform-Switch Profiles"
)

library(dplyr)
library(tidyr)
library(pheatmap)

# 1. Her gen için en güçlü UP ve DOWN isoformunu seç
model_isoforms <- scored_isoforms %>%
  dplyr::mutate(
    direction = ifelse(delta_frac > 0, "UP", "DOWN")
  ) %>%
  dplyr::group_by(gene, direction) %>%
  dplyr::slice_max(
    order_by = switch_score,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

# Seçilen isoformları gör
model_isoforms %>%
  dplyr::select(
    gene, transcript, direction,
    delta_frac, switch_score
  ) %>%
  print(n = Inf)


# 2. Gene-level profile oluştur
gene_profiles <- model_isoforms %>%
  dplyr::select(
    gene,
    direction,
    delta_frac,
    switch_score
  ) %>%
  tidyr::pivot_wider(
    names_from = direction,
    values_from = c(delta_frac, switch_score)
  ) %>%
  as.data.frame()


# 3. Gen isimlerini rowname yap
rownames(gene_profiles) <- gene_profiles$gene

gene_profiles$gene <- NULL


# Kontrol
print(rownames(gene_profiles))
print(gene_profiles)


# 4. Scale
gene_scaled <- scale(gene_profiles)


# 5. Genler arası similarity
gene_cor <- cor(
  t(gene_scaled),
  method = "spearman"
)


# İsimlerin gerçekten bulunduğunu kontrol et
print(rownames(gene_cor))
print(colnames(gene_cor))


# 6. Heatmap
pheatmap::pheatmap(
  gene_cor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,

  labels_row = rownames(gene_cor),
  labels_col = colnames(gene_cor),

  show_rownames = TRUE,
  show_colnames = TRUE,

  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,

  border_color = "grey80",

  main = "Similarity of Isoform-Switch Profiles"
)

library(dplyr)
library(tidyr)
library(pheatmap)

# 1. Sample-level veriyi oku
iso_long <- read.delim(
  gzfile("/Users/lemannur/Downloads/isoform_long.tsv.gz"),
  header = TRUE,
  sep = "\t"
)

dim(iso_long)
head(iso_long)
table(iso_long$source)

selected_transcripts <- all_switches$transcript

iso_selected <- iso_long %>%
  dplyr::filter(
    transcript %in% selected_transcripts,
    source %in% c("ccle_cancer", "gtex_ovary")
  )

dim(iso_selected)
length(unique(iso_selected$transcript))
length(unique(iso_selected$sample))

frac_matrix <- iso_selected %>%
  dplyr::select(sample, transcript, fraction) %>%
  tidyr::pivot_wider(
    names_from = transcript,
    values_from = fraction
  ) %>%
  as.data.frame()

rownames(frac_matrix) <- frac_matrix$sample
frac_matrix$sample <- NULL

iso_cor <- cor(
  frac_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)

pheatmap(
  iso_cor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  fontsize_col = 6,
  main = "Sample-level Isoform Co-variation"
)

# ENST -> Gene + direction annotation
annotation_names <- all_switches %>%
  dplyr::select(gene, transcript, delta_frac) %>%
  dplyr::mutate(
    direction = ifelse(delta_frac > 0, "Cancer UP", "Cancer DOWN"),
    label = paste0(gene, "_", transcript, " [", direction, "]")
  )

# Correlation matrix isimlerini değiştir
new_names <- annotation_names$label[
  match(colnames(iso_cor), annotation_names$transcript)
]

colnames(iso_cor) <- new_names
rownames(iso_cor) <- new_names

# Heatmap
pheatmap::pheatmap(
  iso_cor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 5,
  angle_col = 90,
  border_color = "grey80",
  main = "Sample-level Isoform Co-variation"
)

library(dplyr)
library(igraph)

# ADD3 isoform isimlerini bul
add3_nodes <- grep(
  "^ADD3_",
  rownames(iso_cor),
  value = TRUE
)

print(add3_nodes)

# ADD3 ile tüm isoformların korelasyonlarını çıkar
add3_edges <- do.call(
  rbind,
  lapply(add3_nodes, function(add3) {
    
    data.frame(
      from = add3,
      to = colnames(iso_cor),
      rho = iso_cor[add3, ]
    )
    
  })
) %>%
  dplyr::filter(
    from != to,
    abs(rho) >= 0.60
  ) %>%
  dplyr::distinct(from, to, .keep_all = TRUE) %>%
  dplyr::arrange(desc(abs(rho)))

print(add3_edges, n = Inf)

g <- igraph::graph_from_data_frame(
  add3_edges,
  directed = FALSE
)

plot(
  g,
  vertex.size = 8,
  vertex.label.cex = 0.55,
  edge.width = 1 + 3 * abs(E(g)$rho),
  layout = layout_with_fr(g),
  main = "ADD3-centered Isoform Co-variation Network"
)

# Normal data.frame'e çevir
add3_edges <- as.data.frame(add3_edges)

# Görüntüle
print(add3_edges)

# veya RStudio'da tablo olarak
View(add3_edges)

# Kaç bağlantı var?
nrow(add3_edges)

library(dplyr)
library(igraph)

# ============================================================
# 1. ADD3 isoformlarını bul
# ============================================================

add3_nodes <- grep(
  "^ADD3_",
  rownames(iso_cor),
  value = TRUE
)

cat("ADD3 isoforms:\n")
cat(add3_nodes, sep = "\n")


# ============================================================
# 2. Her ADD3 isoformunun diğer isoformlarla korelasyonunu çıkar
# ============================================================

edge_list <- lapply(add3_nodes, function(add3) {

  data.frame(
    from = add3,
    to   = colnames(iso_cor),
    rho  = as.numeric(iso_cor[add3, ]),
    stringsAsFactors = FALSE
  )

})

add3_edges <- do.call(rbind, edge_list)


# ============================================================
# 3. Kendisiyle olan korelasyonu kaldır
#    ve güçlü bağlantıları seç
# ============================================================

add3_edges <- add3_edges %>%
  dplyr::filter(
    from != to,
    !is.na(rho),
    abs(rho) >= 0.60
  ) %>%
  dplyr::arrange(desc(abs(rho)))

add3_edges <- as.data.frame(add3_edges)


# ============================================================
# 4. Pozitif / negatif ilişki ekle
# ============================================================

add3_edges$relationship <- ifelse(
  add3_edges$rho > 0,
  "Positive",
  "Negative"
)


# ============================================================
# 5. Sonuçları kontrol et
# ============================================================

cat("\nNumber of strong ADD3 connections:",
    nrow(add3_edges), "\n\n")

head(add3_edges, 30)

View(add3_edges)


# ============================================================
# 6. Dosyaya kaydet
# ============================================================

write.csv(
  add3_edges,
  "/Users/lemannur/Downloads/ADD3_isoform_network_edges.csv",
  row.names = FALSE
)


# ============================================================
# 7. Network oluştur
# ============================================================

g <- igraph::graph_from_data_frame(
  add3_edges,
  directed = FALSE
)


# Edge özellikleri
E(g)$width <- 1 + 4 * abs(E(g)$rho)

E(g)$lty <- ifelse(
  E(g)$rho > 0,
  1,   # positive = düz
  2    # negative = kesikli
)


# ADD3 node'larını biraz büyük göster
V(g)$size <- ifelse(
  grepl("^ADD3_", V(g)$name),
  12,
  7
)

V(g)$label.cex <- 0.5


# ============================================================
# 8. Network çiz
# ============================================================

set.seed(123)

plot(
  g,
  layout = igraph::layout_with_fr(g),

  vertex.size = V(g)$size,
  vertex.label.cex = V(g)$label.cex,

  edge.width = E(g)$width,
  edge.lty = E(g)$lty,

  main = "ADD3-centered Sample-level Isoform Network"
)

# Çok güçlü bağlantılar
strong_edges <- add3_edges %>%
  dplyr::filter(abs(rho) >= 0.80)

g2 <- igraph::graph_from_data_frame(
  strong_edges,
  directed = FALSE
)

# Kısa label: Gene + ENST son 5 karakter
V(g2)$label <- paste0(
  sub("_.*", "", V(g2)$name),
  "_",
  substr(
    sub(".*ENST", "ENST", V(g2)$name),
    10, 14
  )
)

# ADD3 daha büyük
V(g2)$size <- ifelse(
  grepl("^ADD3_", V(g2)$name),
  16, 8
)

# Edge kalınlığı = korelasyon gücü
E(g2)$width <- 1 + 5 * abs(E(g2)$rho)

# Negatif = kesikli, pozitif = düz
E(g2)$lty <- ifelse(E(g2)$rho < 0, 2, 1)

set.seed(123)

plot(
  g2,
  layout = layout_with_fr(g2),
  vertex.size = V(g2)$size,
  vertex.label.cex = 0.7,
  edge.width = E(g2)$width,
  edge.lty = E(g2)$lty,
  main = "Strong ADD3-centered Isoform Network (|rho| ≥ 0.8)"
)

# En güçlü ilişkileri tablo olarak da gör
strong_edges %>%
  dplyr::arrange(desc(abs(rho))) %>%
  head(20)

strong_edges <- add3_edges %>%
  dplyr::filter(abs(rho) >= 0.70)

strong_edges %>%
  dplyr::arrange(desc(abs(rho)))

g3 <- igraph::graph_from_data_frame(
  strong_edges,
  directed = FALSE
)

V(g3)$label <- sub(
  "_ENST00000",
  "_",
  sub(" \\[.*", "", V(g3)$name)
)

V(g3)$size <- ifelse(
  grepl("^ADD3_", V(g3)$name),
  16, 8
)

E(g3)$width <- 1 + 4 * abs(E(g3)$rho)
E(g3)$lty <- ifelse(E(g3)$rho < 0, 2, 1)

set.seed(123)

plot(
  g3,
  layout = layout_with_fr(g3),
  vertex.size = V(g3)$size,
  vertex.label.cex = 0.7,
  edge.width = E(g3)$width,
  edge.lty = E(g3)$lty,
  main = "ADD3-centered Isoform Network (|rho| ≥ 0.7)"
)

library(dplyr)
library(tidyr)
library(pheatmap)

# Fonksiyon: source için correlation matrix oluştur
make_cor <- function(data, source_name) {
  
  tmp <- data %>%
    filter(source == source_name) %>%
    select(sample, transcript, fraction) %>%
    pivot_wider(
      names_from = transcript,
      values_from = fraction
    ) %>%
    as.data.frame()
  
  rownames(tmp) <- tmp$sample
  tmp$sample <- NULL
  
  cor(
    tmp,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
}

# CCLE cancer
cor_cancer <- make_cor(
  iso_selected,
  "ccle_cancer"
)

# GTEx normal ovary
cor_normal <- make_cor(
  iso_selected,
  "gtex_ovary"
)

# Aynı isoform sırasına getir
common <- intersect(
  colnames(cor_cancer),
  colnames(cor_normal)
)

cor_cancer <- cor_cancer[common, common]
cor_normal <- cor_normal[common, common]

# Korelasyon değişimi
delta_cor <- cor_cancer - cor_normal

# Gene + transcript isimleri
labels <- annotation_names$label[
  match(common, annotation_names$transcript)
]

rownames(delta_cor) <- labels
colnames(delta_cor) <- labels

# Heatmap
pheatmap(
  delta_cor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 5,
  fontsize_col = 5,
  angle_col = 90,
  main = "Change in Isoform Co-variation: Cancer - Normal"
)

make_cor <- function(data, source_name) {
  
  tmp <- data %>%
    dplyr::filter(source == source_name) %>%
    dplyr::select(sample, transcript, fraction) %>%
    tidyr::pivot_wider(
      names_from = transcript,
      values_from = fraction
    )
  
  tmp <- as.data.frame(tmp)
  
  rownames(tmp) <- tmp$sample
  tmp$sample <- NULL
  
  stats::cor(
    tmp,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
}

# CCLE
cor_cancer <- make_cor(
  iso_selected,
  "ccle_cancer"
)

# GTEx
cor_normal <- make_cor(
  iso_selected,
  "gtex_ovary"
)

dim(cor_cancer)
dim(cor_normal)

# Ortak isoformlar
common <- intersect(
  colnames(cor_cancer),
  colnames(cor_normal)
)

cor_cancer2 <- cor_cancer[common, common]
cor_normal2 <- cor_normal[common, common]

# Cancer - Normal korelasyon farkı
delta_cor <- cor_cancer2 - cor_normal2


# ADD3 isoformlarını bul
add3_transcripts <- all_switches %>%
  dplyr::filter(gene == "ADD3") %>%
  dplyr::pull(transcript)

add3_transcripts <- intersect(
  add3_transcripts,
  rownames(delta_cor)
)


# ADD3 bağlantılarını çıkar
add3_delta <- do.call(
  rbind,
  lapply(add3_transcripts, function(x) {
    
    data.frame(
      ADD3_isoform = x,
      partner = colnames(delta_cor),
      rho_cancer = cor_cancer2[x, ],
      rho_normal = cor_normal2[x, ],
      delta_rho = delta_cor[x, ],
      stringsAsFactors = FALSE
    )
    
  })
)

# Kendisiyle bağlantıyı çıkar
add3_delta <- add3_delta %>%
  dplyr::filter(
    ADD3_isoform != partner
  ) %>%
  dplyr::arrange(desc(abs(delta_rho)))

head(add3_delta, 30)

changed_add3 <- add3_delta %>%
  dplyr::filter(abs(delta_rho) >= 0.5) %>%
  dplyr::mutate(
    change_type = dplyr::case_when(
      rho_normal < 0 & rho_cancer > 0 ~ "Negative → Positive",
      rho_normal > 0 & rho_cancer < 0 ~ "Positive → Negative",
      abs(rho_cancer) > abs(rho_normal) ~ "Strengthened in Cancer",
      TRUE ~ "Weakened in Cancer"
    )
  ) %>%
  dplyr::arrange(desc(abs(delta_rho)))

changed_add3 %>%
  dplyr::select(
    ADD3_isoform, partner,
    rho_normal, rho_cancer,
    delta_rho, change_type
  ) %>%
  head(30)

library(igraph)

# En güçlü değişen bağlantılar
net_edges <- changed_add3 %>%
  dplyr::filter(abs(delta_rho) >= 0.5) %>%
  dplyr::select(
    from = ADD3_isoform,
    to = partner,
    delta_rho,
    rho_normal,
    rho_cancer
  )

g_change <- igraph::graph_from_data_frame(
  net_edges,
  directed = FALSE
)

# ADD3 node'ları büyük
V(g_change)$size <- ifelse(
  V(g_change)$name %in% add3_transcripts,
  15, 7
)

# Edge kalınlığı = değişimin büyüklüğü
E(g_change)$width <- 1 + 4 * abs(E(g_change)$delta_rho)

# Pozitif Δrho düz, negatif kesikli
E(g_change)$lty <- ifelse(
  E(g_change)$delta_rho > 0,
  1, 2
)

set.seed(123)

plot(
  g_change,
  layout = layout_with_fr(g_change),
  vertex.size = V(g_change)$size,
  vertex.label.cex = 0.55,
  edge.width = E(g_change)$width,
  edge.lty = E(g_change)$lty,
  main = "ADD3 Isoform Network Rewiring: Normal → Cancer"
)

library(dplyr)
library(ggplot2)

focus_ADD3 <- iso_long %>%
  dplyr::filter(
    transcript %in% c(
      "ENST00000356080",  # uzun
      "ENST00000360162"   # kısa
    ),
    source %in% c("gtex_ovary", "ccle_cancer")
  ) %>%
  dplyr::mutate(
    isoform = dplyr::recode(
      transcript,
      "ENST00000356080" = "Long (ENST00000356080)",
      "ENST00000360162" = "Short (ENST00000360162)"
    )
  )

# Sayısal özet
focus_ADD3 %>%
  dplyr::group_by(isoform, source) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_fraction = median(fraction, na.rm = TRUE),
    mean_fraction = mean(fraction, na.rm = TRUE),
    sd_fraction = sd(fraction, na.rm = TRUE),
    .groups = "drop"
  )

# Görselleştir
ggplot(
  focus_ADD3,
  aes(x = source, y = fraction, fill = source)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  facet_wrap(~isoform) +
  theme_classic() +
  labs(
    x = NULL,
    y = "ADD3 isoform fraction",
    title = "ADD3 Long vs Short Isoform: Normal → Cancer"
  ) +
  theme(legend.position = "none")

add3_score <- iso_long %>%
  dplyr::filter(
    transcript %in% c(
      "ENST00000356080",
      "ENST00000360162"
    ),
    source %in% c("gtex_ovary", "ccle_cancer")
  ) %>%
  dplyr::select(
    sample, source, transcript, fraction
  ) %>%
  tidyr::pivot_wider(
    names_from = transcript,
    values_from = fraction
  ) %>%
  dplyr::mutate(
    ADD3_transition_score =
      ENST00000356080 - ENST00000360162
  )

ggplot(
  add3_score,
  aes(x = source, y = ADD3_transition_score)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  theme_classic() +
  labs(
    x = NULL,
    y = "ADD3 transition score (Long − Short)",
    title = "ADD3 Isoform Transition Score"
  )

# Grupların score özeti
add3_score %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_score = median(ADD3_transition_score, na.rm = TRUE),
    mean_score = mean(ADD3_transition_score, na.rm = TRUE),
    sd_score = sd(ADD3_transition_score, na.rm = TRUE)
  )

# Normal vs Cancer testi
wilcox.test(
  ADD3_transition_score ~ source,
  data = add3_score,
  exact = FALSE
)

add3_ordered <- add3_score %>%
  dplyr::arrange(ADD3_transition_score) %>%
  dplyr::mutate(
    order = dplyr::row_number()
  )

ggplot(
  add3_ordered,
  aes(
    x = order,
    y = ADD3_transition_score,
    color = source
  )
) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    x = "Samples ordered by ADD3 transition score",
    y = "ADD3 Long − Short score",
    title = "Normal-like → Cancer-like ADD3 Transition"
  )
# ADD3 score'u tüm seçilmiş isoformlarla birleştir
transition_data <- iso_selected %>%
  dplyr::select(sample, transcript, fraction) %>%
  tidyr::pivot_wider(
    names_from = transcript,
    values_from = fraction
  ) %>%
  dplyr::left_join(
    add3_score %>%
      dplyr::select(sample, ADD3_transition_score),
    by = "sample"
  )

# Her isoformun ADD3 transition score ile korelasyonu
isoform_transition_cor <- sapply(
  selected_transcripts,
  function(x) {
    cor(
      transition_data[[x]],
      transition_data$ADD3_transition_score,
      method = "spearman",
      use = "complete.obs"
    )
  }
)

transition_results <- data.frame(
  transcript = names(isoform_transition_cor),
  rho = as.numeric(isoform_transition_cor)
) %>%
  dplyr::left_join(
    all_switches %>% dplyr::select(gene, transcript),
    by = "transcript"
  ) %>%
  dplyr::arrange(desc(abs(rho)))

transition_results

# Sadece CCLE cancer sample'ları
cancer_samples <- add3_score %>%
  dplyr::filter(source == "ccle_cancer") %>%
  dplyr::pull(sample)

cancer_transition <- transition_data %>%
  dplyr::filter(sample %in% cancer_samples)

# Cancer içinde ADD3 score ile korelasyon
cancer_cor <- sapply(
  selected_transcripts,
  function(x) {
    cor(
      cancer_transition[[x]],
      cancer_transition$ADD3_transition_score,
      method = "spearman",
      use = "complete.obs"
    )
  }
)

cancer_results <- data.frame(
  transcript = names(cancer_cor),
  rho_cancer = as.numeric(cancer_cor)
) %>%
  dplyr::left_join(
    transition_results %>%
      dplyr::select(transcript, gene, rho),
    by = "transcript"
  ) %>%
  dplyr::rename(rho_all = rho) %>%
  dplyr::arrange(desc(abs(rho_cancer)))

cancer_results

stable_module <- cancer_results %>%
  dplyr::filter(
    !is.na(rho_all),
    !is.na(rho_cancer),
    sign(rho_all) == sign(rho_cancer),
    abs(rho_all) >= 0.60,
    abs(rho_cancer) >= 0.40
  ) %>%
  dplyr::arrange(desc(abs(rho_cancer)))

stable_module

stable_module %>%
  dplyr::select(
    gene,
    transcript,
    rho_all,
    rho_cancer
  )

non_add3 <- cancer_results %>%
  dplyr::filter(gene != "ADD3") %>%
  dplyr::arrange(desc(abs(rho_cancer)))

non_add3 %>%
  dplyr::select(
    gene,
    transcript,
    rho_all,
    rho_cancer
  ) %>%
  head(20)

# Normal sample'lar
normal_samples <- add3_score %>%
  dplyr::filter(source == "gtex_ovary") %>%
  dplyr::pull(sample)

normal_transition <- transition_data %>%
  dplyr::filter(sample %in% normal_samples)

# Normal içinde ADD3 score korelasyonları
normal_cor <- sapply(
  selected_transcripts,
  function(x) {
    cor(
      normal_transition[[x]],
      normal_transition$ADD3_transition_score,
      method = "spearman",
      use = "complete.obs"
    )
  }
)

# Normal vs Cancer karşılaştır
rewiring <- cancer_results %>%
  dplyr::mutate(
    rho_normal = normal_cor[transcript],
    delta_rho = rho_cancer - rho_normal
  ) %>%
  dplyr::filter(gene != "ADD3") %>%
  dplyr::arrange(desc(abs(delta_rho)))

rewiring %>%
  dplyr::select(
    gene, transcript,
    rho_normal,
    rho_cancer,
    delta_rho
  ) %>%
  head(20)

library(dplyr)
library(igraph)

rewiring_edges <- rewiring %>%
  dplyr::filter(abs(delta_rho) >= 0.25) %>%
  dplyr::mutate(
    from = "ADD3_switch",
    to = paste0(gene, "_", transcript)
  ) %>%
  dplyr::select(
    from, to,
    rho_normal,
    rho_cancer,
    delta_rho
  )

g_rewire <- igraph::graph_from_data_frame(
  rewiring_edges,
  directed = FALSE
)

V(g_rewire)$size <- ifelse(
  V(g_rewire)$name == "ADD3_switch",
  18, 8
)

E(g_rewire)$width <- 2 + 6 * abs(E(g_rewire)$delta_rho)

# düz = ilişki artıyor
# kesikli = ilişki azalıyor
E(g_rewire)$lty <- ifelse(
  E(g_rewire)$delta_rho > 0,
  1, 2
)

set.seed(123)

plot(
  g_rewire,
  layout = layout_with_fr(g_rewire),
  vertex.size = V(g_rewire)$size,
  vertex.label.cex = 0.65,
  edge.width = E(g_rewire)$width,
  edge.lty = E(g_rewire)$lty,
  main = "ADD3-centered Isoform Rewiring: Normal → Cancer"
)

set.seed(123)

bootstrap_delta <- function(transcript, B = 2000) {
  
  normal_x <- normal_transition[[transcript]]
  normal_y <- normal_transition$ADD3_transition_score
  
  cancer_x <- cancer_transition[[transcript]]
  cancer_y <- cancer_transition$ADD3_transition_score
  
  boot_delta <- replicate(B, {
    
    # Sample'ları grup içinde yeniden örnekle
    idx_n <- sample(seq_along(normal_x), replace = TRUE)
    idx_c <- sample(seq_along(cancer_x), replace = TRUE)
    
    rho_n <- suppressWarnings(
      cor(
        normal_x[idx_n],
        normal_y[idx_n],
        method = "spearman",
        use = "complete.obs"
      )
    )
    
    rho_c <- suppressWarnings(
      cor(
        cancer_x[idx_c],
        cancer_y[idx_c],
        method = "spearman",
        use = "complete.obs"
      )
    )
    
    rho_c - rho_n
  })
  
  c(
    CI_low = quantile(boot_delta, 0.025, na.rm = TRUE),
    CI_high = quantile(boot_delta, 0.975, na.rm = TRUE)
  )
}


# Şu an network'te bulunan isoformlar
test_isoforms <- rewiring %>%
  dplyr::filter(
    gene != "ADD3",
    abs(delta_rho) >= 0.25
  ) %>%
  dplyr::pull(transcript)


boot_results <- lapply(
  test_isoforms,
  bootstrap_delta
)

boot_results <- do.call(rbind, boot_results)

boot_results <- data.frame(
  transcript = test_isoforms,
  CI_low = boot_results[,1],
  CI_high = boot_results[,2]
)


# Ana tabloyla birleştir
rewiring_test <- rewiring %>%
  dplyr::inner_join(
    boot_results,
    by = "transcript"
  ) %>%
  dplyr::mutate(
    significant = CI_low > 0 | CI_high < 0
  ) %>%
  dplyr::select(
    gene,
    transcript,
    rho_normal,
    rho_cancer,
    delta_rho,
    CI_low,
    CI_high,
    significant
  ) %>%
  dplyr::arrange(desc(abs(delta_rho)))

rewiring_test

rewiring_test

sig_edges <- rewiring_test %>%
  dplyr::filter(significant == TRUE) %>%
  dplyr::mutate(
    from = "ADD3_switch",
    to = paste0(gene, "_", transcript)
  )

sig_edges

model_df <- transition_data %>%
  dplyr::select(
    sample,
    ADD3_transition_score,
    ENST00000553882
  ) %>%
  dplyr::left_join(
    add3_score %>% dplyr::select(sample, source),
    by = "sample"
  )

ggplot(
  model_df,
  aes(
    x = ADD3_transition_score,
    y = ENST00000553882,
    color = source
  )
) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  theme_classic() +
  labs(
    x = "ADD3 Long − Short transition score",
    y = "ACTN1 ENST00000553882 fraction",
    title = "Loss of ADD3–ACTN1 Coupling in Cancer"
  )

# Model için iki değişken
cluster_df <- model_df %>%
  dplyr::select(
    sample,
    source,
    ADD3_transition_score,
    ENST00000553882
  ) %>%
  tidyr::drop_na()

# İki değişkeni standardize et
X <- scale(
  cluster_df[, c(
    "ADD3_transition_score",
    "ENST00000553882"
  )]
)

# 3 latent-benzeri grup
set.seed(123)

km3 <- kmeans(
  X,
  centers = 3,
  nstart = 100
)

cluster_df$state <- factor(km3$cluster)

# Görselleştir
library(ggplot2)

ggplot(
  cluster_df,
  aes(
    x = ADD3_transition_score,
    y = ENST00000553882,
    color = state,
    shape = source
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  labs(
    x = "ADD3 Long − Short score",
    y = "ACTN1 ENST00000553882 fraction",
    color = "State",
    shape = "Dataset",
    title = "Data-driven ADD3–ACTN1 States"
  )

table(
  cluster_df$state,
  cluster_df$source
)

cluster_df %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(
    n = dplyr::n(),
    ADD3 = mean(ADD3_transition_score),
    ACTN1 = mean(ENST00000553882)
  )

table(
  cluster_df$state,
  cluster_df$source
)

library(ggplot2)

# State isimlerini biyolojik olarak daha anlaşılır yap
cluster_df$state_label <- dplyr::recode(
  cluster_df$state,
  "3" = "Normal-like",
  "2" = "Cancer: Low ADD3 switch",
  "1" = "Cancer: High ADD3 switch"
)

ggplot(
  cluster_df,
  aes(
    x = ADD3_transition_score,
    y = ENST00000553882,
    color = state_label,
    shape = source
  )
) +
  geom_point(size = 3.5, alpha = 0.8) +
  theme_classic() +
  labs(
    x = "ADD3 Long − Short score",
    y = "ACTN1 ENST00000553882 fraction",
    color = "Data-driven state",
    shape = "Source",
    title = "ADD3–ACTN1 State Separation"
  )

state_counts <- cluster_df %>%
  dplyr::count(state_label, source)

ggplot(
  state_counts,
  aes(
    x = state_label,
    y = n,
    fill = source
  )
) +
  geom_col() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of samples",
    fill = "Source",
    title = "Composition of ADD3–ACTN1 States"
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# Cancer sample + state bilgisi
cancer_states <- cluster_df %>%
  dplyr::filter(source == "ccle_cancer") %>%
  dplyr::select(sample, state_label)

# Isoform fraction'ları state ile birleştir
cancer_isoforms <- iso_selected %>%
  dplyr::filter(
    source == "ccle_cancer",
    !gene %in% c("ADD3", "ACTN1")
  ) %>%
  dplyr::inner_join(
    cancer_states,
    by = "sample"
  )

# Low vs High ADD3 switch karşılaştırması
state_comparison <- cancer_isoforms %>%
  dplyr::group_by(gene, transcript) %>%
  dplyr::summarise(
    
    median_high = median(
      fraction[state_label == "Cancer: High ADD3 switch"],
      na.rm = TRUE
    ),
    
    median_low = median(
      fraction[state_label == "Cancer: Low ADD3 switch"],
      na.rm = TRUE
    ),
    
    delta_fraction = median_high - median_low,
    
    p_value = wilcox.test(
      fraction ~ state_label,
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(FDR, desc(abs(delta_fraction)))

state_comparison

state_comparison %>%
  dplyr::filter(
    abs(delta_fraction) >= 0.05
  ) %>%
  dplyr::select(
    gene,
    transcript,
    median_low,
    median_high,
    delta_fraction,
    p_value,
    FDR
  )

sig_state_isoforms <- state_comparison %>%
  dplyr::filter(
    FDR < 0.05,
    abs(delta_fraction) >= 0.05
  ) %>%
  dplyr::arrange(desc(abs(delta_fraction)))

sig_state_isoforms

sig_state_isoforms %>%
  dplyr::select(
    gene, transcript,
    median_low, median_high,
    delta_fraction,
    p_value, FDR
  ) %>%
  print(width = Inf)

SPTAN1_hit <- sig_state_isoforms$transcript[1]

plot_df <- cancer_isoforms %>%
  dplyr::filter(transcript == SPTAN1_hit)

ggplot(
  plot_df,
  aes(
    x = state_label,
    y = fraction,
    fill = state_label
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.7) +
  theme_classic() +
  labs(
    x = NULL,
    y = "SPTAN1 isoform fraction",
    title = paste("SPTAN1", SPTAN1_hit, "across ADD3 cancer states")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

# Gerekirse:
# install.packages("mclust")

library(mclust)
library(dplyr)
library(tidyr)
library(ggplot2)

# SPTAN1 değerini sample bazında al
sptan1_df <- iso_selected %>%
  dplyr::filter(transcript == "ENST00000358161") %>%
  dplyr::select(sample, SPTAN1_fraction = fraction)

# ADD3 + ACTN1 + SPTAN1
prob_df <- model_df %>%
  dplyr::left_join(sptan1_df, by = "sample") %>%
  tidyr::drop_na(
    ADD3_transition_score,
    ENST00000553882,
    SPTAN1_fraction
  )

# Standardize
X <- scale(
  prob_df[, c(
    "ADD3_transition_score",
    "ENST00000553882",
    "SPTAN1_fraction"
  )]
)

# 3-component probabilistic model
set.seed(123)

gmm <- Mclust(
  X,
  G = 3
)

# Cluster assignment
prob_df$cluster <- gmm$classification

# Her cluster'ın özelliklerine bak
cluster_summary <- prob_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n = dplyr::n(),
    ADD3 = mean(ADD3_transition_score),
    ACTN1 = mean(ENST00000553882),
    SPTAN1 = mean(SPTAN1_fraction),
    .groups = "drop"
  )

cluster_summary

# Posterior probabilities
prob_df$P_High <- gmm$z[, 1]
prob_df$P_Low  <- gmm$z[, 2]
prob_df$P_Normal <- gmm$z[, 3]

probabilities <- prob_df %>%
  dplyr::select(
    sample,
    source,
    ADD3_transition_score,
    ENST00000553882,
    SPTAN1_fraction,
    P_Normal,
    P_Low,
    P_High
  )

head(probabilities)

prob_long <- probabilities %>%
  dplyr::select(
    sample, source,
    P_Normal, P_Low, P_High
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("P_"),
    names_to = "State",
    values_to = "Probability"
  )

ggplot(
  prob_long,
  aes(
    x = reorder(sample, P_High),
    y = Probability,
    fill = State
  )
) +
  geom_col() +
  theme_classic() +
  labs(
    x = "Samples",
    y = "Posterior probability",
    title = "Probabilistic ADD3–ACTN1–SPTAN1 States"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

probabilities %>%
  dplyr::select(
    sample,
    source,
    P_Normal,
    P_Low,
    P_High
  ) %>%
  head(20) %>%
  print(width = Inf)

probabilities %>%
  head() %>%
  print(width = Inf)

library(dplyr)
library(tidyr)
library(ggplot2)

# Sample'ları High olasılığına göre sırala
sample_order <- probabilities %>%
  dplyr::arrange(P_High) %>%
  dplyr::pull(sample)

# Long format
prob_long <- probabilities %>%
  dplyr::select(sample, source, P_Normal, P_Low, P_High) %>%
  tidyr::pivot_longer(
    cols = c(P_Normal, P_Low, P_High),
    names_to = "State",
    values_to = "Probability"
  ) %>%
  dplyr::mutate(
    sample = factor(sample, levels = sample_order),
    State = factor(
      State,
      levels = c("P_Normal", "P_Low", "P_High")
    )
  )

# Grafik
ggplot(
  prob_long,
  aes(
    x = sample,
    y = Probability,
    fill = State
  )
) +
  geom_col(width = 1) +
  facet_grid(~source, scales = "free_x", space = "free_x") +
  theme_classic() +
  labs(
    x = "Samples",
    y = "Posterior probability",
    title = "Probabilistic Molecular States"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

library(mclust)
library(dplyr)

set.seed(123)

B <- 200

X_full <- scale(
  prob_df[, c(
    "ADD3_transition_score",
    "ENST00000553882",
    "SPTAN1_fraction"
  )]
)

stability <- matrix(
  NA,
  nrow = nrow(prob_df),
  ncol = B
)

rownames(stability) <- prob_df$sample

for (b in 1:B) {
  
  # %80 sample seç
  idx <- sample(
    seq_len(nrow(prob_df)),
    size = round(0.8 * nrow(prob_df)),
    replace = FALSE
  )
  
  X_boot <- X_full[idx, , drop = FALSE]
  
  fit <- tryCatch(
    Mclust(X_boot, G = 3, verbose = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    stability[idx, b] <- fit$classification
  }
}

co_cluster <- matrix(
  NA,
  nrow = nrow(prob_df),
  ncol = nrow(prob_df)
)

rownames(co_cluster) <- prob_df$sample
colnames(co_cluster) <- prob_df$sample

for (i in seq_len(nrow(prob_df))) {
  for (j in seq_len(nrow(prob_df))) {
    
    valid <- !is.na(stability[i, ]) &
             !is.na(stability[j, ])
    
    if (sum(valid) > 0) {
      co_cluster[i, j] <- mean(
        stability[i, valid] ==
        stability[j, valid]
      )
    }
  }
}

library(pheatmap)

pheatmap(
  co_cluster,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Bootstrap Stability of 3 Molecular States"
)

# Cluster isimleri
prob_df$State <- dplyr::recode(
  as.character(gmm$classification),
  "1" = "High ADD3",
  "2" = "Low ADD3",
  "3" = "Normal-like"
)

# Heatmap sırasıyla eşleştir
annotation <- data.frame(
  State = prob_df$State,
  Source = prob_df$source
)

rownames(annotation) <- prob_df$sample

# co_cluster ile aynı sample'lar
annotation <- annotation[rownames(co_cluster), , drop = FALSE]

pheatmap::pheatmap(
  co_cluster,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annotation,
  annotation_col = annotation,
  main = "Bootstrap Stability of 3 Molecular States"
)

probabilities <- probabilities %>%
  dplyr::mutate(
    Transition_score =
      0 * P_Normal +
      0.5 * P_Low +
      1 * P_High
  )

# Sample'ları sırala
transition_plot <- probabilities %>%
  dplyr::arrange(Transition_score) %>%
  dplyr::mutate(
    order = dplyr::row_number()
  )

ggplot(
  transition_plot,
  aes(
    x = order,
    y = Transition_score,
    color = source
  )
) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = c(0.25, 0.75),
    linetype = "dashed"
  ) +
  theme_classic() +
  labs(
    x = "Samples ordered by transition score",
    y = "Probabilistic transition score",
    title = "ADD3–ACTN1–SPTAN1 Molecular Transition"
  )
probabilities %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median = median(Transition_score),
    mean = mean(Transition_score),
    min = min(Transition_score),
    max = max(Transition_score)
  )

state_map <- prob_df %>%
  dplyr::mutate(
    molecular_state = dplyr::recode(
      as.character(gmm$classification),
      "3" = "Normal",
      "2" = "Low",
      "1" = "High"
    )
  ) %>%
  dplyr::select(sample, molecular_state)
progression_isoforms <- iso_selected %>%
  dplyr::filter(
    !gene %in% c("ADD3", "ACTN1", "SPTAN1")
  ) %>%
  dplyr::inner_join(
    state_map,
    by = "sample"
  ) %>%
  dplyr::group_by(gene, transcript) %>%
  dplyr::summarise(
    Normal = median(
      fraction[molecular_state == "Normal"],
      na.rm = TRUE
    ),
    Low = median(
      fraction[molecular_state == "Low"],
      na.rm = TRUE
    ),
    High = median(
      fraction[molecular_state == "High"],
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    delta_Normal_Low = Low - Normal,
    delta_Low_High = High - Low,

    pattern = dplyr::case_when(
      Normal < Low & Low < High ~ "Increasing",
      Normal > Low & Low > High ~ "Decreasing",
      TRUE ~ "Non-monotonic"
    ),

    total_change = High - Normal
  ) %>%
  dplyr::arrange(desc(abs(total_change)))

progression_isoforms %>%
  dplyr::select(
    gene, transcript,
    Normal, Low, High,
    pattern,
    total_change
  ) %>%
  print(n = 35)

trend_data <- iso_selected %>%
  dplyr::filter(
    transcript %in% progression_isoforms %>%
      dplyr::filter(pattern %in% c("Increasing", "Decreasing")) %>%
      dplyr::pull(transcript)
  ) %>%
  dplyr::inner_join(state_map, by = "sample") %>%
  dplyr::mutate(
    state_numeric = dplyr::recode(
      molecular_state,
      "Normal" = 0,
      "Low" = 1,
      "High" = 2
    )
  )

trend_test <- trend_data %>%
  dplyr::group_by(gene, transcript) %>%
  dplyr::summarise(
    rho = cor(
      fraction,
      state_numeric,
      method = "spearman",
      use = "complete.obs"
    ),
    p_value = cor.test(
      fraction,
      state_numeric,
      method = "spearman",
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(FDR)

trend_test

# Monoton değişen transcriptleri önce seç
trend_transcripts <- progression_isoforms %>%
  dplyr::filter(
    pattern %in% c("Increasing", "Decreasing")
  ) %>%
  dplyr::pull(transcript)

trend_transcripts

trend_data <- iso_selected %>%
  dplyr::filter(transcript %in% trend_transcripts) %>%
  dplyr::inner_join(
    state_map,
    by = "sample"
  ) %>%
  dplyr::mutate(
    state_numeric = dplyr::case_when(
      molecular_state == "Normal" ~ 0,
      molecular_state == "Low" ~ 1,
      molecular_state == "High" ~ 2
    )
  )

trend_test <- trend_data %>%
  dplyr::group_by(gene, transcript) %>%
  dplyr::summarise(
    rho = cor(
      fraction,
      state_numeric,
      method = "spearman",
      use = "complete.obs"
    ),
    
    p_value = cor.test(
      fraction,
      state_numeric,
      method = "spearman",
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(FDR)

trend_test

cancer_trend <- trend_data %>%
  dplyr::filter(
    molecular_state %in% c("Low", "High")
  ) %>%
  dplyr::mutate(
    cancer_state = ifelse(
      molecular_state == "Low", 0, 1
    )
  ) %>%
  dplyr::group_by(gene, transcript) %>%
  dplyr::summarise(
    
    median_low = median(
      fraction[cancer_state == 0],
      na.rm = TRUE
    ),
    
    median_high = median(
      fraction[cancer_state == 1],
      na.rm = TRUE
    ),
    
    delta = median_high - median_low,
    
    p_value = wilcox.test(
      fraction ~ cancer_state,
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(FDR)

cancer_trend

cancer_trend %>%
  dplyr::select(
    gene, transcript,
    median_low, median_high,
    delta, p_value, FDR
  ) %>%
  print(width = Inf)

library(dplyr)
library(tidyr)
library(ggplot2)

# Görselde göstermek istediğimiz temel isoformlar
key_isoforms <- c(
  "ENST00000356080", # ADD3 long
  "ENST00000360162", # ADD3 short
  "ENST00000553882", # ACTN1
  "ENST00000358161", # SPTAN1
  "ENST00000265800"  # DMTN
)

# State bilgisi ile birleştir
trajectory_df <- iso_long %>%
  dplyr::filter(
    transcript %in% key_isoforms,
    source %in% c("gtex_ovary", "ccle_cancer")
  ) %>%
  dplyr::inner_join(
    state_map,
    by = "sample"
  ) %>%
  dplyr::mutate(
    label = dplyr::case_when(
      transcript == "ENST00000356080" ~ "ADD3 Long",
      transcript == "ENST00000360162" ~ "ADD3 Short",
      transcript == "ENST00000553882" ~ "ACTN1",
      transcript == "ENST00000358161" ~ "SPTAN1",
      transcript == "ENST00000265800" ~ "DMTN"
    ),
    molecular_state = factor(
      molecular_state,
      levels = c("Normal", "Low", "High")
    )
  )

# Median fraction
trajectory_summary <- trajectory_df %>%
  dplyr::group_by(label, molecular_state) %>%
  dplyr::summarise(
    median_fraction = median(fraction, na.rm = TRUE),
    .groups = "drop"
  )

# Çiz
ggplot(
  trajectory_summary,
  aes(
    x = molecular_state,
    y = median_fraction,
    group = label,
    color = label
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  theme_classic() +
  labs(
    x = "Data-driven molecular state",
    y = "Median isoform fraction",
    color = "Isoform",
    title = "Isoform Changes Across ADD3-defined Molecular States"
  )

ggplot(
  trajectory_df,
  aes(
    x = molecular_state,
    y = fraction
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.35,
    size = 1
  ) +
  facet_wrap(
    ~label,
    scales = "free_y",
    ncol = 3
  ) +
  theme_classic() +
  labs(
    x = "Molecular state",
    y = "Isoform fraction",
    title = "Sample-level Isoform Patterns Across Molecular States"
  )
all_flow <- iso_selected %>%
  dplyr::inner_join(state_map, by = "sample") %>%
  dplyr::mutate(
    molecular_state = factor(
      molecular_state,
      levels = c("Normal", "Low", "High")
    )
  ) %>%
  dplyr::group_by(gene, transcript, molecular_state) %>%
  dplyr::summarise(
    median_fraction = median(fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = paste0(gene, "_", transcript)
  ) %>%
  tidyr::pivot_wider(
    names_from = molecular_state,
    values_from = median_fraction
  )

mat <- all_flow %>%
  dplyr::select(Normal, Low, High) %>%
  as.matrix()

rownames(mat) <- all_flow$label

pheatmap::pheatmap(
  mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 7,
  main = "Isoform Flow: Normal → Low ADD3 → High ADD3"
)

late_transition <- iso_selected %>%
  dplyr::inner_join(state_map, by = "sample") %>%
  dplyr::filter(molecular_state %in% c("Low", "High")) %>%
  dplyr::group_by(gene, transcript) %>%
  dplyr::summarise(
    median_low = median(
      fraction[molecular_state == "Low"],
      na.rm = TRUE
    ),
    median_high = median(
      fraction[molecular_state == "High"],
      na.rm = TRUE
    ),
    delta_LH = median_high - median_low,

    p_value = wilcox.test(
      fraction ~ molecular_state,
      exact = FALSE
    )$p.value,

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(desc(abs(delta_LH)))

late_transition %>%
  dplyr::select(
    gene, transcript,
    median_low, median_high,
    delta_LH, p_value, FDR
  ) %>%
  print(n = 35, width = Inf)

late_sig <- c(
  "ENST00000358161", # SPTAN1 DOWN
  "ENST00000372731", # SPTAN1 UP
  "ENST00000356080", # ADD3 Long UP
  "ENST00000523782"  # DMTN UP
)

late_sig_df <- iso_selected %>%
  dplyr::filter(transcript %in% late_sig) %>%
  dplyr::inner_join(state_map, by = "sample") %>%
  dplyr::select(sample, molecular_state, transcript, fraction) %>%
  tidyr::pivot_wider(
    names_from = transcript,
    values_from = fraction
  )

late_sig_df <- late_sig_df %>%
  dplyr::mutate(
    Late_transition_score =
      scale(ENST00000372731)[,1] +
      scale(ENST00000356080)[,1] +
      scale(ENST00000523782)[,1] -
      scale(ENST00000358161)[,1]
  )

ggplot(
  late_sig_df,
  aes(
    x = molecular_state,
    y = Late_transition_score
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Molecular state",
    y = "Late-transition signature score",
    title = "Multi-Isoform Late Transition Signature"
  )

library(dplyr)
library(ggplot2)

# Late transition için anlamlı isoformlar
signature_isoforms <- late_transition %>%
  filter(FDR < 0.05) %>%
  mutate(
    weight = abs(delta_LH),
    direction = ifelse(delta_LH > 0, "Increase", "Decrease"),
    label = paste0(gene, "_", transcript)
  ) %>%
  arrange(desc(weight))

signature_isoforms %>%
  select(gene, transcript, delta_LH, FDR, weight)
ggplot(
  signature_isoforms,
  aes(
    x = reorder(label, weight),
    y = weight,
    fill = direction
  )
) +
  geom_col(width = 0.7) +
  coord_flip() +
  labs(
    x = NULL,
    y = "|Low → High change|",
    fill = "Direction",
    title = "Contribution to the Late-Transition Signature"
  ) +
  theme_classic(base_size = 14)
signature_isoforms %>%
  dplyr::select(
    gene,
    transcript,
    delta_LH,
    FDR,
    weight
  )
library(dplyr)
library(ggplot2)

# 1. Late-transition signature'a giren anlamlı isoformları seç
signature_isoforms <- late_transition %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::mutate(
    weight = abs(delta_LH),
    direction = dplyr::if_else(
      delta_LH > 0,
      "Increase",
      "Decrease"
    ),
    label = paste0(gene, "_", transcript)
  ) %>%
  dplyr::arrange(dplyr::desc(weight))


# 2. Sonuçları tablo olarak göster
signature_isoforms %>%
  dplyr::select(
    gene,
    transcript,
    median_low,
    median_high,
    delta_LH,
    FDR,
    weight,
    direction
  ) %>%
  print(n = Inf, width = Inf)


# 3. Her isoformun late-transition signature'a katkısını görselleştir
ggplot2::ggplot(
  signature_isoforms,
  ggplot2::aes(
    x = reorder(label, weight),
    y = weight,
    fill = direction
  )
) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "|Low → High change|",
    fill = "Direction",
    title = "Contribution to the Late-Transition Signature"
  ) +
  ggplot2::theme_classic(base_size = 14)
library(dplyr)
library(tidyr)
library(ggplot2)

# Core signature
core_ids <- c(
  "ENST00000358161",  # SPTAN1 decrease
  "ENST00000372731",  # SPTAN1 increase
  "ENST00000356080",  # ADD3 increase
  "ENST00000523782"   # DMTN increase
)

# Sample x isoform matrix
core_df <- iso_long %>%
  dplyr::filter(transcript %in% core_ids) %>%
  dplyr::select(sample, source, transcript, fraction) %>%
  tidyr::pivot_wider(
    names_from = transcript,
    values_from = fraction
  )

# Önceden oluşturduğumuz molecular state bilgisini ekle
core_df <- core_df %>%
  dplyr::left_join(
    cluster_df %>%
      dplyr::select(sample, state),
    by = "sample"
  )

# State isimlerini mevcut state tanımımıza göre ver
core_df <- core_df %>%
  dplyr::mutate(
    molecular_state = dplyr::case_when(
      state == 1 ~ "High",
      state == 2 ~ "Low",
      state == 3 ~ "Normal",
      TRUE ~ NA_character_
    )
  )

# Isoformları standardize et
core_df <- core_df %>%
  dplyr::mutate(
    z_SPTAN1_down = as.numeric(scale(ENST00000358161)),
    z_SPTAN1_up   = as.numeric(scale(ENST00000372731)),
    z_ADD3_up     = as.numeric(scale(ENST00000356080)),
    z_DMTN_up     = as.numeric(scale(ENST00000523782))
  )

# Core transition score
core_df <- core_df %>%
  dplyr::mutate(
    core_transition_score =
      -z_SPTAN1_down +
       z_SPTAN1_up +
       z_ADD3_up +
       z_DMTN_up
  )


# Sadece cancer Low vs High
core_cancer <- core_df %>%
  dplyr::filter(
    source == "ccle_cancer",
    molecular_state %in% c("Low", "High")
  )


# Özet
core_cancer %>%
  dplyr::group_by(molecular_state) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median = median(core_transition_score, na.rm = TRUE),
    mean = mean(core_transition_score, na.rm = TRUE),
    sd = sd(core_transition_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()


# Low vs High testi
wilcox.test(
  core_transition_score ~ molecular_state,
  data = core_cancer,
  exact = FALSE
)


# Görselleştir
ggplot2::ggplot(
  core_cancer,
  ggplot2::aes(
    x = molecular_state,
    y = core_transition_score
  )
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.10,
    size = 2,
    alpha = 0.7
  ) +
  ggplot2::labs(
    x = "Cancer molecular state",
    y = "4-isoform core transition score",
    title = "Core Isoform Signature: Low → High Transition"
  ) +
  ggplot2::theme_classic(base_size = 14)
# Gerekirse bir kez:
# install.packages("pROC")

library(pROC)

roc_core <- pROC::roc(
  response = core_cancer$molecular_state,
  predictor = core_cancer$core_transition_score,
  levels = c("Low", "High"),
  direction = "<"
)

# AUC
pROC::auc(roc_core)

# %95 CI
pROC::ci.auc(roc_core)

# ROC plot
plot(
  roc_core,
  print.auc = TRUE,
  legacy.axes = TRUE,
  main = "4-Isoform Core Signature: Low vs High"
)
# Gerekirse bir kez:
# install.packages("nnet")

library(dplyr)
library(ggplot2)
library(nnet)
library(tidyr)

# --------------------------------------------------
# 1. ACTN1 bilgisini core dataframe'e ekle
# --------------------------------------------------

actn1_df <- iso_long %>%
  dplyr::filter(transcript == "ENST00000553882") %>%
  dplyr::select(
    sample,
    ACTN1_fraction = fraction
  )


prob_model_df <- core_df %>%
  dplyr::left_join(
    actn1_df,
    by = "sample"
  ) %>%
  dplyr::filter(
    !is.na(molecular_state),
    !is.na(core_transition_score),
    !is.na(ACTN1_fraction)
  )


# --------------------------------------------------
# 2. State sırasını tanımla
# --------------------------------------------------

prob_model_df$molecular_state <- factor(
  prob_model_df$molecular_state,
  levels = c("Normal", "Low", "High")
)


# --------------------------------------------------
# 3. Değişkenleri standardize et
# --------------------------------------------------

prob_model_df <- prob_model_df %>%
  dplyr::mutate(
    z_core = as.numeric(scale(core_transition_score)),
    z_ACTN1 = as.numeric(scale(ACTN1_fraction))
  )


# --------------------------------------------------
# 4. Multinomial probabilistic model
# --------------------------------------------------

prob_model <- nnet::multinom(
  molecular_state ~ z_core + z_ACTN1,
  data = prob_model_df,
  trace = FALSE
)


# Model özeti
summary(prob_model)


# --------------------------------------------------
# 5. Her sample için state probabilities
# --------------------------------------------------

state_prob <- predict(
  prob_model,
  newdata = prob_model_df,
  type = "probs"
)


prob_model_df$P_Normal <- state_prob[, "Normal"]
prob_model_df$P_Low    <- state_prob[, "Low"]
prob_model_df$P_High   <- state_prob[, "High"]


# En olası state
prob_model_df <- prob_model_df %>%
  dplyr::mutate(
    predicted_state = c(
      "Normal",
      "Low",
      "High"
    )[max.col(
      cbind(P_Normal, P_Low, P_High),
      ties.method = "first"
    )]
  )


# --------------------------------------------------
# 6. İlk sonuçlara bak
# --------------------------------------------------

prob_model_df %>%
  dplyr::select(
    sample,
    source,
    molecular_state,
    core_transition_score,
    ACTN1_fraction,
    P_Normal,
    P_Low,
    P_High,
    predicted_state
  ) %>%
  head(15) %>%
  print(width = Inf)


# --------------------------------------------------
# 7. Probability flow görseli
# --------------------------------------------------

prob_plot <- prob_model_df %>%
  dplyr::arrange(P_High) %>%
  dplyr::mutate(
    sample_order = dplyr::row_number()
  ) %>%
  dplyr::select(
    sample_order,
    P_Normal,
    P_Low,
    P_High
  ) %>%
  tidyr::pivot_longer(
    cols = c(P_Normal, P_Low, P_High),
    names_to = "State",
    values_to = "Probability"
  ) %>%
  dplyr::mutate(
    State = factor(
      State,
      levels = c(
        "P_Normal",
        "P_Low",
        "P_High"
      ),
      labels = c(
        "Normal",
        "Low-transition",
        "High-transition"
      )
    )
  )


ggplot2::ggplot(
  prob_plot,
  ggplot2::aes(
    x = sample_order,
    y = Probability,
    fill = State
  )
) +
  ggplot2::geom_area(
    position = "stack"
  ) +
  ggplot2::labs(
    x = "Samples ordered by P(High)",
    y = "State probability",
    fill = "Molecular state",
    title = "Probabilistic ADD3-centered Molecular State Model"
  ) +
  ggplot2::theme_classic(base_size = 14)


# --------------------------------------------------
# 8. Model classification table
# --------------------------------------------------

table(
  Observed = prob_model_df$molecular_state,
  Predicted = prob_model_df$predicted_state
)
library(dplyr)
library(pheatmap)

# --------------------------------------------------
# 1. Heatmap için gerekli değişkenleri hazırla
# --------------------------------------------------

heat_df <- prob_model_df %>%
  dplyr::arrange(P_High) %>%
  dplyr::select(
    sample,
    source,
    molecular_state,
    P_Normal,
    P_Low,
    P_High,
    ACTN1_fraction,
    ENST00000358161,   # SPTAN1 down
    ENST00000372731,   # SPTAN1 up
    ENST00000356080,   # ADD3 up
    ENST00000523782    # DMTN up
  )


# --------------------------------------------------
# 2. Sample isimlerini rowname yap
# --------------------------------------------------

heat_mat <- heat_df %>%
  dplyr::select(
    P_Normal,
    P_Low,
    P_High,
    ACTN1_fraction,
    ENST00000358161,
    ENST00000372731,
    ENST00000356080,
    ENST00000523782
  ) %>%
  as.data.frame()

rownames(heat_mat) <- heat_df$sample


# --------------------------------------------------
# 3. Daha anlaşılır kolon isimleri
# --------------------------------------------------

colnames(heat_mat) <- c(
  "P(Normal)",
  "P(Low)",
  "P(High)",
  "ACTN1_553882",
  "SPTAN1_358161",
  "SPTAN1_372731",
  "ADD3_356080",
  "DMTN_523782"
)


# --------------------------------------------------
# 4. Molecular measurements'ı z-score yap
# Probability kolonlarını değiştirmiyoruz
# --------------------------------------------------

heat_scaled <- heat_mat

heat_scaled[, 4:8] <- scale(
  heat_scaled[, 4:8]
)


# --------------------------------------------------
# 5. Sample annotation
# --------------------------------------------------

annotation_row <- data.frame(
  Source = heat_df$source,
  State = heat_df$molecular_state
)

rownames(annotation_row) <- heat_df$sample


# --------------------------------------------------
# 6. Heatmap
# --------------------------------------------------

pheatmap::pheatmap(
  heat_scaled,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  annotation_row = annotation_row,
  
  show_rownames = FALSE,
  show_colnames = TRUE,
  
  fontsize_col = 10,
  angle_col = 45,
  
  border_color = "grey80",
  
  main = "Molecular Architecture of Probabilistic Transition States"
)
library(dplyr)
library(tidyr)
library(bnlearn)

# ============================================================
# 1. TÜM SAMPLE-LEVEL ISOFORMLARDAN BAŞLA
# ============================================================

all_iso <- iso_long %>%
  dplyr::filter(
    source %in% c("ccle_cancer", "gtex_ovary")
  ) %>%
  dplyr::mutate(
    node = paste0(gene, "_", transcript)
  )

cat("Başlangıç isoform sayısı:",
    dplyr::n_distinct(all_iso$node), "\n")


# ============================================================
# 2. SADECE TEKNİK / UNSUPERVISED FİLTRELEME
#
# prevalence = örneklerin kaçında isoform fraction > 0.01
# variance   = örnekler arasında ne kadar değişiyor
#
# Cancer/normal bilgisi burada KULLANILMIYOR.
# ============================================================

iso_stats <- all_iso %>%
  dplyr::group_by(node) %>%
  dplyr::summarise(
    prevalence = mean(fraction > 0.01, na.rm = TRUE),
    variance = var(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    .groups = "drop"
  )


# En az %20 sample'da görülen ve değişkenliği sıfır olmayanlar
keep_iso <- iso_stats %>%
  dplyr::filter(
    prevalence >= 0.20,
    variance > 0
  )


cat("Prevalence filtresinden sonra:",
    nrow(keep_iso), "\n")


# ============================================================
# 3. BN İÇİN EN DEĞİŞKEN ISOFORMLARI AL
#
# 93 sample olduğu için binlerce node ile BN kurmuyoruz.
# İlk deneme: en değişken 50 isoform.
# ============================================================

top_iso <- keep_iso %>%
  dplyr::arrange(desc(variance)) %>%
  dplyr::slice_head(n = 50)

selected_nodes <- top_iso$node

cat("BN'e girecek isoform sayısı:",
    length(selected_nodes), "\n")

print(top_iso)


# ============================================================
# 4. SAMPLE × ISOFORM MATRIX
# ============================================================

bn_df <- all_iso %>%
  dplyr::filter(node %in% selected_nodes) %>%
  dplyr::select(sample, source, node, fraction) %>%
  tidyr::pivot_wider(
    names_from = node,
    values_from = fraction,
    values_fill = 0
  )


# source'u şimdilik network öğrenmeye SOKMUYORUZ
sample_info <- bn_df %>%
  dplyr::select(sample, source)

X <- bn_df %>%
  dplyr::select(-sample, -source) %>%
  as.data.frame()


# bnlearn'in node isimlerinde problem yaşamaması için
original_names <- colnames(X)
safe_names <- make.names(original_names, unique = TRUE)

colnames(X) <- safe_names

name_map <- data.frame(
  safe = safe_names,
  original = original_names
)

cat("Matrix dimensions:",
    nrow(X), "samples x",
    ncol(X), "isoforms\n")


# ============================================================
# 5. STANDARDIZE
# ============================================================

X_scaled <- as.data.frame(scale(X))

# scale sonrası sorunlu kolon varsa çıkar
valid <- sapply(
  X_scaled,
  function(x)
    all(is.finite(x)) &&
    sd(x, na.rm = TRUE) > 0
)

X_scaled <- X_scaled[, valid, drop = FALSE]

cat("Final BN dimensions:",
    nrow(X_scaled), "x",
    ncol(X_scaled), "\n")


# ============================================================
# 6. DATA-DRIVEN BAYESIAN NETWORK
#
# Gaussian BN:
# Hill-Climbing + BIC
#
# Hiçbir edge'i önceden dayatmıyoruz.
# ============================================================

set.seed(123)

bn_full <- bnlearn::hc(
  X_scaled,
  score = "bic-g"
)

bn_full

cat("Learned edges:",
    nrow(bnlearn::arcs(bn_full)), "\n")


# ============================================================
# 7. BOOTSTRAP STABILITY
#
# Asıl önemli kısım burası.
# Tek bir network yerine 500 bootstrap network öğreniyoruz.
# ============================================================

set.seed(123)

boot_bn <- bnlearn::boot.strength(
  data = X_scaled,
  R = 500,
  algorithm = "hc",
  algorithm.args = list(
    score = "bic-g"
  )
)

head(
  boot_bn[
    order(-boot_bn$strength),
  ],
  20
)


# ============================================================
# 8. STABLE NETWORK
#
# strength >= 0.70:
# bootstrapların en az %70'inde görülen bağlantılar
# ============================================================

stable_bn <- bnlearn::averaged.network(
  boot_bn,
  threshold = 0.70
)

stable_bn

cat(
  "Stable edges:",
  nrow(bnlearn::arcs(stable_bn)),
  "\n"
)


# ============================================================
# 9. EDGE TABLOSUNU GERÇEK GENE_ENST İSİMLERİNE ÇEVİR
# ============================================================

stable_edges <- as.data.frame(
  bnlearn::arcs(stable_bn)
)

stable_edges <- stable_edges %>%
  dplyr::left_join(
    name_map,
    by = c("from" = "safe")
  ) %>%
  dplyr::rename(
    from_original = original
  ) %>%
  dplyr::left_join(
    name_map,
    by = c("to" = "safe")
  ) %>%
  dplyr::rename(
    to_original = original
  )


# Bootstrap strength ekle
boot_table <- boot_bn %>%
  as.data.frame()

stable_edges <- stable_edges %>%
  dplyr::left_join(
    boot_table,
    by = c("from", "to")
  )


stable_edges %>%
  dplyr::select(
    from_original,
    to_original,
    strength,
    direction
  ) %>%
  dplyr::arrange(desc(strength)) %>%
  print(n = 100)


# ============================================================
# 10. ADD3 ÇEVRESİNDE NE ÇIKTI?
# ============================================================

ADD3_edges <- stable_edges %>%
  dplyr::filter(
    grepl("^ADD3_", from_original) |
    grepl("^ADD3_", to_original)
  ) %>%
  dplyr::select(
    from_original,
    to_original,
    strength,
    direction
  ) %>%
  dplyr::arrange(desc(strength))

ADD3_edges %>%
  print(n = Inf)


# ============================================================
# 11. KAYDET
# ============================================================

write.csv(
  stable_edges,
  "/Users/lemannur/Downloads/all_isoform_BN_stable_edges.csv",
  row.names = FALSE
)

write.csv(
  ADD3_edges,
  "/Users/lemannur/Downloads/ADD3_BN_stable_edges.csv",
  row.names = FALSE
)
stable_edges_show <- stable_edges %>%
  dplyr::select(
    from_original,
    to_original,
    strength,
    direction
  ) %>%
  dplyr::arrange(dplyr::desc(strength)) %>%
  as.data.frame()

utils::head(stable_edges_show, 100)
ADD3_edges <- stable_edges_show %>%
  dplyr::filter(
    grepl("^ADD3_", from_original) |
    grepl("^ADD3_", to_original)
  )

View(ADD3_edges)
nrow(keep_iso)
head(top_iso, 20)
ADD3_edges
library(igraph)
library(dplyr)

# Stable edge tablosu
plot_edges <- stable_edges %>%
  dplyr::filter(strength >= 0.70) %>%
  dplyr::select(
    from_original,
    to_original,
    strength
  ) %>%
  as.data.frame()

# igraph
g_all <- igraph::graph_from_data_frame(
  plot_edges,
  directed = TRUE
)

# ENST numaralarını biraz kısalt
V(g_all)$label <- gsub(
  "ENST00000",
  "",
  V(g_all)$name
)

# Edge kalınlığı = bootstrap strength
E(g_all)$width <- 1 + 4 * E(g_all)$strength

# Node büyüklüğü
V(g_all)$size <- 8

# ADD3 / ACTN1 / SPTAN1 / DMTN biraz büyük göster
key_genes <- grepl(
  "^(ADD3|ACTN1|SPTAN1|DMTN)_",
  V(g_all)$name
)

V(g_all)$size[key_genes] <- 13

set.seed(123)

plot(
  g_all,
  layout = igraph::layout_with_fr(g_all),
  vertex.label = V(g_all)$label,
  vertex.label.cex = 0.55,
  vertex.size = V(g_all)$size,
  edge.width = E(g_all)$width,
  edge.arrow.size = 0.35,
  main = "Data-driven Isoform Bayesian Network"
)
# ADD3 node'larını bul
add3_nodes <- V(g_all)$name[
  grepl("^ADD3_", V(g_all)$name)
]

# ADD3'e 2 adım uzaklıktaki bütün node'ları al
near_ADD3 <- unique(
  unlist(
    igraph::ego(
      g_all,
      order = 2,
      nodes = add3_nodes,
      mode = "all"
    )
  )
)

# Subnetwork
g_ADD3 <- igraph::induced_subgraph(
  g_all,
  vids = near_ADD3
)

# Label
V(g_ADD3)$label <- gsub(
  "ENST00000",
  "",
  V(g_ADD3)$name
)

# Key genleri büyük göster
key <- grepl(
  "^(ADD3|ACTN1|SPTAN1|DMTN)_",
  V(g_ADD3)$name
)

V(g_ADD3)$size <- 10
V(g_ADD3)$size[key] <- 17

E(g_ADD3)$width <- 1 + 5 * E(g_ADD3)$strength

set.seed(123)

plot(
  g_ADD3,
  layout = igraph::layout_with_fr(g_ADD3),
  vertex.label.cex = 0.7,
  vertex.size = V(g_ADD3)$size,
  edge.width = E(g_ADD3)$width,
  edge.arrow.size = 0.4,
  main = "ADD3-centered Bayesian Isoform Network"
)
library(dplyr)
library(tidyr)
library(bnlearn)

# ============================================================
# 1. AYNI ISOFORMLARI KULLAN
# ============================================================

selected_nodes <- top_iso$node

make_bn_matrix <- function(source_name) {
  
  tmp <- all_iso %>%
    filter(
      source == source_name,
      node %in% selected_nodes
    ) %>%
    select(sample, node, fraction) %>%
    pivot_wider(
      names_from = node,
      values_from = fraction,
      values_fill = 0
    )
  
  X <- tmp %>%
    select(-sample) %>%
    as.data.frame()
  
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  
  # constant/problematic variables remove
  valid <- sapply(
    X,
    function(x) sd(x, na.rm = TRUE) > 0
  )
  
  X <- X[, valid, drop = FALSE]
  
  as.data.frame(scale(X))
}


X_normal <- make_bn_matrix("gtex_ovary")
X_cancer <- make_bn_matrix("ccle_cancer")

dim(X_normal)
dim(X_cancer)
set.seed(123)

boot_normal <- boot.strength(
  X_normal,
  R = 500,
  algorithm = "hc",
  algorithm.args = list(score = "bic-g")
)

set.seed(123)

boot_cancer <- boot.strength(
  X_cancer,
  R = 500,
  algorithm = "hc",
  algorithm.args = list(score = "bic-g")
)


bn_normal <- averaged.network(
  boot_normal,
  threshold = 0.70
)

bn_cancer <- averaged.network(
  boot_cancer,
  threshold = 0.70
)
cat(
  "Normal stable edges:",
  nrow(arcs(bn_normal)),
  "\n"
)

cat(
  "Cancer stable edges:",
  nrow(arcs(bn_cancer)),
  "\n"
)
normal_edges <- as.data.frame(arcs(bn_normal)) %>%
  mutate(
    edge = paste(
      pmin(from, to),
      pmax(from, to),
      sep = " -- "
    )
  )

cancer_edges <- as.data.frame(arcs(bn_cancer)) %>%
  mutate(
    edge = paste(
      pmin(from, to),
      pmax(from, to),
      sep = " -- "
    )
  )


common_edges <- intersect(
  normal_edges$edge,
  cancer_edges$edge
)

normal_only <- setdiff(
  normal_edges$edge,
  cancer_edges$edge
)

cancer_only <- setdiff(
  cancer_edges$edge,
  normal_edges$edge
)


cat("Common:", length(common_edges), "\n")
cat("Normal only:", length(normal_only), "\n")
cat("Cancer only:", length(cancer_only), "\n")
cat("Common:", length(common_edges), "\n")
cat("Normal only:", length(normal_only), "\n")
cat("Cancer only:", length(cancer_only), "\n")
make_bn_matrix <- function(source_name) {
  
  tmp <- all_iso %>%
    dplyr::filter(
      source == source_name,
      node %in% selected_nodes
    ) %>%
    dplyr::select(sample, node, fraction) %>%
    tidyr::pivot_wider(
      names_from = node,
      values_from = fraction,
      values_fill = 0
    )
  
  X <- tmp %>%
    dplyr::select(-sample) %>%
    as.data.frame()
  
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  
  valid <- sapply(
    X,
    function(x) sd(x, na.rm = TRUE) > 0
  )
  
  X <- X[, valid, drop = FALSE]
  
  as.data.frame(scale(X))
}

X_normal <- make_bn_matrix("gtex_ovary")
X_cancer <- make_bn_matrix("ccle_cancer")

dim(X_normal)
dim(X_cancer)
set.seed(123)

boot_normal <- bnlearn::boot.strength(
  X_normal,
  R = 500,
  algorithm = "hc",
  algorithm.args = list(score = "bic-g")
)

set.seed(123)

boot_cancer <- bnlearn::boot.strength(
  X_cancer,
  R = 500,
  algorithm = "hc",
  algorithm.args = list(score = "bic-g")
)

bn_normal <- bnlearn::averaged.network(
  boot_normal,
  threshold = 0.70
)

bn_cancer <- bnlearn::averaged.network(
  boot_cancer,
  threshold = 0.70
)

cat("Normal stable edges:", nrow(bnlearn::arcs(bn_normal)), "\n")
cat("Cancer stable edges:", nrow(bnlearn::arcs(bn_cancer)), "\n")
normal_edges <- as.data.frame(bnlearn::arcs(bn_normal)) %>%
  dplyr::mutate(
    edge = paste(
      pmin(from, to),
      pmax(from, to),
      sep = " -- "
    )
  )

cancer_edges <- as.data.frame(bnlearn::arcs(bn_cancer)) %>%
  dplyr::mutate(
    edge = paste(
      pmin(from, to),
      pmax(from, to),
      sep = " -- "
    )
  )

common_edges <- intersect(normal_edges$edge, cancer_edges$edge)
normal_only  <- setdiff(normal_edges$edge, cancer_edges$edge)
cancer_only  <- setdiff(cancer_edges$edge, normal_edges$edge)

cat("Common:", length(common_edges), "\n")
cat("Normal only:", length(normal_only), "\n")
cat("Cancer only:", length(cancer_only), "\n")
library(igraph)
library(dplyr)

# ============================================================
# 1. EDGE'LERİ SINIFLANDIR
# ============================================================

rewire_edges <- dplyr::bind_rows(

  normal_edges %>%
    dplyr::filter(edge %in% normal_only) %>%
    dplyr::transmute(from, to, type = "Normal only"),

  cancer_edges %>%
    dplyr::filter(edge %in% cancer_only) %>%
    dplyr::transmute(from, to, type = "Cancer only"),

  cancer_edges %>%
    dplyr::filter(edge %in% common_edges) %>%
    dplyr::transmute(from, to, type = "Common")
)


# ============================================================
# 2. NETWORK
# ============================================================

g_rewire <- igraph::graph_from_data_frame(
  rewire_edges,
  directed = FALSE
)


# ============================================================
# 3. NODE LABEL'LARINI DÜZELT
# ============================================================

# make.names() ile değişmiş isimleri tekrar okunabilir yap
V(g_rewire)$label <- V(g_rewire)$name

V(g_rewire)$label <- gsub(
  "\\.ENST",
  "_ENST",
  V(g_rewire)$label
)

V(g_rewire)$label <- gsub(
  "ENST00000",
  "",
  V(g_rewire)$label
)


# ============================================================
# 4. KEY GENLERİ BÜYÜT
# ============================================================

key_nodes <- grepl(
  "^(ADD3|ACTN1|SPTAN1|DMTN)",
  V(g_rewire)$name
)

V(g_rewire)$size <- 8
V(g_rewire)$size[key_nodes] <- 15


# ============================================================
# 5. EDGE TİPLERİ
# ============================================================

edge_cols <- c(
  "Normal only" = "steelblue",
  "Cancer only" = "firebrick",
  "Common" = "grey40"
)

E(g_rewire)$color <- edge_cols[
  E(g_rewire)$type
]

E(g_rewire)$width <- ifelse(
  E(g_rewire)$type == "Common",
  3,
  2
)


# ============================================================
# 6. PLOT
# ============================================================

set.seed(123)

plot(
  g_rewire,
  layout = igraph::layout_with_fr(g_rewire),

  vertex.label = V(g_rewire)$label,
  vertex.label.cex = 0.55,

  vertex.size = V(g_rewire)$size,

  edge.color = E(g_rewire)$color,
  edge.width = E(g_rewire)$width,

  main = "Bayesian Network Rewiring: Normal vs Cancer"
)

legend(
  "topleft",
  legend = c(
    "Normal only",
    "Cancer only",
    "Common"
  ),
  col = c(
    "steelblue",
    "firebrick",
    "grey40"
  ),
  lwd = c(2,2,3),
  bty = "n"
)
# Daha temiz layout
set.seed(123)

lay <- igraph::layout_with_kk(g_rewire)

# label'ları daha kısa yap
V(g_rewire)$label <- sub(
  ".*_ENST00000",
  "",
  V(g_rewire)$name
)

# gene adını ekle
V(g_rewire)$label <- paste0(
  sub("_.*", "", V(g_rewire)$name),
  "_",
  V(g_rewire)$label
)

# Node büyüklüğü
V(g_rewire)$size <- 7

key <- grepl(
  "^(ADD3|ACTN1|SPTAN1|DMTN)",
  V(g_rewire)$name
)

V(g_rewire)$size[key] <- 11

plot(
  g_rewire,
  layout = lay,
  vertex.size = V(g_rewire)$size,
  vertex.label = V(g_rewire)$label,
  vertex.label.cex = 0.55,
  vertex.label.dist = 0.8,
  edge.color = E(g_rewire)$color,
  edge.width = 2,
  edge.curved = 0.08,
  main = "Isoform Network Rewiring"
)

legend(
  "topleft",
  legend = c("Normal only", "Cancer only", "Common"),
  col = c("steelblue", "firebrick", "grey40"),
  lwd = 3,
  bty = "n"
)
library(igraph)

# ------------------------------------------------------------
# Ortak node evreni
# ------------------------------------------------------------

all_nodes <- union(
  unique(c(normal_edges$from, normal_edges$to)),
  unique(c(cancer_edges$from, cancer_edges$to))
)

# İki network
gN <- graph_from_data_frame(
  normal_edges[, c("from","to")],
  directed = FALSE,
  vertices = all_nodes
)

gC <- graph_from_data_frame(
  cancer_edges[, c("from","to")],
  directed = FALSE,
  vertices = all_nodes
)

# ------------------------------------------------------------
# ORTAK LAYOUT
# İki networkün union'ından koordinat öğren
# ------------------------------------------------------------

union_edges <- unique(
  rbind(
    normal_edges[, c("from","to")],
    cancer_edges[, c("from","to")]
  )
)

g_union <- graph_from_data_frame(
  union_edges,
  directed = FALSE,
  vertices = all_nodes
)

set.seed(123)

coords <- layout_with_fr(
  g_union,
  niter = 3000
)

rownames(coords) <- V(g_union)$name

# ------------------------------------------------------------
# Kısa label
# ------------------------------------------------------------

short_label <- function(x) {

  gene <- sub("\\..*$", "", x)

  enst <- sub("^.*ENST00000", "", x)

  paste0(gene, "\n", enst)
}

# Sadece önemli genleri label'la
important <- c(
  "ADD3",
  "ACTN1",
  "SPTAN1",
  "DMTN",
  "DIAPH1",
  "ESYT1"
)

make_labels <- function(g) {

  gene <- sub("\\..*$", "", V(g)$name)

  ifelse(
    gene %in% important,
    short_label(V(g)$name),
    ""
  )
}

# ------------------------------------------------------------
# SIDE-BY-SIDE
# ------------------------------------------------------------

par(
  mfrow = c(1,2),
  mar = c(1,1,3,1)
)

plot(
  gN,
  layout = coords[V(gN)$name, ],
  vertex.size = 8,
  vertex.label = make_labels(gN),
  vertex.label.cex = 0.7,
  edge.width = 1.5,
  edge.arrow.mode = 0,
  main = "Normal"
)

plot(
  gC,
  layout = coords[V(gC)$name, ],
  vertex.size = 8,
  vertex.label = make_labels(gC),
  vertex.label.cex = 0.7,
  edge.width = 1.5,
  edge.arrow.mode = 0,
  main = "Cancer"
)

par(mfrow = c(1,1))
# Her node'un bağlantı sayısı
deg_normal <- igraph::degree(gN)
deg_cancer <- igraph::degree(gC)

rewiring_score <- data.frame(
  node = all_nodes,
  degree_normal = deg_normal[all_nodes],
  degree_cancer = deg_cancer[all_nodes]
) %>%
  dplyr::mutate(
    delta_degree = degree_cancer - degree_normal,
    abs_rewiring = abs(delta_degree)
  ) %>%
  dplyr::arrange(desc(abs_rewiring))

head(rewiring_score, 15)
library(igraph)
library(dplyr)

# Kısa label fonksiyonu
short_iso_label <- function(x) {
  
  gene <- sub("_ENST.*$", "", x)
  enst <- sub("^.*ENST", "", x)
  
  # ENST00000372731 -> son 6 rakam
  enst_short <- substr(
    enst,
    nchar(enst) - 5,
    nchar(enst)
  )
  
  paste0(gene, "_", enst_short)
}


# ------------------------------------------------------------
# ORTAK NODE SETİ
# ------------------------------------------------------------

all_nodes <- union(
  unique(c(normal_edges$from, normal_edges$to)),
  unique(c(cancer_edges$from, cancer_edges$to))
)


gN <- igraph::graph_from_data_frame(
  normal_edges[, c("from","to")],
  directed = FALSE,
  vertices = all_nodes
)

gC <- igraph::graph_from_data_frame(
  cancer_edges[, c("from","to")],
  directed = FALSE,
  vertices = all_nodes
)


# ------------------------------------------------------------
# ORTAK LAYOUT
# ------------------------------------------------------------

union_edges <- unique(
  rbind(
    normal_edges[, c("from","to")],
    cancer_edges[, c("from","to")]
  )
)

g_union <- igraph::graph_from_data_frame(
  union_edges,
  directed = FALSE,
  vertices = all_nodes
)

set.seed(123)

coords <- igraph::layout_with_kk(g_union)

rownames(coords) <- V(g_union)$name


# ------------------------------------------------------------
# LABEL
# ------------------------------------------------------------

V(gN)$label <- short_iso_label(V(gN)$name)
V(gC)$label <- short_iso_label(V(gC)$name)


# degree'ye göre node size
V(gN)$size <- 6 + sqrt(igraph::degree(gN)) * 3
V(gC)$size <- 6 + sqrt(igraph::degree(gC)) * 3


# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

par(
  mfrow = c(1,2),
  mar = c(1,1,4,1)
)

plot(
  gN,
  layout = coords[V(gN)$name, ],
  vertex.size = V(gN)$size,
  vertex.label = V(gN)$label,
  vertex.label.cex = 0.55,
  vertex.label.dist = 0.8,
  edge.width = 1.5,
  edge.arrow.mode = 0,
  main = "Normal Bayesian Network"
)

plot(
  gC,
  layout = coords[V(gC)$name, ],
  vertex.size = V(gC)$size,
  vertex.label = V(gC)$label,
  vertex.label.cex = 0.55,
  vertex.label.dist = 0.8,
  edge.width = 1.5,
  edge.arrow.mode = 0,
  main = "Cancer Bayesian Network"
)

par(mfrow = c(1,1))
library(ggplot2)
library(dplyr)
library(tidyr)

top_rewired <- rewiring_score %>%
  dplyr::slice_max(
    order_by = abs_rewiring,
    n = 15,
    with_ties = FALSE
  ) %>%
  dplyr::mutate(
    label = short_iso_label(node)
  )


rewire_long <- top_rewired %>%
  dplyr::select(
    label,
    degree_normal,
    degree_cancer
  ) %>%
  tidyr::pivot_longer(
    cols = c(degree_normal, degree_cancer),
    names_to = "condition",
    values_to = "degree"
  ) %>%
  dplyr::mutate(
    condition = dplyr::recode(
      condition,
      degree_normal = "Normal",
      degree_cancer = "Cancer"
    )
  )


ggplot() +
  
  # Normal-Cancer arasındaki çizgi
  geom_segment(
    data = top_rewired,
    aes(
      x = degree_normal,
      xend = degree_cancer,
      y = reorder(label, abs_rewiring),
      yend = reorder(label, abs_rewiring)
    ),
    linewidth = 1,
    alpha = 0.5
  ) +
  
  # Noktalar
  geom_point(
    data = rewire_long,
    aes(
      x = degree,
      y = reorder(label, top_rewired$abs_rewiring),
      shape = condition
    ),
    size = 4
  ) +
  
  scale_x_continuous(
    breaks = 0:6
  ) +
  
  labs(
    title = "Top Rewired Isoforms",
    subtitle = "Bayesian network connectivity: Normal vs Cancer",
    x = "Number of stable Bayesian network connections",
    y = NULL,
    shape = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    legend.position = "top"
  )
# Sıralamayı ÖNCEDEN belirle
label_order <- top_rewired %>%
  dplyr::arrange(abs_rewiring) %>%
  dplyr::pull(label)

top_rewired <- top_rewired %>%
  dplyr::mutate(
    label = factor(label, levels = label_order)
  )

rewire_long <- top_rewired %>%
  dplyr::select(
    label,
    degree_normal,
    degree_cancer
  ) %>%
  tidyr::pivot_longer(
    cols = c(degree_normal, degree_cancer),
    names_to = "condition",
    values_to = "degree"
  ) %>%
  dplyr::mutate(
    condition = dplyr::recode(
      condition,
      degree_normal = "Normal",
      degree_cancer = "Cancer"
    )
  )


ggplot() +
  geom_segment(
    data = top_rewired,
    aes(
      x = degree_normal,
      xend = degree_cancer,
      y = label,
      yend = label
    ),
    linewidth = 1,
    alpha = 0.5
  ) +
  
  geom_point(
    data = rewire_long,
    aes(
      x = degree,
      y = label,
      shape = condition
    ),
    size = 4
  ) +
  
  scale_x_continuous(
    breaks = 0:6
  ) +
  
  labs(
    title = "Top Rewired Isoforms",
    subtitle = "Bayesian network connectivity: Normal vs Cancer",
    x = "Number of stable Bayesian network connections",
    y = NULL,
    shape = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(face = "bold", size = 18),
    axis.text.y = element_text(face = "bold"),
    legend.position = "top"
  )
all_rewired <- rewiring_score %>%
  dplyr::mutate(
    label = short_iso_label(node)
  ) %>%
  dplyr::arrange(abs_rewiring)

# Sıralama
all_rewired$label <- factor(
  all_rewired$label,
  levels = all_rewired$label
)

all_rewire_long <- all_rewired %>%
  dplyr::select(
    label,
    degree_normal,
    degree_cancer
  ) %>%
  tidyr::pivot_longer(
    cols = c(degree_normal, degree_cancer),
    names_to = "condition",
    values_to = "degree"
  ) %>%
  dplyr::mutate(
    condition = dplyr::recode(
      condition,
      degree_normal = "Normal",
      degree_cancer = "Cancer"
    )
  )

ggplot() +
  geom_segment(
    data = all_rewired,
    aes(
      x = degree_normal,
      xend = degree_cancer,
      y = label,
      yend = label
    ),
    linewidth = 0.7,
    alpha = 0.5
  ) +
  
  geom_point(
    data = all_rewire_long,
    aes(
      x = degree,
      y = label,
      shape = condition
    ),
    size = 3
  ) +
  
  scale_x_continuous(
    breaks = 0:max(all_rewired$degree_normal,
                   all_rewired$degree_cancer)
  ) +
  
  labs(
    title = "Isoform Network Rewiring",
    subtitle = "All isoforms: Normal vs Cancer Bayesian networks",
    x = "Number of stable Bayesian network connections",
    y = NULL,
    shape = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(face = "bold", size = 18),
    axis.text.y = element_text(size = 8),
    legend.position = "top"
  )
library(igraph)
library(dplyr)

# ============================================================
# 1. CANCER STABLE EDGES + BOOTSTRAP STRENGTH
# ============================================================

cancer_strength <- as.data.frame(boot_cancer) %>%
  dplyr::select(from, to, strength, direction)

cancer_plot_edges <- cancer_edges %>%
  dplyr::select(from, to) %>%
  dplyr::left_join(
    cancer_strength,
    by = c("from", "to")
  )

# Normal'de edge var mı?
normal_edge_ids <- normal_edges$edge

cancer_plot_edges <- cancer_plot_edges %>%
  dplyr::mutate(
    edge_id = paste(
      pmin(from, to),
      pmax(from, to),
      sep = " -- "
    ),
    edge_type = ifelse(
      edge_id %in% normal_edge_ids,
      "Common",
      "Cancer-specific"
    )
  )


# ============================================================
# 2. NETWORK
# ============================================================

g <- igraph::graph_from_data_frame(
  cancer_plot_edges,
  directed = FALSE
)


# ============================================================
# 3. GENE INFORMATION
# ============================================================

V(g)$gene <- sub(
  "\\..*$",
  "",
  V(g)$name
)

# Eğer isimlerde "_" kullanılıyorsa
V(g)$gene <- sub(
  "_ENST.*$",
  "",
  V(g)$name
)


# ============================================================
# 4. SHORT ISOFORM LABEL
# ============================================================

short_label <- function(x) {

  gene <- sub("_ENST.*$", "", x)

  enst <- sub("^.*ENST", "", x)

  last6 <- substr(
    enst,
    nchar(enst) - 5,
    nchar(enst)
  )

  paste0(gene, "_", last6)
}

V(g)$label <- short_label(V(g)$name)


# ============================================================
# 5. NODE SIZE = DEGREE
# ============================================================

deg <- igraph::degree(g)

V(g)$size <- 6 + deg * 2.5


# ============================================================
# 6. GENE-GROUPED ORDER
#
# Aynı gene ait isoformlar yan yana gelecek
# ============================================================

node_order <- order(
  V(g)$gene,
  -deg
)

g <- igraph::permute(
  g,
  node_order
)


# ============================================================
# 7. EDGE STYLE
# ============================================================

E(g)$width <- 1 + 4 * E(g)$strength

E(g)$color <- ifelse(
  E(g)$edge_type == "Cancer-specific",
  "firebrick",
  "grey70"
)


# ============================================================
# 8. NODE STYLE
#
# Önceki analizde önemli çıkan genleri vurgula
# ============================================================

key_genes <- c(
  "ADD3",
  "SPTAN1",
  "DMTN",
  "ACTN1",
  "DIAPH1",
  "SLC2A1"
)

V(g)$color <- ifelse(
  V(g)$gene %in% key_genes,
  "goldenrod1",
  "grey90"
)

V(g)$frame.color <- "grey30"


# ============================================================
# 9. CIRCULAR PLOT
# ============================================================

lay <- igraph::layout_in_circle(g)

plot(
  g,
  layout = lay,

  vertex.size = V(g)$size,

  vertex.label = V(g)$label,
  vertex.label.cex = 0.55,
  vertex.label.dist = 1.2,

  vertex.color = V(g)$color,
  vertex.frame.color = V(g)$frame.color,

  edge.color = E(g)$color,
  edge.width = E(g)$width,

  edge.curved = 0.15,

  main = "Cancer Isoform Bayesian Network"
)

legend(
  "topleft",

  legend = c(
    "Cancer-specific edge",
    "Shared with Normal",
    "Key isoform genes"
  ),

  col = c(
    "firebrick",
    "grey70",
    "goldenrod1"
  ),

  lwd = c(3, 3, NA),
  pch = c(NA, NA, 21),
  pt.bg = c(NA, NA, "goldenrod1"),

  bty = "n"
)
# Gerekirse:
# install.packages("circlize")

library(circlize)
library(dplyr)

# ============================================================
# 1. CANCER EDGE TABLOSU
# ============================================================

circ_edges <- as.data.frame(boot_cancer) %>%
  dplyr::filter(strength >= 0.70) %>%
  dplyr::mutate(
    from_gene = sub("\\..*$", "", from),
    from_gene = sub("_ENST.*$", "", from_gene),
    
    to_gene = sub("\\..*$", "", to),
    to_gene = sub("_ENST.*$", "", to_gene)
  )

# Sadece averaged network'te gerçekten kalan edge'ler
stable_ids <- cancer_edges$edge

circ_edges <- circ_edges %>%
  dplyr::mutate(
    edge = paste(
      pmin(from, to),
      pmax(from, to),
      sep = " -- "
    )
  ) %>%
  dplyr::filter(edge %in% stable_ids)


# ============================================================
# 2. KISA ISOFORM İSİMLERİ
# ============================================================

short_name <- function(x) {
  
  gene <- sub("_ENST.*$", "", x)
  
  enst <- sub("^.*ENST", "", x)
  
  last6 <- substr(
    enst,
    nchar(enst) - 5,
    nchar(enst)
  )
  
  paste0(gene, "_", last6)
}

circ_edges$from_short <- short_name(circ_edges$from)
circ_edges$to_short   <- short_name(circ_edges$to)


# ============================================================
# 3. CHORD INPUT
# strength = chord kalınlığı
# ============================================================

chord_df <- circ_edges %>%
  dplyr::select(
    from_short,
    to_short,
    strength
  )


# ============================================================
# 4. CIRCOS / CHORD DIAGRAM
# ============================================================

circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = 2,
  track.margin = c(0.01, 0.01)
)

chordDiagram(
  chord_df,
  directional = 0,
  transparency = 0.55,
  annotationTrack = "grid",
  preAllocateTracks = 1
)


# ============================================================
# 5. LABEL'LAR
# ============================================================

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    
    sector.name <- CELL_META$sector.index
    
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.55
    )
  },
  bg.border = NA
)

title(
  "Cancer Isoform Bayesian Network",
  cex.main = 1.5,
  font.main = 2
)
library(circlize)
library(dplyr)

# ============================================================
# 1. STABLE CANCER EDGES
# ============================================================

circ_edges <- as.data.frame(boot_cancer) %>%
  dplyr::filter(strength >= 0.70) %>%
  dplyr::mutate(
    edge = paste(
      pmin(from, to),
      pmax(from, to),
      sep = " -- "
    )
  ) %>%
  dplyr::filter(edge %in% cancer_edges$edge)


# ============================================================
# 2. COMMON vs CANCER-SPECIFIC
# ============================================================

circ_edges <- circ_edges %>%
  dplyr::mutate(
    edge_type = ifelse(
      edge %in% common_edges,
      "Common",
      "Cancer-specific"
    )
  )


# ============================================================
# 3. SHORT LABEL
# ============================================================

short_name <- function(x) {
  
  gene <- sub("_ENST.*$", "", x)
  enst <- sub("^.*ENST", "", x)
  
  last6 <- substr(
    enst,
    nchar(enst) - 5,
    nchar(enst)
  )
  
  paste0(gene, "_", last6)
}

circ_edges$from_short <- short_name(circ_edges$from)
circ_edges$to_short   <- short_name(circ_edges$to)


# ============================================================
# 4. NODE INFORMATION
# ============================================================

nodes <- unique(
  c(circ_edges$from_short,
    circ_edges$to_short)
)

node_info <- data.frame(
  node = nodes
) %>%
  dplyr::mutate(
    gene = sub("_.*$", "", node)
  ) %>%
  dplyr::arrange(gene, node)


# ============================================================
# 5. GENE COLORS
# ============================================================

genes <- unique(node_info$gene)

gene_cols <- grDevices::hcl.colors(
  length(genes),
  palette = "Dynamic"
)

names(gene_cols) <- genes

grid_cols <- gene_cols[node_info$gene]
names(grid_cols) <- node_info$node


# ============================================================
# 6. NODE ORDER
# aynı gene ait isoformlar yan yana
# ============================================================

node_order <- node_info$node


# ============================================================
# 7. EDGE COLORS
# ============================================================

link_cols <- ifelse(
  circ_edges$edge_type == "Cancer-specific",
  "#B2182B",
  "#BDBDBD"
)

# transparency
link_cols <- grDevices::adjustcolor(
  link_cols,
  alpha.f = 0.65
)


# ============================================================
# 8. CHORD DATA
# strength = bağlantı kalınlığı
# ============================================================

chord_df <- data.frame(
  from = circ_edges$from_short,
  to = circ_edges$to_short,
  value = circ_edges$strength
)


# ============================================================
# 9. PLOT
# ============================================================

circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = 1.5,
  track.margin = c(0.01, 0.01),
  points.overflow.warning = FALSE
)

chordDiagram(
  x = chord_df,
  
  order = node_order,
  
  grid.col = grid_cols,
  
  col = link_cols,
  
  transparency = 0,
  
  directional = 0,
  
  annotationTrack = "grid",
  
  preAllocateTracks = list(
    track.height = 0.12
  )
)


# ============================================================
# 10. ISOFORM LABELS
# ============================================================

circos.trackPlotRegion(
  track.index = 1,
  
  panel.fun = function(x, y) {
    
    sector <- CELL_META$sector.index
    
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      sector,
      
      facing = "clockwise",
      niceFacing = TRUE,
      
      adj = c(0, 0.5),
      cex = 0.52
    )
  },
  
  bg.border = NA
)


# ============================================================
# 11. TITLE
# ============================================================

title(
  "Cancer Isoform Bayesian Network",
  cex.main = 1.6,
  font.main = 2
)


# ============================================================
# 12. LEGEND — EDGE TYPE
# ============================================================

legend(
  "topleft",
  
  legend = c(
    "Cancer-specific",
    "Shared with Normal"
  ),
  
  col = c(
    "#B2182B",
    "#BDBDBD"
  ),
  
  lwd = 4,
  bty = "n",
  cex = 0.9
)
# ============================================================
# GENE GRUPLARI ARASINDA FARKLI BOŞLUK
# ============================================================

ordered_info <- node_info %>%
  dplyr::filter(node %in% node_order) %>%
  dplyr::arrange(match(node, node_order))

# Default: aynı gene ait isoformlar arasında çok küçük boşluk
gap_vector <- rep(0.7, nrow(ordered_info))

# Bir sonraki node farklı gene aitse büyük boşluk bırak
gene_end <- c(
  ordered_info$gene[-1] != ordered_info$gene[-nrow(ordered_info)],
  TRUE
)

gap_vector[gene_end] <- 5
circos.par(
  start.degree = 90,
  gap.degree = 1.5,
  track.margin = c(0.01, 0.01),
  points.overflow.warning = FALSE
)
circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = gap_vector,
  track.margin = c(0.01, 0.01),
  points.overflow.warning = FALSE
)
chordDiagram(
  x = chord_df,
  order = node_order,
  grid.col = grid_cols,
  col = link_cols,
  transparency = 0,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = list(
    track.height = 0.12
  )
)
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {

    sector <- CELL_META$sector.index

    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      sector,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.52
    )
  },
  bg.border = NA
)

title(
  "Cancer Isoform Bayesian Network",
  cex.main = 1.6,
  font.main = 2
)
# ============================================================
# 1. EDGE'LERİ 4 GRUBA AYIR
# ============================================================

circ_edges <- circ_edges %>%
  dplyr::mutate(
    
    from_gene = sub("_ENST.*$", "", from),
    to_gene   = sub("_ENST.*$", "", to),
    
    gene_relation = ifelse(
      from_gene == to_gene,
      "Intra-gene",
      "Inter-gene"
    ),
    
    edge_class = dplyr::case_when(
      
      edge_type == "Cancer-specific" &
        gene_relation == "Intra-gene" ~
        "Cancer-specific | Intra-gene",
      
      edge_type == "Cancer-specific" &
        gene_relation == "Inter-gene" ~
        "Cancer-specific | Inter-gene",
      
      edge_type == "Common" &
        gene_relation == "Intra-gene" ~
        "Shared | Intra-gene",
      
      TRUE ~
        "Shared | Inter-gene"
    )
  )
# ============================================================
# 2. EDGE COLOR PALETTE
# ============================================================

edge_palette <- c(
  "Cancer-specific | Intra-gene" = "#8B0015",
  "Cancer-specific | Inter-gene" = "#D95F70",
  
  "Shared | Intra-gene" = "#4D4D4D",
  "Shared | Inter-gene" = "#C7C7C7"
)

link_cols <- edge_palette[circ_edges$edge_class]

link_cols <- grDevices::adjustcolor(
  link_cols,
  alpha.f = 0.75
)
chord_df <- data.frame(
  from = circ_edges$from_short,
  to = circ_edges$to_short,
  value = circ_edges$strength
)
circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = gap_vector,
  track.margin = c(0.01, 0.01),
  points.overflow.warning = FALSE
)

chordDiagram(
  x = chord_df,
  order = node_order,
  grid.col = grid_cols,
  col = link_cols,
  transparency = 0,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = list(
    track.height = 0.12
  )
)
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    
    sector <- CELL_META$sector.index
    
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      sector,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.52
    )
  },
  bg.border = NA
)

title(
  "Cancer Isoform Bayesian Network",
  cex.main = 1.6,
  font.main = 2
)
legend(
  "topleft",
  legend = c(
    "Cancer-specific: intra-gene",
    "Cancer-specific: inter-gene",
    "Shared: intra-gene",
    "Shared: inter-gene"
  ),
  col = c(
    "#8B0015",
    "#D95F70",
    "#4D4D4D",
    "#C7C7C7"
  ),
  lwd = 4,
  bty = "n",
  cex = 0.85
)
library(circlize)
library(dplyr)

# ============================================================
# 1. HER GENE İÇİN ANA RENK
# ============================================================

genes <- unique(node_info$gene)

gene_base_cols <- grDevices::hcl.colors(
  length(genes),
  palette = "Dark 3"
)

names(gene_base_cols) <- genes


# ============================================================
# 2. AYNI GENE AİT ISOFORMLARA AYNI RENK AİLESİNDEN TONLAR
# ============================================================

make_isoform_colors <- function(gene_name) {
  
  iso <- node_info$node[node_info$gene == gene_name]
  n <- length(iso)
  
  base_col <- gene_base_cols[gene_name]
  
  # Ana renkten açık tonlara
  pal <- grDevices::colorRampPalette(
    c(
      grDevices::adjustcolor(base_col, alpha.f = 0.45),
      base_col
    )
  )(max(n, 2))
  
  pal <- pal[seq_len(n)]
  names(pal) <- iso
  
  pal
}

isoform_cols <- unlist(
  lapply(genes, make_isoform_colors)
)

# node_order ile aynı sıraya getir
grid_cols <- isoform_cols[node_order]
names(grid_cols) <- node_order


# ============================================================
# 3. İKİ RENGİ KARIŞTIRAN FONKSİYON
# ============================================================

mix_colors <- function(col1, col2, alpha = 1) {
  
  rgb1 <- grDevices::col2rgb(col1)
  rgb2 <- grDevices::col2rgb(col2)
  
  mixed <- (rgb1 + rgb2) / 2
  
  grDevices::rgb(
    mixed[1],
    mixed[2],
    mixed[3],
    maxColorValue = 255,
    alpha = alpha * 255
  )
}


# ============================================================
# 4. HER EDGE İÇİN ISOFORM RENKLERİNİ BUL
# ============================================================

circ_edges <- circ_edges %>%
  dplyr::mutate(
    
    from_col = isoform_cols[from_short],
    to_col   = isoform_cols[to_short],
    
    from_gene = sub("_.*$", "", from_short),
    to_gene   = sub("_.*$", "", to_short),
    
    gene_relation = ifelse(
      from_gene == to_gene,
      "Intra-gene",
      "Inter-gene"
    )
  )


# ============================================================
# 5. CHORD RENGİ
#
# intra-gene  -> aynı renk ailesi
# inter-gene  -> iki gene/isoform renginin karışımı
#
# cancer-specific -> güçlü/doygun
# shared          -> transparan/soluk
# ============================================================

link_cols <- mapply(
  
  function(c1, c2, edge_type, relation) {
    
    if (relation == "Intra-gene") {
      
      # Aynı gene ait olduğundan renk ailesi korunur
      mixed <- mix_colors(c1, c2, alpha = 1)
      
    } else {
      
      # Farklı genler: iki rengin ortalaması
      mixed <- mix_colors(c1, c2, alpha = 1)
    }
    
    
    if (edge_type == "Cancer-specific") {
      
      # Cancer-specific: güçlü
      grDevices::adjustcolor(
        mixed,
        alpha.f = 0.85
      )
      
    } else {
      
      # Shared: daha soluk
      grDevices::adjustcolor(
        mixed,
        alpha.f = 0.25
      )
    }
    
  },
  
  circ_edges$from_col,
  circ_edges$to_col,
  circ_edges$edge_type,
  circ_edges$gene_relation
)


# ============================================================
# 6. CHORD DATA
# ============================================================

chord_df <- data.frame(
  from  = circ_edges$from_short,
  to    = circ_edges$to_short,
  value = circ_edges$strength
)


# ============================================================
# 7. PLOT
# ============================================================

circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = gap_vector,
  track.margin = c(0.01, 0.01),
  points.overflow.warning = FALSE
)

chordDiagram(
  x = chord_df,
  
  order = node_order,
  
  # dış halka = isoform rengi
  grid.col = grid_cols,
  
  # chord = iki isoform renginin karışımı
  col = link_cols,
  
  transparency = 0,
  
  directional = 0,
  
  annotationTrack = "grid",
  
  link.sort = TRUE,
  link.decreasing = FALSE,
  
  preAllocateTracks = list(
    track.height = 0.13
  )
)


# ============================================================
# 8. ISOFORM LABEL
# ============================================================

circos.trackPlotRegion(
  track.index = 1,
  
  panel.fun = function(x, y) {
    
    sector <- CELL_META$sector.index
    
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      sector,
      
      facing = "clockwise",
      niceFacing = TRUE,
      
      adj = c(0, 0.5),
      cex = 0.50
    )
  },
  
  bg.border = NA
)


# ============================================================
# 9. TITLE
# ============================================================

title(
  "Cancer Isoform Bayesian Network",
  cex.main = 1.6,
  font.main = 2
)


# ============================================================
# 10. LEGEND
# ============================================================

legend(
  "topleft",
  
  legend = c(
    "Cancer-specific relationship",
    "Shared with Normal"
  ),
  
  col = c(
    "grey25",
    "grey75"
  ),
  
  lwd = c(5, 5),
  bty = "n",
  cex = 0.9
)
library(circlize)
library(dplyr)

# ============================================================
# 1. STABLE CANCER EDGES
# ============================================================

circ_edges <- as.data.frame(boot_cancer) %>%
  dplyr::filter(strength >= 0.70) %>%
  dplyr::mutate(
    edge = paste(
      pmin(from, to),
      pmax(from, to),
      sep = " -- "
    )
  ) %>%
  dplyr::filter(edge %in% cancer_edges$edge)


# ============================================================
# 2. COMMON vs CANCER-SPECIFIC
# ============================================================

circ_edges <- circ_edges %>%
  dplyr::mutate(
    edge_type = ifelse(
      edge %in% common_edges,
      "Common",
      "Cancer-specific"
    )
  )


# ============================================================
# 3. SHORT ISOFORM NAMES
# ============================================================

short_name <- function(x) {
  
  gene <- sub("_ENST.*$", "", x)
  enst <- sub("^.*ENST", "", x)
  
  last6 <- substr(
    enst,
    nchar(enst) - 5,
    nchar(enst)
  )
  
  paste0(gene, "_", last6)
}

circ_edges$from_short <- short_name(circ_edges$from)
circ_edges$to_short   <- short_name(circ_edges$to)


# ============================================================
# 4. NODE INFORMATION
# ============================================================

nodes <- unique(
  c(
    circ_edges$from_short,
    circ_edges$to_short
  )
)

node_info <- data.frame(
  node = nodes
) %>%
  dplyr::mutate(
    gene = sub("_.*$", "", node)
  ) %>%
  dplyr::arrange(gene, node)

genes <- unique(node_info$gene)


# ============================================================
# 5. BASE COLOR FOR EACH GENE
# ============================================================

gene_base_cols <- grDevices::hcl.colors(
  length(genes),
  palette = "Dark 3"
)

names(gene_base_cols) <- genes


# ============================================================
# 6. ISOFORM COLORS
# Same gene = same color family
# Different isoforms = different shades
# ============================================================

make_isoform_colors <- function(gene_name) {
  
  iso <- node_info$node[
    node_info$gene == gene_name
  ]
  
  n <- length(iso)
  
  base_col <- gene_base_cols[gene_name]
  
  pal <- grDevices::colorRampPalette(
    c(
      grDevices::adjustcolor(
        base_col,
        alpha.f = 0.45
      ),
      base_col
    )
  )(max(n, 2))
  
  pal <- pal[seq_len(n)]
  names(pal) <- iso
  
  pal
}

isoform_cols <- unlist(
  lapply(
    genes,
    make_isoform_colors
  )
)


# ============================================================
# 7. NODE ORDER
# Same-gene isoforms stay together
# ============================================================

node_order <- node_info$node

grid_cols <- isoform_cols[node_order]
names(grid_cols) <- node_order


# ============================================================
# 8. GAPS
#
# Small gap = within gene
# Large gap = between genes
# ============================================================

ordered_info <- node_info %>%
  dplyr::filter(node %in% node_order) %>%
  dplyr::arrange(
    match(node, node_order)
  )

gap_vector <- rep(
  0.7,
  nrow(ordered_info)
)

gene_end <- c(
  ordered_info$gene[-1] !=
    ordered_info$gene[-nrow(ordered_info)],
  TRUE
)

gap_vector[gene_end] <- 5


# ============================================================
# 9. FUNCTION TO MIX TWO ISOFORM COLORS
# ============================================================

mix_colors <- function(
  col1,
  col2,
  alpha = 1
) {
  
  rgb1 <- grDevices::col2rgb(col1)
  rgb2 <- grDevices::col2rgb(col2)
  
  mixed <- (rgb1 + rgb2) / 2
  
  grDevices::rgb(
    mixed[1],
    mixed[2],
    mixed[3],
    maxColorValue = 255,
    alpha = alpha * 255
  )
}


# ============================================================
# 10. EDGE INFORMATION
# ============================================================

circ_edges <- circ_edges %>%
  dplyr::mutate(
    
    from_col = isoform_cols[from_short],
    to_col   = isoform_cols[to_short],
    
    from_gene = sub(
      "_.*$",
      "",
      from_short
    ),
    
    to_gene = sub(
      "_.*$",
      "",
      to_short
    ),
    
    gene_relation = ifelse(
      from_gene == to_gene,
      "Intra-gene",
      "Inter-gene"
    )
  )


# ============================================================
# 11. EDGE COLORS
#
# Color = connected isoform color families
#
# Dark/opaque = Cancer-specific
# Pale/transparent = Shared with Normal
# ============================================================

link_cols <- mapply(
  
  function(
    c1,
    c2,
    edge_type,
    relation
  ) {
    
    mixed <- mix_colors(
      c1,
      c2,
      alpha = 1
    )
    
    if (edge_type == "Cancer-specific") {
      
      grDevices::adjustcolor(
        mixed,
        alpha.f = 0.85
      )
      
    } else {
      
      grDevices::adjustcolor(
        mixed,
        alpha.f = 0.25
      )
    }
  },
  
  circ_edges$from_col,
  circ_edges$to_col,
  circ_edges$edge_type,
  circ_edges$gene_relation
)


# ============================================================
# 12. CHORD DATA
#
# Bootstrap strength controls chord width
# ============================================================

chord_df <- data.frame(
  from  = circ_edges$from_short,
  to    = circ_edges$to_short,
  value = circ_edges$strength
)


# ============================================================
# 13. DRAW CIRCOS
# ============================================================

circos.clear()

circos.par(
  start.degree = 90,
  gap.degree = gap_vector,
  track.margin = c(0.01, 0.01),
  points.overflow.warning = FALSE
)

chordDiagram(
  x = chord_df,
  
  order = node_order,
  
  # outer ring = isoform colors
  grid.col = grid_cols,
  
  # links = mixed isoform colors
  col = link_cols,
  
  transparency = 0,
  
  directional = 0,
  
  annotationTrack = "grid",
  
  link.sort = TRUE,
  link.decreasing = FALSE,
  
  preAllocateTracks = list(
    track.height = 0.13
  )
)


# ============================================================
# 14. ISOFORM LABELS
# ============================================================

circos.trackPlotRegion(
  track.index = 1,
  
  panel.fun = function(x, y) {
    
    sector <- CELL_META$sector.index
    
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      sector,
      
      facing = "clockwise",
      niceFacing = TRUE,
      
      adj = c(0, 0.5),
      cex = 0.50
    )
  },
  
  bg.border = NA
)


# ============================================================
# 15. TITLE
# ============================================================

title(
  "Cancer Isoform Bayesian Network",
  cex.main = 1.6,
  font.main = 2
)


# ============================================================
# 16. ONLY ONE LEGEND — TOP RIGHT
#
# Color itself = gene/isoform identity
# Opacity = relationship type
# ============================================================

legend(
  "topright",
  
  legend = c(
    "Cancer-specific relationship",
    "Shared with Normal"
  ),
  
  col = c(
    grDevices::adjustcolor(
      "black",
      alpha.f = 0.85
    ),
    grDevices::adjustcolor(
      "black",
      alpha.f = 0.25
    )
  ),
  
  lwd = 5,
  bty = "n",
  cex = 0.85,
  inset = c(0.02, 0.02)
)
library(dplyr)
library(igraph)

# Cancer network edges
inter_edges <- cancer_edges %>%
  mutate(
    from_gene = sub("_ENST.*$", "", from),
    to_gene   = sub("_ENST.*$", "", to)
  ) %>%
  filter(from_gene != to_gene)

# Undirected graph
g_inter <- graph_from_data_frame(
  inter_edges[, c("from", "to")],
  directed = FALSE
)

# Inter-gene degree
inter_hubs <- data.frame(
  node = names(degree(g_inter)),
  inter_gene_degree = degree(g_inter)
) %>%
  mutate(
    gene = sub("_ENST.*$", "", node)
  ) %>%
  arrange(desc(inter_gene_degree))

inter_hubs
deg_inter <- igraph::degree(g_inter)

inter_hubs <- data.frame(
  node = names(deg_inter),
  inter_gene_degree = as.numeric(deg_inter)
) %>%
  dplyr::mutate(
    gene = sub("_ENST.*$", "", node)
  ) %>%
  dplyr::arrange(desc(inter_gene_degree))

head(inter_hubs, 15)
library(ggplot2)
library(dplyr)

top_hubs <- inter_hubs %>%
  slice_head(n = 15) %>%
  mutate(
    label = paste0(
      gene, "_",
      substr(node, nchar(node) - 5, nchar(node))
    ),
    label = reorder(label, inter_gene_degree)
  )

ggplot(
  top_hubs,
  aes(
    x = inter_gene_degree,
    y = label
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = inter_gene_degree),
    hjust = -0.3,
    size = 4
  ) +
  labs(
    title = "Top Inter-Gene Isoform Hubs in Cancer",
    subtitle = "Intra-gene connections excluded",
    x = "Inter-gene degree",
    y = NULL
  ) +
  expand_limits(
    x = max(top_hubs$inter_gene_degree) + 1
  ) +
  theme_classic(base_size = 13)
all_hubs <- inter_hubs %>%
  filter(inter_gene_degree > 0) %>%
  mutate(
    label = paste0(
      gene, "_",
      substr(node, nchar(node) - 5, nchar(node))
    ),
    label = reorder(label, inter_gene_degree)
  )

ggplot(
  all_hubs,
  aes(
    x = inter_gene_degree,
    y = label
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = inter_gene_degree),
    hjust = -0.3,
    size = 3.5
  ) +
  labs(
    title = "Inter-Gene Isoform Hubs in Cancer",
    subtitle = "Intra-gene connections excluded",
    x = "Inter-gene degree",
    y = NULL
  ) +
  expand_limits(
    x = max(all_hubs$inter_gene_degree) + 1
  ) +
  theme_classic(base_size = 12)
library(igraph)

# İki networkteki bütün node'lar
all_nodes <- union(V(gN)$name, V(gC)$name)

# Aynı node setini kullan
gN2 <- add_vertices(
  gN,
  nv = length(setdiff(all_nodes, V(gN)$name)),
  name = setdiff(all_nodes, V(gN)$name)
)

gC2 <- add_vertices(
  gC,
  nv = length(setdiff(all_nodes, V(gC)$name)),
  name = setdiff(all_nodes, V(gC)$name)
)

# Aynı sıraya getir
gN2 <- permute(gN2, match(all_nodes, V(gN2)$name))
gC2 <- permute(gC2, match(all_nodes, V(gC2)$name))


# ------------------------------------------------------------
# ORTAK LAYOUT
# Normal + Cancer birleşik network üzerinden
# ------------------------------------------------------------

g_union <- union(gN2, gC2)

set.seed(123)

common_layout <- layout_with_fr(
  g_union,
  niter = 3000
)


# ------------------------------------------------------------
# KISA LABEL
# ------------------------------------------------------------

short_label <- function(x) {
  
  gene <- sub("_ENST.*$", "", x)
  enst <- sub("^.*ENST", "", x)
  
  paste0(
    gene, "_",
    substr(enst, nchar(enst)-5, nchar(enst))
  )
}

labels <- short_label(all_nodes)


# ------------------------------------------------------------
# NODE COLORS BY GENE
# ------------------------------------------------------------

genes <- sub("_ENST.*$", "", all_nodes)

gene_levels <- unique(genes)

gene_cols <- hcl.colors(
  length(gene_levels),
  "Dark 3"
)

names(gene_cols) <- gene_levels

node_cols <- gene_cols[genes]


# ------------------------------------------------------------
# SIDE-BY-SIDE
# ------------------------------------------------------------

par(
  mfrow = c(1,2),
  mar = c(1,1,4,1)
)


# NORMAL
plot(
  gN2,
  layout = common_layout,
  
  vertex.color = node_cols,
  vertex.frame.color = "white",
  vertex.size = 12,
  
  vertex.label = labels,
  vertex.label.cex = 0.55,
  vertex.label.color = "black",
  vertex.label.dist = 0.7,
  
  edge.color = adjustcolor("grey40", 0.45),
  edge.width = 1.3,
  
  main = "Normal Bayesian Network"
)


# CANCER
plot(
  gC2,
  layout = common_layout,
  
  vertex.color = node_cols,
  vertex.frame.color = "white",
  vertex.size = 12,
  
  vertex.label = labels,
  vertex.label.cex = 0.55,
  vertex.label.color = "black",
  vertex.label.dist = 0.7,
  
  edge.color = adjustcolor("grey40", 0.45),
  edge.width = 1.3,
  
  main = "Cancer Bayesian Network"
)

par(mfrow = c(1,1))
set.seed(123)

common_layout <- layout_with_fr(
  g_union,
  niter = 5000,
  repulserad = vcount(g_union)^3
)

# Networkü biraz daha genişlet
common_layout <- common_layout * 1.35
par(
  mfrow = c(1,2),
  mar = c(1,1,4,1)
)

# =========================
# NORMAL
# =========================

plot(
  gN2,
  layout = common_layout,

  vertex.color = node_cols,
  vertex.frame.color = "white",

  vertex.size = 9,

  vertex.label = labels,
  vertex.label.cex = 0.42,
  vertex.label.color = "grey15",
  vertex.label.dist = 1.1,

  edge.color = adjustcolor("grey40", alpha.f = 0.35),
  edge.width = 1,

  main = "Normal Bayesian Network"
)


# =========================
# CANCER
# =========================

plot(
  gC2,
  layout = common_layout,

  vertex.color = node_cols,
  vertex.frame.color = "white",

  vertex.size = 9,

  vertex.label = labels,
  vertex.label.cex = 0.42,
  vertex.label.color = "grey15",
  vertex.label.dist = 1.1,

  edge.color = adjustcolor("grey40", alpha.f = 0.35),
  edge.width = 1,

  main = "Cancer Bayesian Network"
)

par(mfrow = c(1,1))
library(igraph)

# ============================================================
# 1. ALL NODES
# ============================================================

all_nodes <- union(
  V(gN)$name,
  V(gC)$name
)


# ============================================================
# 2. ADD MISSING NODES TO EACH NETWORK
# Böylece Normal ve Cancer aynı node setine sahip olur
# ============================================================

missing_N <- setdiff(all_nodes, V(gN)$name)
missing_C <- setdiff(all_nodes, V(gC)$name)

gN2 <- gN
gC2 <- gC

if (length(missing_N) > 0) {
  gN2 <- add_vertices(
    gN2,
    nv = length(missing_N),
    name = missing_N
  )
}

if (length(missing_C) > 0) {
  gC2 <- add_vertices(
    gC2,
    nv = length(missing_C),
    name = missing_C
  )
}


# ============================================================
# 3. SAME NODE ORDER
# ============================================================

gN2 <- permute(
  gN2,
  match(all_nodes, V(gN2)$name)
)

gC2 <- permute(
  gC2,
  match(all_nodes, V(gC2)$name)
)


# ============================================================
# 4. UNION NETWORK
# Layout bunun üzerinden hesaplanacak
# ============================================================

g_union <- union(
  gN2,
  gC2
)


# ============================================================
# 5. COMMON LAYOUT
# Aynı isoform iki panelde de aynı yerde olacak
# ============================================================

set.seed(123)

common_layout <- layout_with_fr(
  g_union,
  niter = 5000
)

# Biraz daha ferah hale getir
common_layout <- common_layout * 1.7


# ============================================================
# 6. SHORT ISOFORM LABELS
# Örn:
# SPTAN1_ENST00000372731 -> SPTAN1_372731
# ============================================================

short_label <- function(x) {

  gene <- sub(
    "_ENST.*$",
    "",
    x
  )

  enst <- sub(
    "^.*ENST",
    "",
    x
  )

  last6 <- substr(
    enst,
    nchar(enst) - 5,
    nchar(enst)
  )

  paste0(
    gene,
    "_",
    last6
  )
}

labels <- short_label(all_nodes)


# ============================================================
# 7. GENE INFORMATION
# ============================================================

genes <- sub(
  "_ENST.*$",
  "",
  all_nodes
)

gene_levels <- unique(genes)


# ============================================================
# 8. ONE COLOR FAMILY PER GENE
# ============================================================

gene_cols <- grDevices::hcl.colors(
  length(gene_levels),
  palette = "Dark 3"
)

names(gene_cols) <- gene_levels

node_cols <- gene_cols[genes]


# ============================================================
# 9. NODE DEGREE
# Hafif degree-based node size
# ============================================================

deg_N <- degree(gN2)
deg_C <- degree(gC2)

# Aynı node'un iki networkteki maksimum degree'si
max_degree <- pmax(
  deg_N,
  deg_C
)

# Node boyutu çok değişmesin
node_size <- 7 + (max_degree * 1.3)

# Çok büyük olmasını engelle
node_size[node_size > 14] <- 14


# ============================================================
# 10. PLOT SETTINGS
# ============================================================

par(
  mfrow = c(1, 2),
  mar = c(1, 1, 4, 1)
)


# ============================================================
# 11. NORMAL NETWORK
# ============================================================

plot(
  gN2,

  layout = common_layout,

  # Nodes
  vertex.color = node_cols,
  vertex.frame.color = "white",
  vertex.frame.width = 1,

  vertex.size = node_size,

  # Labels
  vertex.label = labels,
  vertex.label.cex = 0.42,
  vertex.label.color = "grey15",
  vertex.label.dist = 1.1,

  # Edges
  edge.color = adjustcolor(
    "grey35",
    alpha.f = 0.35
  ),

  edge.width = 1.2,

  # No arrows — structure-focused comparison
  edge.arrow.size = 0,

  main = "Normal Bayesian Network"
)


# ============================================================
# 12. CANCER NETWORK
# ============================================================

plot(
  gC2,

  layout = common_layout,

  # Nodes
  vertex.color = node_cols,
  vertex.frame.color = "white",
  vertex.frame.width = 1,

  vertex.size = node_size,

  # Labels
  vertex.label = labels,
  vertex.label.cex = 0.42,
  vertex.label.color = "grey15",
  vertex.label.dist = 1.1,

  # Edges
  edge.color = adjustcolor(
    "grey35",
    alpha.f = 0.35
  ),

  edge.width = 1.2,

  edge.arrow.size = 0,

  main = "Cancer Bayesian Network"
)


# ============================================================
# 13. RESET GRAPHICS
# ============================================================

par(
  mfrow = c(1, 1)
)
library(igraph)

# ============================================================
# 1. ALL NODES
# ============================================================

all_nodes <- union(
  V(gN)$name,
  V(gC)$name
)


# ============================================================
# 2. ADD MISSING NODES
# Normal ve Cancer aynı node setine sahip olsun
# ============================================================

missing_N <- setdiff(all_nodes, V(gN)$name)
missing_C <- setdiff(all_nodes, V(gC)$name)

gN2 <- gN
gC2 <- gC

if (length(missing_N) > 0) {
  gN2 <- add_vertices(
    gN2,
    nv = length(missing_N),
    name = missing_N
  )
}

if (length(missing_C) > 0) {
  gC2 <- add_vertices(
    gC2,
    nv = length(missing_C),
    name = missing_C
  )
}


# ============================================================
# 3. SAME NODE ORDER
# ============================================================

gN2 <- permute(
  gN2,
  match(all_nodes, V(gN2)$name)
)

gC2 <- permute(
  gC2,
  match(all_nodes, V(gC2)$name)
)


# ============================================================
# 4. UNION NETWORK
# Ortak layout bunun üzerinden hesaplanacak
# ============================================================

g_union <- union(
  gN2,
  gC2
)


# ============================================================
# 5. COMMON LAYOUT
# İki panelde aynı isoform aynı yerde
# ============================================================

set.seed(123)

common_layout <- layout_with_fr(
  g_union,
  niter = 5000
)

# Görsel olarak biraz aç
common_layout <- common_layout * 1.7


# ============================================================
# 6. SHORT LABELS
# Örnek:
# SPTAN1_ENST00000372731 -> SPTAN1_372731
# ============================================================

short_label <- function(x) {

  gene <- sub(
    "_ENST.*$",
    "",
    x
  )

  enst <- sub(
    "^.*ENST",
    "",
    x
  )

  last6 <- substr(
    enst,
    nchar(enst) - 5,
    nchar(enst)
  )

  paste0(
    gene,
    "_",
    last6
  )
}

labels <- short_label(all_nodes)


# ============================================================
# 7. GENE INFORMATION
# ============================================================

genes <- sub(
  "_ENST.*$",
  "",
  all_nodes
)

gene_levels <- unique(genes)


# ============================================================
# 8. COLOR BY GENE
# Aynı gene ait isoformlar aynı renkte
# ============================================================

gene_cols <- grDevices::hcl.colors(
  length(gene_levels),
  palette = "Dark 3"
)

names(gene_cols) <- gene_levels

node_cols <- gene_cols[genes]


# ============================================================
# 9. NODE DEGREE
# FIX: degree class -> numeric
# ============================================================

deg_N <- as.numeric(
  igraph::degree(gN2)
)

deg_C <- as.numeric(
  igraph::degree(gC2)
)

# Normal veya Cancer'daki maksimum bağlantı sayısı
max_degree <- pmax(
  deg_N,
  deg_C
)

# Degree arttıkça node biraz büyüsün
node_size <- 7 + (max_degree * 1.3)

# Aşırı büyük node oluşmasın
node_size <- pmin(
  node_size,
  14
)


# ============================================================
# 10. GLOBAL PLOT SETTINGS
# ============================================================

par(
  mfrow = c(1, 2),
  mar = c(1, 1, 4, 1)
)


# ============================================================
# 11. NORMAL BAYESIAN NETWORK
# ============================================================

plot(
  gN2,

  layout = common_layout,

  # -----------------
  # NODE
  # -----------------

  vertex.color = node_cols,

  vertex.frame.color = "white",

  vertex.size = node_size,

  # -----------------
  # LABEL
  # -----------------

  vertex.label = labels,

  vertex.label.cex = 0.42,

  vertex.label.color = "grey15",

  vertex.label.dist = 1.1,

  # -----------------
  # EDGE
  # -----------------

  edge.color = grDevices::adjustcolor(
    "grey35",
    alpha.f = 0.35
  ),

  edge.width = 1.2,

  # Directionı görselde göstermiyoruz
  edge.arrow.size = 0,

  # -----------------
  # TITLE
  # -----------------

  main = "Normal Bayesian Network"
)


# ============================================================
# 12. CANCER BAYESIAN NETWORK
# ============================================================

plot(
  gC2,

  layout = common_layout,

  # -----------------
  # NODE
  # -----------------

  vertex.color = node_cols,

  vertex.frame.color = "white",

  vertex.size = node_size,

  # -----------------
  # LABEL
  # -----------------

  vertex.label = labels,

  vertex.label.cex = 0.42,

  vertex.label.color = "grey15",

  vertex.label.dist = 1.1,

  # -----------------
  # EDGE
  # -----------------

  edge.color = grDevices::adjustcolor(
    "grey35",
    alpha.f = 0.35
  ),

  edge.width = 1.2,

  edge.arrow.size = 0,

  # -----------------
  # TITLE
  # -----------------

  main = "Cancer Bayesian Network"
)


# ============================================================
# 13. RESET
# ============================================================

par(
  mfrow = c(1, 1)
)
