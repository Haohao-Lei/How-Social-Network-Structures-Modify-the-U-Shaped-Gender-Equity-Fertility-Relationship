# Full ABM + visualization script
# - Builds ideal and buffered small-world matching scenarios
# - Produces two visualizations:
#     (A) Three-line: Ideal, Proximity-by-Similarity, Proximity-by-Dissimilarity (with jittered dots)
#     (B) Two-line: Ideal vs Buffered (Buffered = Proximity-by-Similarity) (with jittered dots)
# - WARNING: N = 100000 can be memory- and time-intensive for graph neighborhood computations.
#            Run on a machine with sufficient RAM and increase R's memory limits if necessary.
#
# Requirements: igraph, ggplot2, viridis, dplyr
# Save as e.g. "abm_full_100k.R" and run: source("abm_full_100k.R")
#
# Author: assistant (code generated for user's ABM)

library(igraph)
library(ggplot2)
library(viridis)
library(dplyr)

# ----------------------------
# PARAMETERS
# ----------------------------
set.seed(123)
N <- 100000L                  # total agents (must be even) as requested
if (N %% 2L != 0L) stop("N must be even")
N_half <- N / 2L
distance_threshold <- 1L      # 1-hop neighborhood for buffered matching
sw_nei <- 2L                  # small-world nearest neighbors on ring
sw_p <- 0.10                  # rewiring probability

# For plotting: subsample points shown as dots to avoid huge plotting time/memory
plot_point_sample <- 20000L   # number of jittered points to show (per combined dataset)

# ----------------------------
# MICRO-LEVEL FERTILITY FUNCTION (U-shaped cubic)
# f(GE) = -0.012*GE^3 + 0.32*GE^2 - 2.3*GE + 6
# ----------------------------
fertility_intention <- function(ge) {
  -0.012 * ge^3 + 0.32 * ge^2 - 2.3 * ge + 6
}

# ----------------------------
# Helper: match_via_network
# male-proposing: each male selects uniformly at random one unmatched female
# from his 1-hop neighborhood (excluding self). Returns pairs and diagnostics.
# NOTE: computing neighborhoods for all vertices at N=100k may be expensive.
# We compute male neighborhoods only to reduce memory/time.
# ----------------------------
match_via_network <- function(g, agents_df, distance_threshold = 1L) {
  # safety checks
  stopifnot(vcount(g) == nrow(agents_df))
  male_vertices <- which(agents_df$gender == "M")
  female_vertices <- which(agents_df$gender == "F")
  # compute ego neighborhoods for male vertices only (includes self)
  neighs_males <- ego(g, order = distance_threshold, nodes = male_vertices, mode = "all")
  # mapping: neighs_males is a list aligned with male_vertices
  
  # diagnostics: proportion female neighbors for each male
  prop_female_neighbors <- vapply(seq_along(male_vertices), function(i) {
    v <- male_vertices[i]
    neigh <- setdiff(as.integer(neighs_males[[i]]), v)
    if (length(neigh) == 0L) return(0.0)
    mean(agents_df$gender[neigh] == "F")
  }, numeric(1))
  
  # perform matching
  male_order <- sample(male_vertices, length(male_vertices))
  paired_females <- integer(0)
  pairs_list <- vector("list", 0L)
  # we need quick lookup mapping from male vertex to its ego list index
  male_to_index <- match(male_order, male_vertices)
  for (k in seq_along(male_order)) {
    m <- male_order[k]
    idx <- male_to_index[k]
    neigh <- setdiff(as.integer(neighs_males[[idx]]), m)
    if (length(neigh) == 0L) next
    cand_fem <- neigh[agents_df$gender[neigh] == "F"]
    if (length(cand_fem) == 0L) next
    # exclude already paired females
    cand_fem <- setdiff(cand_fem, paired_females)
    if (length(cand_fem) == 0L) next
    selected <- sample(cand_fem, 1L)
    paired_females <- c(paired_females, selected)
    pairs_list[[length(pairs_list) + 1L]] <- list(male = m, female = selected)
  }
  
  if (length(pairs_list) == 0L) {
    pairs_df <- data.frame()
  } else {
    pairs_df <- do.call(rbind, lapply(pairs_list, function(x) {
      data.frame(male_id = x$male,
                 female_id = x$female,
                 male_ge = agents_df$ge[x$male],
                 female_ge = agents_df$ge[x$female],
                 stringsAsFactors = FALSE)
    }))
    pairs_df$avg_ge <- (pairs_df$male_ge + pairs_df$female_ge) / 2
    pairs_df$fertility <- pmin(fertility_intention(pairs_df$male_ge),
                               fertility_intention(pairs_df$female_ge))
  }
  
  diagnostics <- list(
    mean_prop_female_neighbors = mean(prop_female_neighbors),
    prop_female_neighbors_summary = summary(prop_female_neighbors),
    matching_rate = ifelse(length(male_vertices) > 0L, nrow(pairs_df) / length(male_vertices), NA_real_),
    num_pairs = ifelse(is.data.frame(pairs_df), nrow(pairs_df), 0L)
  )
  
  return(list(pairs = pairs_df, diagnostics = diagnostics))
}

# ----------------------------
# 1) Ideal matching scenario
# - as in manuscript: both sexes receive same sorted GE sequence and are paired one-to-one
# ----------------------------
# We draw a single GE vector and assign identical to males and females so each couple
# has exactly equal GE (ideal matching)
ge_seq <- sort(runif(N_half, 0, 10))
males_ideal <- data.frame(id = 1:N_half, gender = "M", ge = ge_seq, stringsAsFactors = FALSE)
females_ideal <- data.frame(id = 1:N_half, gender = "F", ge = ge_seq, stringsAsFactors = FALSE)

ideal_pairs <- data.frame(
  male_ge = males_ideal$ge,
  female_ge = females_ideal$ge,
  stringsAsFactors = FALSE
)
ideal_pairs$avg_ge <- (ideal_pairs$male_ge + ideal_pairs$female_ge) / 2
ideal_pairs$fertility <- pmin(fertility_intention(ideal_pairs$male_ge),
                              fertility_intention(ideal_pairs$female_ge))

# ----------------------------
# 2) Buffered: Proximity by Similarity
# - place agents on ring such that neighbors have similar GE (sorted GE)
# - but randomize sex assignment across vertices (to avoid sex blocking)
# ----------------------------
agents_caseA <- data.frame(
  id = 1L:N,
  gender = sample(rep(c("M", "F"), each = N_half)), # randomized sex order across vertices
  ge = sort(runif(N, 0, 10)),
  stringsAsFactors = FALSE
)
gA <- sample_smallworld(dim = 1, size = N, nei = sw_nei, p = sw_p)
stopifnot(vcount(gA) == nrow(agents_caseA))
V(gA)$ge <- agents_caseA$ge
V(gA)$gender <- agents_caseA$gender

res_A <- match_via_network(gA, agents_caseA, distance_threshold)
pairs_A <- res_A$pairs
diag_A <- res_A$diagnostics
message("Similarity scenario: matched pairs = ", diag_A$num_pairs,
        ", matching rate = ", signif(diag_A$matching_rate, 4),
        ", mean_prop_female_neighbors = ", signif(diag_A$mean_prop_female_neighbors, 4))

# ----------------------------
# 3) Buffered: Proximity by Dissimilarity
# - create permuted GE ordering to maximize local dissimilarity:
#   e.g. positions 1,N,2,N-1,3,N-2,...
# - randomize sex assignment across vertices
# ----------------------------
ge_all <- sort(runif(N, 0, 10))
# create permutation 1, N, 2, N-1, ...
half_seq <- seq_len(N / 2L)
perm <- as.vector(rbind(half_seq, N - half_seq + 1L))
ge_perm <- ge_all[perm]
agents_caseB <- data.frame(
  id = 1L:N,
  gender = sample(rep(c("M", "F"), each = N_half)), # randomize sexes
  ge = ge_perm,
  stringsAsFactors = FALSE
)
gB <- sample_smallworld(dim = 1, size = N, nei = sw_nei, p = sw_p)
stopifnot(vcount(gB) == nrow(agents_caseB))
V(gB)$ge <- agents_caseB$ge
V(gB)$gender <- agents_caseB$gender

res_B <- match_via_network(gB, agents_caseB, distance_threshold)
pairs_B <- res_B$pairs
diag_B <- res_B$diagnostics
message("Dissimilarity scenario: matched pairs = ", diag_B$num_pairs,
        ", matching rate = ", signif(diag_B$matching_rate, 4),
        ", mean_prop_female_neighbors = ", signif(diag_B$mean_prop_female_neighbors, 4))

# ----------------------------
# Aggregate data frame for plotting
# ----------------------------
# Convert scenario labels to friendly names
get_plot_df <- function(pairs_df, scenario_label) {
  if (is.null(pairs_df) || nrow(pairs_df) == 0L) {
    return(data.frame(scenario = character(0), avg_ge = numeric(0), fertility = numeric(0)))
  }
  data.frame(scenario = scenario_label,
             avg_ge = pairs_df$avg_ge,
             fertility = pairs_df$fertility,
             stringsAsFactors = FALSE)
}

plot_df_full <- bind_rows(
  data.frame(scenario = "Ideal", avg_ge = ideal_pairs$avg_ge, fertility = ideal_pairs$fertility, stringsAsFactors = FALSE),
  get_plot_df(pairs_A, "Proximity by Similarity"),
  get_plot_df(pairs_B, "Proximity by Dissimilarity")
)

# ----------------------------
# For dot clouds: sample a subset of pair-level points (to keep plotting feasible)
# We'll sample up to plot_point_sample points from combined buffered pairs (both A and B)
# ----------------------------
buffered_points_all <- bind_rows(
  get_plot_df(pairs_A, "Proximity by Similarity"),
  get_plot_df(pairs_B, "Proximity by Dissimilarity")
)

set.seed(456)
if (nrow(buffered_points_all) > plot_point_sample) {
  buffered_points_sample <- buffered_points_all %>% sample_n(plot_point_sample)
} else {
  buffered_points_sample <- buffered_points_all
}

# For the two-line figure (Ideal vs Buffered (Similarity)), sample buffered similarity points only
set.seed(789)
buffered_similarity_points <- get_plot_df(pairs_A, "Proximity by Similarity")
if (nrow(buffered_similarity_points) > plot_point_sample) {
  buffered_similarity_points_sample <- buffered_similarity_points %>% sample_n(plot_point_sample)
} else {
  buffered_similarity_points_sample <- buffered_similarity_points
}

# ----------------------------
# Colors and theme (match user's preference)
# ----------------------------
three_color_map <- c(
  "Ideal" = "#FF4D4D",                      # red for Ideal in three-line
  "Proximity by Dissimilarity" = "#00B050", # green
  "Proximity by Similarity" = "#0070C0"     # blue
)
two_color_map <- c(
  "Ideal" = "#0070C0",       # blue for Ideal in two-line figure
  "Buffered" = "#FF4D4D"     # red for Buffered line/dots (we map Similarity -> Buffered)
)

# LOESS smoothing span
loess_span <- 0.5

# ----------------------------
# FIGURE A: Three-line version (Ideal + Similarity + Dissimilarity), with jittered dots overlay
# ----------------------------
p_three <- ggplot() +
  # jittered dots for buffered points (both A & B), light alpha
  geom_point(data = buffered_points_sample, aes(x = avg_ge, y = fertility, color = scenario),
             alpha = 0.1, size = 0.6, position = position_jitter(width = 0.02, height = 0.02)) +
  # LOESS smooth lines per scenario
  stat_smooth(data = plot_df_full, aes(x = avg_ge, y = fertility, color = scenario),
              method = "loess", se = FALSE, span = loess_span, size = 1.0) +
  # Add ideal curve as emphasized line (we already included it above) - keep consistent color mapping
  scale_color_manual(values = three_color_map) +
  labs(x = "Average Gender Equity Score", y = "Fertility Outcome", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

print(p_three)

# Save figure A
ggsave("Figure_A_three_line.png", plot = p_three, width = 9, height = 5, dpi = 200)
message("Saved Figure_A_three_line.png")

# ----------------------------
# FIGURE B: Two-line version (Ideal vs Buffered)
# - Buffered will be the Proximity-by-Similarity case (renamed to 'Buffered')
# - Overlay jittered dots from similarity buffered sample (red)
# ----------------------------
# Prepare data: rename similarity scenario to "Buffered" for this plot
ideal_df_plot <- data.frame(scenario = "Ideal", avg_ge = ideal_pairs$avg_ge, fertility = ideal_pairs$fertility, stringsAsFactors = FALSE)
similarity_df_plot <- get_plot_df(pairs_A, "Buffered")  # rename

plot_df_two <- bind_rows(ideal_df_plot, similarity_df_plot)
buffered_similarity_points_sample <- buffered_similarity_points_sample %>% mutate(scenario = "Buffered")

p_two <- ggplot() +
  geom_point(data = buffered_similarity_points_sample, aes(x = avg_ge, y = fertility, color = scenario),
             alpha = 0.1, size = 0.6, position = position_jitter(width = 0.02, height = 0.02)) +
  stat_smooth(data = plot_df_two, aes(x = avg_ge, y = fertility, color = scenario),
              method = "loess", se = FALSE, span = loess_span, size = 1.0) +
  scale_color_manual(values = two_color_map) +
  labs(x = "Average Gender Equity Score", y = "Fertility Outcome", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

print(p_two)

# Save figure B
ggsave("Figure_B_two_line.png", plot = p_two, width = 9, height = 5, dpi = 200)
message("Saved Figure_B_two_line.png")

# ----------------------------
# OPTIONAL: Original network visualization block (unchanged in style)
# NOTE: plotting entire 100k network is not practical; subsample nodes for plotting.
# We'll produce two PNGs that sample up to 1000 vertices for readability.
# ----------------------------
plot_network_sample <- function(g, agents_df, sample_size = 1000L, title = "Network sample", file = NULL) {
  Nlocal <- vcount(g)
  s <- min(sample_size, Nlocal)
  set.seed(42)
  samp_vertices <- sort(sample(seq_len(Nlocal), s))
  sg <- induced_subgraph(g, vids = samp_vertices)
  ge_vals <- V(sg)$ge
  if (is.null(ge_vals)) ge_vals <- agents_df$ge[samp_vertices]
  col_idx <- as.numeric(cut(ge_vals, breaks = 100))
  vcols <- viridis(100)[col_idx]
  layout_circ <- layout_in_circle(sg)
  png(file, width = 900, height = 900, res = 150)
  par(mar = c(0,0,2,0))
  plot(sg,
       layout = layout_circ,
       vertex.color = adjustcolor(vcols, alpha.f = 0.8),
       vertex.size = 6,
       vertex.label = NA,
       vertex.frame.color = NA,
       main = title)
  dev.off()
}

# Save network samples
plot_network_sample(gA, agents_caseA, sample_size = 900L, title = "Case A: Proximity by Similarity (sample)", file = "network_caseA_sample.png")
plot_network_sample(gB, agents_caseB, sample_size = 900L, title = "Case B: Proximity by Dissimilarity (sample)", file = "network_caseB_sample.png")
message("Saved network sample PNGs (may be large).")

# ----------------------------
# Save pair-level csvs (optional, large files)
# WARNING: these files may be large if many pairs exist. You can comment out if not needed.
# ----------------------------
# write.csv(plot_df_full, "abm_pairlevel_plot_df_full.csv", row.names = FALSE)
# write.csv(get_plot_df(pairs_A, "Proximity by Similarity"), "abm_pairs_similarity.csv", row.names = FALSE)
# write.csv(get_plot_df(pairs_B, "Proximity by Dissimilarity"), "abm_pairs_dissimilarity.csv", row.names = FALSE)

# ----------------------------
# Print summary diagnostics
# ----------------------------
cat("\n--- Diagnostics summary ---\n")
cat("Similarity scenario: mean_prop_female_neighbors =", signif(diag_A$mean_prop_female_neighbors, 4),
    ", matching_rate =", signif(diag_A$matching_rate, 4), "\n")
cat("Dissimilarity scenario: mean_prop_female_neighbors =", signif(diag_B$mean_prop_female_neighbors, 4),
    ", matching_rate =", signif(diag_B$matching_rate, 4), "\n")

cat("\nDone. Figures produced:\n - Figure_A_three_line.png\n - Figure_B_two_line.png\n - network_caseA_sample.png\n - network_caseB_sample.png\n")