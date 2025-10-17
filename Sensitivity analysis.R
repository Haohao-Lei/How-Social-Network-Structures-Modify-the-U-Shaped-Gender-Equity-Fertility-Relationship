# ============================================================
# Full sensitivity + visualization R script
# Titles now show the explicit fertility equations;
# X-axis label simplified to "Average GE"
# ============================================================

library(igraph)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(patchwork)

# ----------------------------
# PARAMETERS
# ----------------------------
set.seed(123)
N <- 10000
N_half <- N / 2
sw_nei <- 2
sw_p <- 0.1
distance_threshold <- 1
n_reps <- 50
bin_width <- 0.5
ge_breaks <- seq(0, 10, by = bin_width)
ge_values_plot <- seq(0, 10, length.out = 301)

# ----------------------------
# MICRO-LEVEL FERTILITY FUNCTIONS
# ----------------------------
fert_funcs <- list(
  sym_sharp = function(ge) (-3/250)*ge^3 + (8/25)*ge^2 - (23/10)*ge + 6,
  sym_flat  = function(ge) (-2/250)*ge^3 + (5/25)*ge^2 - (10/10)*ge + 5.5,
  asym_lr   = function(ge) {
    base <- (-3/250)*ge^3 + (8/25)*ge^2 - (23/10)*ge + 6
    base + (-0.15)*(ge - 5)
  }
)

# explicit equation titles for each function
eq_titles <- list(
  sym_sharp = expression(f(GE)==-0.012*GE^3 + 0.32*GE^2 - 2.3*GE + 6),
  sym_flat  = expression(f(GE)==-0.008*GE^3 + 0.20*GE^2 - 1.0*GE + 5.5),
  asym_lr   = expression(f(GE)==-0.012*GE^3 + 0.32*GE^2 - 2.3*GE + 6 - 0.15*(GE-5))
)

# ----------------------------
# NETWORK-MATCHING HELPER
# ----------------------------
match_via_network <- function(g, agents_df, distance_threshold=1){
  Nlocal <- vcount(g)
  neighs_all <- neighborhood(g, order=distance_threshold, nodes=V(g))
  male_vertices <- which(agents_df$gender=="M")
  prop_female_neighbors <- sapply(male_vertices,function(v){
    neigh<-setdiff(neighs_all[[v]],v)
    if(length(neigh)==0)return(0)
    mean(agents_df$gender[neigh]=="F")
  })
  male_order <- sample(male_vertices,length(male_vertices))
  paired_females <- integer(0)
  pairs_list <- vector("list",length=0)
  for(m in male_order){
    neigh<-setdiff(neighs_all[[m]],m)
    if(length(neigh)==0)next
    cand_fem<-neigh[agents_df$gender[neigh]=="F"]
    if(length(cand_fem)==0)next
    cand_fem<-setdiff(cand_fem,paired_females)
    if(length(cand_fem)==0)next
    selected<-sample(cand_fem,1)
    paired_females<-c(paired_females,selected)
    pairs_list[[length(pairs_list)+1]]<-list(male=m,female=selected)
  }
  if(length(pairs_list)==0){
    pairs_df<-data.frame()
  }else{
    pairs_df<-do.call(rbind,lapply(pairs_list,function(x){
      data.frame(male_id=x$male,female_id=x$female,
                 male_ge=agents_df$ge[x$male],
                 female_ge=agents_df$ge[x$female])
    }))
  }
  list(pairs=pairs_df,
       mean_prop_female_neighbors=mean(prop_female_neighbors),
       matching_rate=ifelse(length(male_vertices)>0,nrow(pairs_df)/length(male_vertices),NA))
}

# ----------------------------
# PAIR COMPUTATION HELPERS
# ----------------------------
compute_ideal_pairs <- function(agents_males, agents_females, fert_fun){
  males<-agents_males[order(agents_males$ge),]
  females<-agents_females[order(agents_females$ge),]
  n_pairs<-min(nrow(males),nrow(females))
  df<-data.frame(male_ge=males$ge[1:n_pairs],female_ge=females$ge[1:n_pairs])
  df$avg_ge<-(df$male_ge+df$female_ge)/2
  df$fertility<-pmin(fert_fun(df$male_ge),fert_fun(df$female_ge))
  df
}

compute_buffered_pairs <- function(agents_all,g,fert_fun){
  res<-match_via_network(g,agents_all,distance_threshold)
  pairs<-res$pairs
  if(nrow(pairs)==0)return(list(pairs=pairs,diag=res))
  pairs$avg_ge<-(pairs$male_ge+pairs$female_ge)/2
  pairs$fertility<-pmin(fert_fun(pairs$male_ge),fert_fun(pairs$female_ge))
  list(pairs=pairs,diag=res)
}

bin_and_summarize <- function(df){
  if(nrow(df)==0)return(data.frame())
  df<-df%>%mutate(bin=cut(avg_ge,breaks=ge_breaks,include.lowest=TRUE,right=FALSE))
  df%>%group_by(bin)%>%
    summarise(n=n(),
              mean_fert=mean(fertility,na.rm=TRUE),
              sd_fert=sd(fertility,na.rm=TRUE),
              .groups="drop")%>%
    mutate(bin_center=as.numeric(sub("\\[(.*),(.*)\\)","\\1",bin))+bin_width/2)
}

# ----------------------------
# MAIN SIMULATION
# ----------------------------
all_results<-list()
for(fname in names(fert_funcs)){
  fert_fun<-fert_funcs[[fname]]
  rep_bins<-list(); diag_list<-list()
  for(rep in seq_len(n_reps)){
    males<-data.frame(id=1:N_half,gender="M",ge=runif(N_half,0,10))
    females<-data.frame(id=1:N_half,gender="F",ge=runif(N_half,0,10))
    ideal_df<-compute_ideal_pairs(males,females,fert_fun)
    agents_all<-data.frame(id=1:N,gender=sample(rep(c("M","F"),each=N_half)),ge=runif(N,0,10))
    g<-sample_smallworld(dim=1,size=N,nei=sw_nei,p=sw_p)
    V(g)$ge<-agents_all$ge; V(g)$gender<-agents_all$gender
    buff_res<-compute_buffered_pairs(agents_all,g,fert_fun)
    ideal_bins<-bin_and_summarize(ideal_df)%>%mutate(model="ideal",rep=rep)
    buffer_bins<-bin_and_summarize(buff_res$pairs)%>%mutate(model="buffered",rep=rep)
    rep_bins[[rep]]<-bind_rows(ideal_bins,buffer_bins)
    diag_list[[rep]]<-buff_res$diag
  }
  rep_bins_df<-bind_rows(rep_bins)
  agg<-rep_bins_df%>%filter(!is.na(bin))%>%
    group_by(model,bin_center)%>%
    summarise(mean_of_mean_fert=mean(mean_fert,na.rm=TRUE),
              sd_of_mean_fert=sd(mean_fert,na.rm=TRUE),
              .groups="drop")%>%
    mutate(fert_func=fname)
  all_results[[fname]]<-list(agg=agg,rep_bins=rep_bins_df,diags=diag_list)
}

# ----------------------------
# PLOTTING
# ----------------------------
plot_list <- list()
for (fname in names(all_results)) {
  agg <- all_results[[fname]]$agg
  fert_fun <- fert_funcs[[fname]]
  micro_df <- data.frame(ge = ge_values_plot, micro_fert = fert_fun(ge_values_plot))
  
  p <- ggplot() +
    geom_ribbon(
      data = agg,
      aes(
        x = bin_center,
        ymin = mean_of_mean_fert - sd_of_mean_fert,
        ymax = mean_of_mean_fert + sd_of_mean_fert,
        fill = model
      ),
      alpha = 0.15,
      color = NA
    ) +
    geom_line(
      data = agg,
      aes(x = bin_center, y = mean_of_mean_fert, color = model),
      size = 1
    ) +
    geom_line(
      data = micro_df,
      aes(x = ge, y = micro_fert),
      color = "black",
      linetype = "dashed",
      size = 0.9
    ) +
    scale_color_viridis_d(option = "D") +
    scale_fill_viridis_d(option = "D") +
    labs(
      title = eq_titles[[fname]],
      x = "Average GE",
      y = "Mean couple fertility"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold")
    )
  
  # For Figure 1 and 3: keep y-axis ticks but hide label
  if (fname %in% c("sym_sharp", "asym_lr")) {
    p <- p + theme(axis.title.y = element_blank())
  }
  
  plot_list[[fname]] <- p
}

combined_plot <- (plot_list[["sym_sharp"]] / plot_list[["sym_flat"]] / plot_list[["asym_lr"]]) +
  plot_annotation(title = "")

print(combined_plot)
ggsave("combined_three_equations_titles.png", plot = combined_plot, width = 8, height = 12, dpi = 300)

# ----------------------------
# SAVE OUTPUTS
# ----------------------------
saveRDS(all_results,file="sensitivity_results_equations.rds")
write.csv(do.call(rbind,lapply(names(all_results),function(fn){
  df<-all_results[[fn]]$agg; df$fert_func<-fn; df
})),file="sensitivity_agg_equations.csv",row.names=FALSE)

cat("Done. Figure saved as combined_three_equations_titles.png\n")
# ---------------- END ----------------