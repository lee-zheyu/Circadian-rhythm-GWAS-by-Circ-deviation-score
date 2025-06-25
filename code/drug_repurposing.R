library("readxl")
library("writexl")
library('dbplyr')
library('ggplot2')
library('ggrepel')

drug_df_complete <- read_excel("Desktop/circ_gene/data/drugbank_lookup_complete.xlsx")

condition_groups_final2 <- read_excel("Desktop/circ_gene/data/condition_groups_final3.xlsx")

hit_groups <- vector('list', length = nrow(drug_df_complete))
names(hit_groups) <- drug_df_complete$drug_id
for (i in 1:length(hit_groups)) {
  conditions_indiv <- unlist(strsplit(drug_df_complete$condition[i],';'))
  if (sum(is.na(conditions_indiv))>0) {
    next
  }
  condition_groups <- c()
  for (j in 1:length(conditions_indiv)) {
    condition_group <- condition_groups_final2$group[condition_groups_final2$condition == 
                                                       conditions_indiv[j]]
    if (length(condition_group) > 0) {
      condition_groups <- c(condition_groups, unlist(strsplit(condition_group, ';')))
    }
  }
  hit_groups[[i]] <- condition_groups
}
hit_groups_pool <- unlist(hit_groups)
hit_groups_table <- data.frame(table(hit_groups_pool))

all_groups <- as.character(hit_groups_table$hit_groups_pool)
hit_drug_by_groups <- vector('list', length = length(all_groups))
names(hit_drug_by_groups) <- all_groups
for (i in 1:length(hit_groups)) {
  individual_group <- hit_groups[i]
  drug_id <- names(individual_group)
  groups <- unlist(individual_group)
  if (length(groups) == 0) {
    next
  }
  for (j in 1:length(groups)) {
    hit_drug_by_groups[[groups[j]]] <- c(hit_drug_by_groups[[groups[j]]], drug_id)
  }
}

condition_df <- c()
for (i in 1:length(hit_drug_by_groups)) {
  individual_group <- hit_drug_by_groups[i]
  group <- names(individual_group)
  freq <- length(unlist(individual_group))
  condition_df <- rbind(condition_df, 
                        data.frame(group=group, frequency=freq))
}
group_seq <- order(condition_df$frequency, decreasing = TRUE)
condition_df_ordered <- condition_df[group_seq,]
condition_df_clean <- condition_df_ordered[condition_df_ordered$group!='discard',]

indication_complete <- read_excel("Desktop/circ_gene/data/indication_final.xlsx")
all_groups_clean <- condition_df_clean$group
all_drug_by_groups <- vector('list', length = length(all_groups_clean))
names(all_drug_by_groups) <- all_groups_clean
for (i in 1:nrow(indication_complete)) {
  groups_indiv <- unlist(strsplit(indication_complete$indication_group[i],';'))
  drug_id_indiv <- unlist(strsplit(indication_complete$indication_drug_id[i],';'))
  for (j in 1:length(groups_indiv)) {
    all_drug_by_groups[[groups_indiv[j]]] <- union(all_drug_by_groups[[groups_indiv[j]]], drug_id_indiv)
  }
}

drug_df <- c()
for (i in 1:length(all_drug_by_groups)) {
  individual_group <- all_drug_by_groups[i]
  group <- names(individual_group)
  freq <- length(unlist(individual_group))
  drug_df <- rbind(drug_df, 
                   data.frame(group=group, frequency=freq))
}
drug_df_clean <- drug_df[drug_df$group!='discard',]

p_val_df <- c()
hit_tot <- sum(condition_df_clean$frequency)
all_tot <- sum(drug_df_clean$frequency)
for (i in 1:length(all_groups_clean)) {
  group <- all_groups_clean[i]

  hit_freq <- condition_df_clean$frequency[condition_df_clean$group == group]
  all_freq <- drug_df_clean$frequency[drug_df_clean$group == group]
  all_perc <- all_freq / all_tot

  hit_prop_test <- prop.test(x = c(hit_freq, all_freq) , n = c(hit_tot, all_tot), #p = all_perc,
                             alternative = 'greater', correct = TRUE)

  group_df <- data.frame(condition_group = group,
                         hit_freq = hit_freq, hit_perc = hit_prop_test$estimate[1],
                         all_freq = all_freq, all_perc = all_perc,
                         p_val = hit_prop_test$p.value)
  p_val_df <- rbind(p_val_df, group_df)
}

write_xlsx(p_val_df,
           "Desktop/circ_gene/data/drug_p_val.xlsx")

condition_perc <- round(condition_df_clean$frequency/sum(condition_df_clean$frequency), digits = 4)
condition_df2 <- cbind(condition_df_clean,condition_perc)
non_prominent_groups <- condition_df2$group[condition_df2$frequency<=40]
other_perc <- condition_df2$condition_perc[condition_df2$group=='Other Conditions']
other_freq <- condition_df2$frequency[condition_df2$group=='Other Conditions']
row_to_be_removed <- c()
for (i in 1:length(non_prominent_groups)) {
  row_ind <- match(non_prominent_groups[i], condition_df2$group)
  other_perc <- other_perc + condition_df2$condition_perc[row_ind]
  other_freq <- other_freq + condition_df2$frequency[row_ind]
  row_to_be_removed <- c(row_to_be_removed, row_ind)
}
condition_df2$condition_perc[condition_df2$group=='Other Conditions'] <- other_perc
condition_df2$frequency[condition_df2$group=='Other Conditions'] <- other_freq
condition_df2 <- condition_df2[-row_to_be_removed,]
orders <- order(condition_df2$condition_perc, decreasing = TRUE)
new_group_order <- c(orders[-1], orders[1])
condition_df2 <- condition_df2[new_group_order,]
group_names <- condition_df2$group 
save(condition_df2, file = 'Desktop/circ_gene/data/condition_df2.RData')

condition_df3 <- condition_df2 %>% 
  mutate(csum = rev(cumsum(rev(condition_perc))), 
         pos = condition_perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), condition_perc/2, pos))

condition_df3$group <- factor(condition_df3$group, levels = group_names)

ggplot(condition_df3, aes(x = "" , y = condition_perc, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = condition_df3,
                   aes(y = pos, label = paste0(condition_perc*100, "%")),
                   size = 3, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "group")) +
  theme_void() + theme(legend.position = 'none') 

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'condition_groups.png',
  plot = last_plot(),
  width = 1.55,
  height = 1.55,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
)

drug_perc <- round(drug_df_clean$frequency / sum(drug_df_clean$frequency), digits = 4)
drug_df2 <- cbind(drug_df_clean, drug_perc)

other_perc2 <- drug_df2$drug_perc[drug_df2$group=='Other Conditions']
other_freq2 <- drug_df2$frequency[drug_df2$group=='Other Conditions']
row_to_be_removed2 <- c()
for (i in 1:length(non_prominent_groups)) {
  row_ind <- match(non_prominent_groups[i], drug_df2$group)
  other_perc2 <- other_perc2 + drug_df2$drug_perc[row_ind]
  other_freq2 <- other_freq2 + drug_df2$frequency[row_ind]
  row_to_be_removed2 <- c(row_to_be_removed2, row_ind)
}
drug_df2$drug_perc[drug_df2$group=='Other Conditions'] <- other_perc2
drug_df2$frequency[drug_df2$group=='Other Conditions'] <- other_freq2
drug_df2 <- drug_df2[-row_to_be_removed,]
new_group_order2 <- vector('numeric', length = nrow(drug_df2))
for (i in 1:nrow(drug_df2)) {
  new_group_order2[i] <- match(group_names[i], drug_df2$group)
}
drug_df2 <- drug_df2[new_group_order2,]
save(drug_df2, file = 'Desktop/circ_gene/data/drug_df2.RData')
load('Desktop/circ_gene/data/drug_df2.RData')

drug_df3 <- drug_df2 %>% 
  mutate(csum = rev(cumsum(rev(drug_perc))), 
         pos = drug_perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), drug_perc/2, pos))

ggplot(drug_df3, aes(x = "" , y = drug_perc, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = drug_df3,
                   aes(y = pos, label = paste0(drug_perc*100, "%")),
                   size = 3.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "group")) +
  theme_void() + theme(legend.title = element_blank(), legend.text = element_text(size = 8))

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'condition_groups.png',
  plot = last_plot(),
  width = 3.5,
  height = 3.5,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
)

p_val_df2 <- c()
hit_tot <- sum(condition_df2$frequency)
all_tot <- sum(drug_df2$frequency)
for (i in 1:nrow(condition_df2)) {
  group <- condition_df2$group[i]
  
  hit_freq <- condition_df2$frequency[condition_df2$group == group]
  all_freq <- drug_df2$frequency[drug_df2$group == group]
  all_perc <- all_freq / all_tot
  
  hit_prop_test <- prop.test(x = c(hit_freq, all_freq) , n = c(hit_tot, all_tot), #p = all_perc, 
                             alternative = 'greater', correct = TRUE)
  
  group_df <- data.frame(condition_group = group, 
                         hit_freq = hit_freq, hit_perc = hit_prop_test$estimate[1],
                         all_freq = all_freq, all_perc = all_perc,
                         p_val = hit_prop_test$p.value)
  
  p_val_df2 <- rbind(p_val_df2, group_df)
}

