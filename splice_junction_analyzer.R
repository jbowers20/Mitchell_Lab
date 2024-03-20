library(ggseqlogo)
library(ggpubr)
library(dplyr)
library(stats)
library(stringr)
library(tidyr)

######################### Collects data output #########################
setwd("~/Documents/splice_junction_project/analytic_docs")

RDS_files <- list.files(pattern = "*.rds")
RDS_list <- lapply(RDS_files, readRDS)
splice_sites <- do.call(rbind, RDS_list)

acceptor_sites <- splice_sites[!is.na(splice_sites$acceptor),]
donor_sites <- splice_sites[!is.na(splice_sites$donor),]

acceptor_sites <- acceptor_sites[(nchar(acceptor_sites$acceptor) == 40),]
donor_sites <- donor_sites[(nchar(donor_sites$donor) == 40),]
######################################################################

# Methods used in code below (findSpread, donor_quantifier, acceptor_quantifier, p.signif
findSpread <- function(sequences){
  seq_length <- nchar(sequences[1])
  spread <- data.frame(position = numeric(), a = numeric(), g = numeric(), t = numeric(), c = numeric())
  
  for (position in 1:seq_length) {
    Ade <- 0
    Gua <- 0
    Thy <- 0
    Cyt <- 0
    N <- 0
    for (sequence in sequences) {
      if (substr(sequence, position, position) == "A") {
        Ade <- Ade + 1
      } else if (substr(sequence, position, position) == "G") {
        Gua <- Gua + 1
      }
      else if (substr(sequence, position, position) == "T") {
        Thy <- Thy + 1
      }
      else if (substr(sequence, position, position) == "C") {
        Cyt <- Cyt + 1
      } else {
        N <- N + 1
      }
      
    }
    total = Ade + Gua + Thy + Cyt + N
    spread[position,] <- c(position, Ade / total, Gua / total, Thy / total, Cyt / total, N / total)
    #spread[[position]] <- c(Ade / total, Gua / total, Thy / total, Cyt / total, N / total)
  }
  return( spread )
}

donor_quantifier <- function(spread, sequence){
  score <- 0
  
  ## Position 1 ##
  # Only award points if this is a G
  pos1 <- substr(sequence, 1, 1)
  
  if (pos1 == "G") {
    score <- score + 1
  }
  
  ## Position 2 ##
  # only award points if this is a T
  pos2 <- substr(sequence, 2, 2)
  
  if (pos2 == "T"){
    score <- score + 1
  }
  
  
  ## Position 3 ##
  # Award points based on percentage of that nucleotide at the given position
  pos_3_spread <- spread[donor_spread$position == 3,]
  pos3 <- substr(sequence, 3, 3)
  if (pos3 == "A") {
    score <- score + pos_3_spread$a
  } else if (pos3 == "G") {
    score <- score + pos_3_spread$g
  } else if (pos3 == "T") {
    score <- score + pos_3_spread$t
  } else if (pos3 == "C") {
    score <- score + pos_3_spread$c
  }
  
  ## Position 4 ##
  # Award points based on percentage of that nucleotide at the given position
  pos_4_spread <- spread[donor_spread$position == 4,]
  pos4 <- substr(sequence, 4, 4)
  if (pos4 == "A") {
    score <- score + pos_4_spread$a
  } else if (pos4 == "G") {
    score <- score + pos_4_spread$g
  } else if (pos4 == "T") {
    score <- score + pos_4_spread$t
  } else if (pos4 == "C") {
    score <- score + pos_4_spread$c
  }
  
  ## Position 5 ##
  # Award points based on percentage of that nucleotide at the given position
  pos_5_spread <- spread[donor_spread$position == 5,]
  pos5 <- substr(sequence, 5, 5)
  if (pos5 == "A") {
    score <- score + pos_5_spread$a
  } else if (pos5 == "G") {
    score <- score + pos_5_spread$g
  } else if (pos5 == "T") {
    score <- score + pos_5_spread$t
  } else if (pos5 == "C") {
    score <- score + pos_5_spread$c
  }
  
  ## Position 6 ##
  # Award points based on percentage of that nucleotide at the given position
  pos_6_spread <- spread[spread$position == 6,]
  pos6 <- substr(sequence, 6, 6)
  if (pos6 == "A") {
    score <- score + pos_6_spread$a
  } else if (pos6 == "G") {
    score <- score + pos_6_spread$g
  } else if (pos6 == "T") {
    score <- score + pos_6_spread$t
  } else if (pos6 == "C") {
    score <- score + pos_6_spread$c
  }
  
  return(score)
}

acceptor_quantifier <- function(spread, sequence){
  score <- 0
  
  ## Pos 40 ##
  if (substr(sequence, 40, 40) == "G") {
    score <- score + 1 
  }
  ## Pos 39 ##
  if (substr(sequence, 39, 39) == "A") {
    score <- score + 1
  }
  ## Pos 38 ##
  pos_38_spread <- spread[spread$position == 38,]
  pos38 <- substr(sequence, 38, 38)
  if (pos38 == "A") {
    score <- score + pos_38_spread$a
  } else if (pos38 == "G") {
    score <- score + pos_38_spread$g
  } else if (pos38 == "T") {
    score <- score + pos_38_spread$t
  } else if (pos38 == "C") {
    score <- score + pos_38_spread$c
  }
  
  ## Pos 37 ##
  pos_37_spread <- spread[spread$position == 37,]
  pos37 <- substr(sequence, 37, 37)
  if (pos38 == "A") {
    score <- score + pos_37_spread$a
  } else if (pos37 == "G") {
    score <- score + pos_37_spread$g
  } else if (pos37 == "T") {
    score <- score + pos_37_spread$t
  } else if (pos37 == "C") {
    score <- score + pos_37_spread$c
  }
  
  ## Pos 36 ##
  pos_36_spread <- spread[spread$position == 36,]
  pos36 <- substr(sequence, 36, 36)
  if (pos36 == "A") {
    score <- score + pos_36_spread$a
  } else if (pos36 == "G") {
    score <- score + pos_36_spread$g
  } else if (pos36 == "T") {
    score <- score + pos_36_spread$t
  } else if (pos36 == "C") {
    score <- score + pos_36_spread$c
  }
  
  ## Pos 35 ##
  pos_35_spread <- spread[spread$position == 35,]
  pos35 <- substr(sequence, 35, 35)
  if (pos35 == "A") {
    score <- score + pos_35_spread$a
  } else if (pos35 == "G") {
    score <- score + pos_35_spread$g
  } else if (pos35 == "T") {
    score <- score + pos_35_spread$t
  } else if (pos35 == "C") {
    score <- score + pos_35_spread$c
  }
  
  ## Poly T track ##
  
  # find expected number of T's in 9 bases in the non-splicing intronic region (pos 1-20)
  expected_val <- 0
  for (i in 1:20) {
    expected_val <- expected_val + spread$t[i]
  }
  expected_val <- expected_val / 20 * 9
  
  # finds expected number of T's in poly T track (positions 26-34)
  expected_junction_val <- 0
  for (i in 26:34) {
    expected_junction_val<- expected_junction_val + spread$t[i]
  }
  
  # finds actual number of T's in poly T track (positions 26-34)
  actual_val <- 0
  for (i in 26:34) {
    if (substr(sequence, i, i) == "T") {
      actual_val <- actual_val + 1
    }
  }
  
  # dif <- expected_junction_val - expected_val # number of T's needed to score 1 point
  
  if (actual_val == 4) {
    score <- score + 1
  } else if (actual_val > 4) {
    score <- score + 1 + ((1 / 5) * (actual_val - 4))
  }
  
  
  return(score)
}

p.signif <- function(p_val) {
  vtr <- c()
  for (i in 1:length(p_val)) {
    vtr[i] <- 'ns'
    if (p_val[i] <= 0.05 & p_val[i] > 0.01) {
      vtr[i] <- '*'
    } else if (p_val[i] <= 0.01 & p_val[i] > 0.001) {
      vtr[i] <- '**'
    } else if (p_val[i] <= 0.001) {
      vtr[i] <- '***'
    }
  }
  return(vtr)
}      

# All sequences #
ggseqlogo(list(acceptor_sites$acceptor, donor_sites$donor), ncol = 1, method = "prob")

## Manipulating position 38 ##
manipulation_list <- list(acceptor_sites$acceptor, 
                          acceptor_sites$acceptor[substr(acceptor_sites$acceptor, 38, 38) == "A"],
                          acceptor_sites$acceptor[substr(acceptor_sites$acceptor, 38, 38) == "G"],
                          acceptor_sites$acceptor[substr(acceptor_sites$acceptor, 38, 38) == "T"],
                          acceptor_sites$acceptor[substr(acceptor_sites$acceptor, 38, 38) == "C"])
names(manipulation_list) <- c("All splice sites", "Position 38 == A", "Position 38 == G", 
                              "Position 38 == T", "Position 38 == C")

ggseqlogo(manipulation_list, ncol = 1, method = "prob")

setwd("~/Documents/splice_junction_project/Output")

# This is very inefficient but it works
temp <- manipulation_list$`All splice sites`
print(paste0("n=", length(temp)))
x <- findSpread(temp)
write.csv(x, "Standard.csv", row.names = F)

temp <- manipulation_list$`Position 38 == A`
print(paste0("n=", length(temp)))
x <- findSpread(temp)
write.csv(x, "Pos38A.csv", row.names = F)

temp <- manipulation_list$`Position 38 == G`
print(paste0("n=", length(temp)))
x <- findSpread(temp)
write.csv(x, "Pos38G.csv", row.names = F)

temp <- manipulation_list$`Position 38 == T`
print(paste0("n=", length(temp)))
x <- findSpread(temp)
write.csv(x, "Pos38T.csv", row.names = F)

temp <- manipulation_list$`Position 38 == C`
print(paste0("n=", length(temp)))
x <- findSpread(temp)
write.csv(x, "Pos38C.csv", row.names = F)




## Manipulating position 37 ##
manipulation_list <- list(acceptor_sites$acceptor, 
                          acceptor_sites$acceptor[substr(acceptor_sites$acceptor, 37, 37) == "A"],
                          acceptor_sites$acceptor[substr(acceptor_sites$acceptor, 37, 37) == "G"],
                          acceptor_sites$acceptor[substr(acceptor_sites$acceptor, 37, 37) == "T"],
                          acceptor_sites$acceptor[substr(acceptor_sites$acceptor, 37, 37) == "C"])
names(manipulation_list) <- c("All splice sites", "Position 37 == A", "Position 37 == G", 
                              "Position 37 == T", "Position 37 == C")

ggseqlogo(manipulation_list, ncol = 1, method = "prob")

temp <- manipulation_list$`All splice sites`
print(paste0)
x <- findSpread(temp)
write.csv(x, "Standard.csv", row.names = F)

temp <- manipulation_list$`Position 37 == A`
print(paste0("n=", length(temp)))
x <- findSpread(temp)
write.csv(x, "Pos37A.csv", row.names = F)

temp <- manipulation_list$`Position 37 == G`
print(paste0("n=", length(temp)))
x <- findSpread(temp)
write.csv(x, "Pos37G.csv", row.names = F)

temp <- manipulation_list$`Position 37 == T`
print(paste0("n=", length(temp)))
x <- findSpread(temp)
write.csv(x, "Pos37T.csv", row.names = F)

temp <- manipulation_list$`Position 37 == C`
print(paste0("n=", length(temp)))
x <- findSpread(temp)
write.csv(x, "Pos37C.csv", row.names = F)



#################### Splice junction quantifier ####################
setwd("~/Documents/splice_junction_project//Output")
donor_spread <- findSpread(donor_sites$donor)





# Test for donor quantifier. This should give the max value of the standards (~4.69)
print(donor_quantifier(donor_spread, "GTAAGT"))


## Acceptor Sequence Analysis ## 
# donor and acceptor point values are non-comparable
acceptor_spread <- findSpread(acceptor_sites$acceptor)
acceptor_sites <- na.omit(acceptor_sites[c("gene", "exon", "acceptor")])
for (i in 1:dim(acceptor_sites)[1]) {
  acceptor_sites$score[i] <- acceptor_quantifier(acceptor_spread, acceptor_sites$acceptor[i])
}

B_sites <- na.omit(acceptor_sites[acceptor_sites$exon == "Exon B", c("gene", "exon", "acceptor", "score")])
C_sites <- na.omit(acceptor_sites[acceptor_sites$exon == "Exon C", c("gene", "exon", "acceptor", "score")])
D_sites <- na.omit(acceptor_sites[acceptor_sites$exon == "Exon D", c("gene", "exon", "acceptor", "score")])
E_sites <- na.omit(acceptor_sites[acceptor_sites$exon == "Exon E", c("gene", "exon", "acceptor", "score")])

test <- rbind(B_sites, C_sites, D_sites, E_sites)
test <- test[test$score != 0,]

## Produces the Tukey Test values for th
test$factor <- as.factor(test$exon)
aov <- aov(test$score ~ test$factor, data = test)
tukey_table <- TukeyHSD(aov)

tukey_vals <- c()
i <- 1
for (x in (length(tukey_table$`test$factor`) - (length(tukey_table$`test$factor`) / 4) + 1):length(tukey_table$`test$factor`)) {
  tukey_vals[i] <- tukey_table$`test$factor`[x]
  i <- i + 1
}



tukey <- data.frame("group1" = c("Exon C", "Exon D", "Exon E", "Exon D", "Exon E", "Exon E"), 
                    "group2" = c("Exon B", "Exon B", "Exon B", "Exon C", "Exon C", "Exon D"), 
                    "p" = tukey_vals, 
                    "y.position" = c(7, 7.5, 8, 8.5, 9, 9.5),
                    "p.signif" = p.signif(tukey_vals))


ggboxplot(test, x = "exon", y = "score", fill = "exon") + 
  stat_compare_means(method = "anova", label.x = 1, label.y = 10) + 
  stat_pvalue_manual(tukey, label = "p.signif") + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = .5),
        legend.position = "none",
        axis.title.x = element_blank()) + 
  ylab("Score") + ggtitle("Acceptor Sites")

all_acceptor <- rbind(B_sites, C_sites, D_sites, E_sites)
acceptor_gene_scores <- all_acceptor %>%
  group_by(gene) %>% 
  filter(n() == 4) %>%
  summarize(total_score = sum(score)) %>%
  ungroup()




## Donor Sequence Analysis ##
# donor and acceptor point values are non-comparable
donor_spread <- findSpread(donor_sites$donor)
donor_sites <- na.omit(donor_sites[c("gene", "exon", "donor")])
for (i in 1:dim(donor_sites)[1]) {
  donor_sites$score[i] <- donor_quantifier(donor_spread, donor_sites$donor[i])
}

## Find the last A value for every gene
A_sites <- na.omit(donor_sites[grep("Exon A", donor_sites$exon), c("gene", "exon", "donor", "score")])
## Removes any odd values that don't really fit our motif
A_sites <- A_sites[nchar(A_sites$exon) <= 7, ]
## This finds the value of the highest A exon for every gene and manipulates the A_sites variable accordingly
A_sites <- A_sites %>%
  group_by(gene) %>%
  mutate(exon_num = as.numeric(str_extract(exon, "\\d+")),
         exon_num = if_else(is.na(exon_num), 1, exon_num)) %>%
  slice(which.max(exon_num)) %>%
  select(-exon_num)

## Makes all of the exons 'Exon A'.  This is for the comparison step.  We need to do this because we're lumping all of the A values together
A_sites$exon <- "Exon A"

B_sites <- na.omit(donor_sites[donor_sites$exon == "Exon B", c("gene", "exon", "donor", "score")])
C_sites <- na.omit(donor_sites[donor_sites$exon == "Exon C", c("gene", "exon", "donor", "score")])
D_sites <- na.omit(donor_sites[donor_sites$exon == "Exon D", c("gene", "exon", "donor", "score")])

test <- rbind(A_sites, B_sites, C_sites, D_sites)
test <- test[test$score != 0,]

## Produces the Tukey Test values for th
test$factor <- as.factor(test$exon)
aov <- aov(test$score ~ test$factor, data = test)
tukey_table <- TukeyHSD(aov)

tukey_vals <- c()
i <- 1
for (x in (length(tukey_table$`test$factor`) - (length(tukey_table$`test$factor`) / 4) + 1):length(tukey_table$`test$factor`)) {
  tukey_vals[i] <- tukey_table$`test$factor`[x]
  i <- i + 1
}
   

tukey <- data.frame("group1" = c("Exon B", "Exon C", "Exon D", "Exon C", "Exon D", "Exon D"), 
                    "group2" = c("Exon A", "Exon A", "Exon A", "Exon B", "Exon B", "Exon C"), 
                    "p" = tukey_vals, 
                    "y.position" = c(5, 5.5, 6, 6.5, 7, 7.5),
                    "p.signif" = p.signif(tukey_vals))


ggboxplot(test, x = "exon", y = "score", fill = "exon") + 
  stat_compare_means(method = "anova", label.x = 1, label.y = 8) + 
  stat_pvalue_manual(tukey, label = "p.signif") + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = .5),
        legend.position = "none",
        axis.title.x = element_blank()) + 
  ylab("Score") + ggtitle("Donor Sites")



all_donor <- rbind(A_sites, B_sites, C_sites, D_sites)
donor_gene_scores <- all_donor %>%
  group_by(gene) %>% 
  filter(n() == 4) %>%
  summarize(total_score = sum(score)) %>%
  ungroup()

merged_df <- merge(donor_gene_scores, acceptor_gene_scores, by = "gene")
colnames(merged_df)[2:3] <- c("donor_score", "acceptor_score")

merged_df <- merged_df %>%
  mutate(species = str_split(gene, 'OR')[[1]][1],
         receptor_num = paste0('OR', str_split(gene, 'OR')[[1]][2]))

merged_df["species"] <- NA
merged_df["receptor_num"] <- NA
for (i in 1:dim(merged_df)[1]) {
  sts <- str_split(merged_df[i, 'gene'], "OR")
  print(sts)
  merged_df[i, "species"] <- sts[[1]][1]
  merged_df[i, "receptor_num"] <- paste0('OR' ,sts[[1]][2])
}

ggplot(merged_df, aes(x = species, y = acceptor_score, fill = species)) + 
  geom_boxplot(alpha = 0.7, outlier.alpha = 0) + 
  theme_bw() + 
  ggtitle("Acceptor Score by Species") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  
ggplot(merged_df, aes(x = species, y = donor_score, fill = species)) + 
  geom_boxplot(alpha = 0.7, outlier.alpha = 0) + 
  theme_bw() + 
  ggtitle("Donor Score by Species") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")


write.csv(merged_df[1:3], "~/Documents/Splice_Junction_Project/gene_scores.csv")
