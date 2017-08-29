library(readr)
library(reshape2)
library(ggplot2)


# Load the tables
albacore_v1_1_2_reads <- read_tsv("results/albacore_v1.1.2_reads.tsv", skip = 1, col_names = c("Name", "Length_1.1.2", "Identity_1.1.2", "Rel_len_1.1.2"))
albacore_v1_2_6_reads <- read_tsv("results/albacore_v1.2.6_reads.tsv", skip = 1, col_names = c("Name", "Length_1.2.6", "Identity_1.2.6", "Rel_len_1.2.6"))
albacore_v2_0_0_reads <- read_tsv("results/albacore_v2.0.0_reads.tsv", skip = 1, col_names = c("Name", "Length_2.0.0", "Identity_2.0.0", "Rel_len_2.0.0"))

albacore_v1_1_2_assembly <- read_tsv("results/albacore_v1.1.2_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.1.2", "Identity_1.1.2", "Rel_len_1.1.2"))
albacore_v1_2_6_assembly <- read_tsv("results/albacore_v1.2.6_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.2.6", "Identity_1.2.6", "Rel_len_1.2.6"))
albacore_v2_0_0_assembly <- read_tsv("results/albacore_v2.0.0_assembly.tsv", skip = 1, col_names = c("Name", "Length_2.0.0", "Identity_2.0.0", "Rel_len_2.0.0"))

albacore_v1_1_2_nanopolish <- read_tsv("results/albacore_v1.1.2_nanopolished_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.1.2", "Identity_1.1.2", "Rel_len_1.1.2"))
albacore_v1_2_6_nanopolish <- read_tsv("results/albacore_v1.2.6_nanopolished_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.2.6", "Identity_1.2.6", "Rel_len_1.2.6"))
albacore_v2_0_0_nanopolish <- read_tsv(Name = numeric(0), Length_2.0.0 = numeric(0), Identity_2.0.0 = numeric(0), Rel_len_2.0.0 = numeric(0))  # temporary empty table, until Nanopolish supports Albacore v2



# Merge the tables together
all_reads <- albacore_v1_1_2_reads
all_reads <- merge(all_reads, albacore_v1_2_6_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, albacore_v2_0_0_reads, by=1, all=TRUE)

all_assemblies <- albacore_v1_1_2_assembly
all_assemblies <- merge(all_assemblies, albacore_v1_2_6_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, albacore_v2_0_0_assembly, by=1, all=TRUE)

all_nanopolish <- albacore_v1_1_2_nanopolish
all_nanopolish <- merge(all_nanopolish, albacore_v1_2_6_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, albacore_v2_0_0_nanopolish, by=1, all=TRUE)



# Prepare data frames for ggplot
read_identities <- all_reads[,c("Name", "Length_2.0.0", "Identity_1.1.2", "Identity_1.2.6", "Identity_2.0.0")]
colnames(read_identities) <- c("Name", "Length", "Albacore v1.1.2", "Albacore v1.2.6", "Albacore v2.0.0")
read_identities <- melt(read_identities, id=c("Name", "Length"))
colnames(read_identities) <- c("Read_name", "Length", "Basecaller", "Identity")

read_rel_lengths <- all_reads[,c("Name", "Length_2.0.0", "Rel_len_1.1.2", "Rel_len_1.2.6", "Rel_len_2.0.0")]
colnames(read_rel_lengths) <- c("Name", "Length", "Albacore v1.1.2", "Albacore v1.2.6", "Albacore v2.0.0")
read_rel_lengths <- melt(read_rel_lengths, id=c("Name", "Length"))
colnames(read_rel_lengths) <- c("Read_name", "Length", "Basecaller", "Relative_length")

assembly_identities <- all_assemblies[,c("Name", "Length_2.0.0", "Identity_1.1.2", "Identity_1.2.6", "Identity_2.0.0")]
colnames(assembly_identities) <- c("Name", "Length", "Albacore v1.1.2", "Albacore v1.2.6", "Albacore v2.0.0")
assembly_identities <- melt(assembly_identities, id=c("Name", "Length"))
colnames(assembly_identities) <- c("Read_name", "Length", "Basecaller", "Identity")

nanopolish_identities <- all_nanopolish[,c("Name", "Length_1.2.6", "Identity_1.1.2", "Identity_1.2.6", "Identity_2.0.0")]
colnames(nanopolish_identities) <- c("Name", "Length", "Albacore v1.1.2", "Albacore v1.2.6")
nanopolish_identities <- melt(nanopolish_identities, id=c("Name", "Length"))
colnames(nanopolish_identities) <- c("Read_name", "Length", "Basecaller", "Identity")



# Violin plots
ggplot(read_identities, aes(x = factor(Basecaller), y = Identity, weight = Length, fill = factor(Basecaller))) + 
  geom_violin(draw_quantiles = c(0.5)) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 100, 5),
                     minor_breaks = seq(0, 100, 1),
                     labels = scales::unit_format("%")) +
  coord_cartesian(ylim=c(65, 100))

ggplot(read_rel_lengths, aes(x = factor(Basecaller), y = Relative_length, weight = Length, fill = factor(Basecaller))) + 
  geom_violin(draw_quantiles = c(0.5)) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 200, 2),
                     minor_breaks = seq(0, 200, 1),
                     labels = scales::unit_format("%")) +
  coord_cartesian(ylim=c(92, 108))

ggplot(assembly_identities, aes(x = factor(Basecaller), y = Identity, weight = Length, fill = factor(Basecaller))) + 
  geom_violin(draw_quantiles = c(0.5)) +
  theme_bw() + 
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 100, 0.5),
                     minor_breaks = seq(0, 100, 0.1),
                     labels = scales::unit_format("%")) +
  coord_cartesian(ylim=c(98, 100))

ggplot(nanopolish_identities, aes(x = factor(Basecaller), y = Identity, weight = Length, fill = factor(Basecaller))) + 
  geom_violin(draw_quantiles = c(0.5)) +
  theme_bw() + 
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 100, 0.5),
                     minor_breaks = seq(0, 100, 0.1),
                     labels = scales::unit_format("%")) +
  coord_cartesian(ylim=c(98.5, 100))
