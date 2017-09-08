library(readr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)


# Prepare colours for ggplot
basecaller_names <- c("Nanonet v2.0.0",
                      "Albacore v0.8.4", "Albacore v0.9.1", "Albacore v1.0.4", "Albacore v1.1.2", "Albacore v1.2.6", "Albacore v2.0.2",
                      "Scrappie events v1.0.0", "Scrappie raw v1.0.0", "Scrappie events v1.1.0", "Scrappie raw v1.1.0 (rgrgr_r94)",
                      "Chiron (847ad10)")
basecaller_colours <- c("#73B165",                                   # Nanonet colour
                        brewer.pal(9, "Reds")[2:7],                  # Albacore colours
                        "#9E9AC8", "#C09AC8", "#7F77B1", "#A56BB1",  # Scrappie colours
                        "#639CB1")                                   # Chiron colour
names(basecaller_colours) <- basecaller_names
fill_scale <- scale_fill_manual(name = "Basecaller", values = basecaller_colours)



# Load the tables
nanonet_reads <- read_tsv("results/nanonet_reads.tsv", skip = 1, col_names = c("Name", "Length_Nanonet", "Identity_Nanonet", "Rel_len_Nanonet"))
albacore_v0.8.4_reads <- read_tsv("results/albacore_v0.8.4_reads.tsv", skip = 1, col_names = c("Name", "Length_0.8.4", "Identity_0.8.4", "Rel_len_0.8.4"))
albacore_v0.9.1_reads <- read_tsv("results/albacore_v0.9.1_reads.tsv", skip = 1, col_names = c("Name", "Length_0.9.1", "Identity_0.9.1", "Rel_len_0.9.1"))
albacore_v1.0.4_reads <- read_tsv("results/albacore_v1.0.4_reads.tsv", skip = 1, col_names = c("Name", "Length_1.0.4", "Identity_1.0.4", "Rel_len_1.0.4"))
albacore_v1.1.2_reads <- read_tsv("results/albacore_v1.1.2_reads.tsv", skip = 1, col_names = c("Name", "Length_1.1.2", "Identity_1.1.2", "Rel_len_1.1.2"))
albacore_v1.2.6_reads <- read_tsv("results/albacore_v1.2.6_reads.tsv", skip = 1, col_names = c("Name", "Length_1.2.6", "Identity_1.2.6", "Rel_len_1.2.6"))
albacore_v2.0.2_reads <- read_tsv("results/albacore_v2.0.2_reads.tsv", skip = 1, col_names = c("Name", "Length_2.0.2", "Identity_2.0.2", "Rel_len_2.0.2"))
scrappie_events_v1.0.0_reads <- data.frame(Name = numeric(), Length_Scrappie_events_1.0.0 = numeric(), Identity_Scrappie_events_1.0.0 = numeric(), Rel_len_Scrappie_events_1.0.0 = numeric())  # temporary empty table
scrappie_raw_v1.0.0_reads <- data.frame(Name = numeric(), Length_Scrappie_raw_1.0.0 = numeric(), Identity_Scrappie_raw_1.0.0 = numeric(), Rel_len_Scrappie_raw_1.0.0 = numeric())  # temporary empty table
scrappie_events_v1.1.0_reads <- read_tsv("results/scrappie_v1.1.0_events_reads.tsv", skip = 1, col_names = c("Name", "Length_Scrappie_events_1.1.0", "Identity_Scrappie_events_1.1.0", "Rel_len_Scrappie_events_1.1.0"))
scrappie_raw_rgrgr_r94_v1.1.0_reads <- read_tsv("results/scrappie_v1.1.0_raw_rgrgr_r94_reads.tsv", skip = 1, col_names = c("Name", "Length_Scrappie_raw_rgrgr_r94_1.1.0", "Identity_Scrappie_raw_rgrgr_r94_1.1.0", "Rel_len_Scrappie_raw_rgrgr_r94_1.1.0"))
chiron_reads <- data.frame(Name = numeric(), Length_Chiron = numeric(), Identity_Chiron = numeric(), Rel_len_Chiron = numeric())  # temporary empty table

nanonet_assembly <- read_tsv("results/nanonet_assembly.tsv", skip = 1, col_names = c("Name", "Length_Nanonet", "Identity_Nanonet", "Rel_len_Nanonet"))
albacore_v0.8.4_assembly <- read_tsv("results/albacore_v0.8.4_assembly.tsv", skip = 1, col_names = c("Name", "Length_0.8.4", "Identity_0.8.4", "Rel_len_0.8.4"))
albacore_v0.9.1_assembly <- read_tsv("results/albacore_v0.9.1_assembly.tsv", skip = 1, col_names = c("Name", "Length_0.9.1", "Identity_0.9.1", "Rel_len_0.9.1"))
albacore_v1.0.4_assembly <- read_tsv("results/albacore_v1.0.4_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.0.4", "Identity_1.0.4", "Rel_len_1.0.4"))
albacore_v1.1.2_assembly <- read_tsv("results/albacore_v1.1.2_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.1.2", "Identity_1.1.2", "Rel_len_1.1.2"))
albacore_v1.2.6_assembly <- read_tsv("results/albacore_v1.2.6_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.2.6", "Identity_1.2.6", "Rel_len_1.2.6"))
albacore_v2.0.2_assembly <- read_tsv("results/albacore_v2.0.2_assembly.tsv", skip = 1, col_names = c("Name", "Length_2.0.2", "Identity_2.0.2", "Rel_len_2.0.2"))
scrappie_events_v1.0.0_assembly <- data.frame(Name = numeric(), Length_Scrappie_events_1.0.0 = numeric(), Identity_Scrappie_events_1.0.0 = numeric(), Rel_len_Scrappie_events_1.0.0 = numeric())  # temporary empty table
scrappie_raw_v1.0.0_assembly <- data.frame(Name = numeric(), Length_Scrappie_raw_1.0.0 = numeric(), Identity_Scrappie_raw_1.0.0 = numeric(), Rel_len_Scrappie_raw_1.0.0 = numeric())  # temporary empty table
scrappie_events_v1.1.0_assembly <- data.frame(Name = numeric(), Length_Scrappie_events_1.1.0 = numeric(), Identity_Scrappie_events_1.1.0 = numeric(), Rel_len_Scrappie_events_1.1.0 = numeric())  # temporary empty table
scrappie_raw_rgrgr_r94_v1.1.0_assembly <- data.frame(Name = numeric(), Length_Scrappie_raw_rgrgr_r94_1.1.0 = numeric(), Identity_Scrappie_raw_rgrgr_r94_1.1.0 = numeric(), Rel_len_Scrappie_raw_rgrgr_r94_1.1.0 = numeric())  # temporary empty table
chiron_assembly <- data.frame(Name = numeric(), Length_Chiron = numeric(), Identity_Chiron = numeric(), Rel_len_Chiron = numeric())  # temporary empty table

nanonet_nanopolish <- data.frame(Name = numeric(), Length_Nanonet = numeric(), Identity_Nanonet = numeric(), Rel_len_Nanonet = numeric())  # temporary empty table
albacore_v0.8.4_nanopolish <- data.frame(Name = numeric(), Length_0.8.4 = numeric(), Identity_0.8.4 = numeric(), Rel_len_0.8.4 = numeric())  # temporary empty table
albacore_v0.9.1_nanopolish <- data.frame(Name = numeric(), Length_0.9.1 = numeric(), Identity_0.9.1 = numeric(), Rel_len_0.9.1 = numeric())  # temporary empty table
albacore_v1.0.4_nanopolish <- read_tsv("results/albacore_v1.0.4_nanopolished_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.0.4", "Identity_1.0.4", "Rel_len_1.0.4"))
albacore_v1.1.2_nanopolish <- read_tsv("results/albacore_v1.1.2_nanopolished_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.1.2", "Identity_1.1.2", "Rel_len_1.1.2"))
albacore_v1.2.6_nanopolish <- read_tsv("results/albacore_v1.2.6_nanopolished_assembly.tsv", skip = 1, col_names = c("Name", "Length_1.2.6", "Identity_1.2.6", "Rel_len_1.2.6"))
albacore_v2.0.2_nanopolish <- data.frame(Name = numeric(), Length_2.0.2 = numeric(), Identity_2.0.2 = numeric(), Rel_len_2.0.2 = numeric())  # temporary empty table, until Nanopolish supports Albacore v2
scrappie_events_v1.0.0_nanopolish <- data.frame(Name = numeric(), Length_Scrappie_events_1.0.0 = numeric(), Identity_Scrappie_events_1.0.0 = numeric(), Rel_len_Scrappie_events_1.0.0 = numeric())  # temporary empty table
scrappie_raw_v1.0.0_nanopolish <- data.frame(Name = numeric(), Length_Scrappie_raw_1.0.0 = numeric(), Identity_Scrappie_raw_1.0.0 = numeric(), Rel_len_Scrappie_raw_1.0.0 = numeric())  # temporary empty table
scrappie_events_v1.1.0_nanopolish <- data.frame(Name = numeric(), Length_Scrappie_events_1.1.0 = numeric(), Identity_Scrappie_events_1.1.0 = numeric(), Rel_len_Scrappie_events_1.1.0 = numeric())  # temporary empty table
scrappie_raw_rgrgr_r94_v1.1.0_nanopolish <- data.frame(Name = numeric(), Length_Scrappie_raw_rgrgr_r94_1.1.0 = numeric(), Identity_Scrappie_raw_rgrgr_r94_1.1.0 = numeric(), Rel_len_Scrappie_raw_rgrgr_r94_1.1.0 = numeric())  # temporary empty table
chiron_nanopolish <- data.frame(Name = numeric(), Length_Chiron = numeric(), Identity_Chiron = numeric(), Rel_len_Chiron = numeric())  # temporary empty table



# Merge the tables together
all_reads <- data.frame(Name = character())
all_reads <- merge(all_reads, nanonet_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, albacore_v0.8.4_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, albacore_v0.9.1_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, albacore_v1.0.4_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, albacore_v1.1.2_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, albacore_v1.2.6_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, albacore_v2.0.2_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, scrappie_events_v1.0.0_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, scrappie_raw_v1.0.0_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, scrappie_events_v1.1.0_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, scrappie_raw_rgrgr_r94_v1.1.0_reads, by=1, all=TRUE)
all_reads <- merge(all_reads, chiron_reads, by=1, all=TRUE)

all_assemblies <- data.frame(Name = numeric())
all_assemblies <- merge(all_assemblies, nanonet_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, albacore_v0.8.4_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, albacore_v0.9.1_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, albacore_v1.0.4_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, albacore_v1.1.2_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, albacore_v1.2.6_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, albacore_v2.0.2_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, scrappie_events_v1.0.0_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, scrappie_raw_v1.0.0_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, scrappie_events_v1.1.0_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, scrappie_raw_rgrgr_r94_v1.1.0_assembly, by=1, all=TRUE)
all_assemblies <- merge(all_assemblies, chiron_assembly, by=1, all=TRUE)

all_nanopolish <- data.frame(Name = numeric())
all_nanopolish <- merge(all_nanopolish, nanonet_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, albacore_v0.8.4_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, albacore_v0.9.1_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, albacore_v1.0.4_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, albacore_v1.1.2_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, albacore_v1.2.6_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, albacore_v2.0.2_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, scrappie_events_v1.0.0_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, scrappie_raw_v1.0.0_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, scrappie_events_v1.1.0_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, scrappie_raw_rgrgr_r94_v1.1.0_nanopolish, by=1, all=TRUE)
all_nanopolish <- merge(all_nanopolish, chiron_nanopolish, by=1, all=TRUE)



# Each read's length is the median value of the different basecallers' lengths.
read_lengths <- all_reads[grepl("Length_", names(all_reads))]
all_reads["Length"] <- round(apply(read_lengths, 1, median, na.rm = TRUE))

assembly_lengths <- all_assemblies[grepl("Length_", names(all_assemblies))]
all_assemblies["Length"] <- round(apply(assembly_lengths, 1, median, na.rm = TRUE))

nanopolish_lengths <- all_nanopolish[grepl("Length_", names(all_nanopolish))]
all_nanopolish["Length"] <- round(apply(nanopolish_lengths, 1, median, na.rm = TRUE))



# If a read lacks an identity, then it's an unaligned read.
# Replace NA with 0.0 in identity columns (unless the whole column is NA, which implies that it's a basecaller for which I don't have data).
for(col_name in names(all_reads)){
  if (startsWith(col_name, "Identity_")) {
    identities <- all_reads[,col_name]
    if (!(all(is.na(identities)))) {
      identities[is.na(identities)] <- 0
      all_reads[,col_name] <- identities
    }
  }
}



# Reformat data frames for ggplot
basecaller_identities <- c("Identity_Nanonet",
                           "Identity_0.8.4", "Identity_0.9.1", "Identity_1.0.4", "Identity_1.1.2", "Identity_1.2.6", "Identity_2.0.2",
                           "Identity_Scrappie_events_1.0.0", "Identity_Scrappie_raw_1.0.0", "Identity_Scrappie_events_1.1.0", "Identity_Scrappie_raw_rgrgr_r94_1.1.0",
                           "Identity_Chiron")

basecaller_rel_lengths <- c("Rel_len_Nanonet",
                            "Rel_len_0.8.4", "Rel_len_0.9.1", "Rel_len_1.0.4", "Rel_len_1.1.2", "Rel_len_1.2.6", "Rel_len_2.0.2",
                            "Rel_len_Scrappie_events_1.0.0", "Rel_len_Scrappie_raw_1.0.0", "Rel_len_Scrappie_events_1.1.0", "Rel_len_Scrappie_raw_rgrgr_r94_1.1.0",
                            "Rel_len_Chiron")

read_identities <- all_reads[,c("Name", "Length", basecaller_identities)]
colnames(read_identities) <- c("Name", "Length", basecaller_names)
read_identities <- melt(read_identities, id=c("Name", "Length"))
colnames(read_identities) <- c("Read_name", "Length", "Basecaller", "Identity")

read_rel_lengths <- all_reads[,c("Name", "Length", basecaller_rel_lengths)]
colnames(read_rel_lengths) <- c("Name", "Length", basecaller_names)
read_rel_lengths <- melt(read_rel_lengths, id=c("Name", "Length"))
colnames(read_rel_lengths) <- c("Read_name", "Length", "Basecaller", "Relative_length")

assembly_identities <- all_assemblies[,c("Name", "Length", basecaller_identities)]
colnames(assembly_identities) <- c("Name", "Length", basecaller_names)
assembly_identities <- melt(assembly_identities, id=c("Name", "Length"))
colnames(assembly_identities) <- c("Read_name", "Length", "Basecaller", "Identity")

assembly_rel_lengths <- all_assemblies[,c("Name", "Length", basecaller_rel_lengths)]
colnames(assembly_rel_lengths) <- c("Name", "Length", basecaller_names)
assembly_rel_lengths <- melt(assembly_rel_lengths, id=c("Name", "Length"))
colnames(assembly_rel_lengths) <- c("Read_name", "Length", "Basecaller", "Relative_length")

nanopolish_identities <- all_nanopolish[,c("Name", "Length", basecaller_identities)]
colnames(nanopolish_identities) <- c("Name", "Length", basecaller_names)
nanopolish_identities <- melt(nanopolish_identities, id=c("Name", "Length"))
colnames(nanopolish_identities) <- c("Read_name", "Length", "Basecaller", "Identity")



# Make aligned proportion data frame
total_length <- sum(all_reads$Length)

nanonet_aligned <- 100.0 * sum(all_reads[all_reads$Identity_Nanonet != 0.0,]$Length) / total_length
nanonet_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_Nanonet == 0.0,]$Length) / total_length

albacore_v0.8.4_aligned <- 100.0 * sum(all_reads[all_reads$Identity_0.8.4 != 0.0,]$Length) / total_length
albacore_v0.8.4_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_0.8.4 == 0.0,]$Length) / total_length

albacore_v0.9.1_aligned <- 100.0 * sum(all_reads[all_reads$Identity_0.9.1 != 0.0,]$Length) / total_length
albacore_v0.9.1_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_0.9.1 == 0.0,]$Length) / total_length

albacore_v1.0.4_aligned <- 100.0 * sum(all_reads[all_reads$Identity_1.0.4 != 0.0,]$Length) / total_length
albacore_v1.0.4_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_1.0.4 == 0.0,]$Length) / total_length

albacore_v1.1.2_aligned <- 100.0 * sum(all_reads[all_reads$Identity_1.1.2 != 0.0,]$Length) / total_length
albacore_v1.1.2_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_1.1.2 == 0.0,]$Length) / total_length

albacore_v1.2.6_aligned <- 100.0 * sum(all_reads[all_reads$Identity_1.2.6 != 0.0,]$Length) / total_length
albacore_v1.2.6_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_1.2.6 == 0.0,]$Length) / total_length

albacore_v2.0.2_aligned <- 100.0 * sum(all_reads[all_reads$Identity_2.0.2 != 0.0,]$Length) / total_length
albacore_v2.0.2_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_2.0.2 == 0.0,]$Length) / total_length

scrappie_events_v1.0.0_aligned <- 100.0 * sum(all_reads[all_reads$Identity_Scrappie_events_1.0.0 != 0.0,]$Length) / total_length
scrappie_events_v1.0.0_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_Scrappie_events_1.0.0 == 0.0,]$Length) / total_length

scrappie_raw_v1.0.0_aligned <- 100.0 * sum(all_reads[all_reads$Identity_Scrappie_raw_1.0.0 != 0.0,]$Length) / total_length
scrappie_raw_v1.0.0_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_Scrappie_raw_1.0.0 == 0.0,]$Length) / total_length

scrappie_events_v1.1.0_aligned <- 100.0 * sum(all_reads[all_reads$Identity_Scrappie_events_1.1.0 != 0.0,]$Length) / total_length
scrappie_events_v1.1.0_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_Scrappie_events_1.1.0 == 0.0,]$Length) / total_length

scrappie_raw_rgrgr_r94_v1.1.0_aligned <- 100.0 * sum(all_reads[all_reads$Identity_Scrappie_raw_rgrgr_r94_1.1.0 != 0.0,]$Length) / total_length
scrappie_raw_rgrgr_r94_v1.1.0_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_Scrappie_raw_rgrgr_r94_1.1.0 == 0.0,]$Length) / total_length

chiron_aligned <- 100.0 * sum(all_reads[all_reads$Identity_Chiron != 0.0,]$Length) / total_length
chiron_unaligned <- 100.0 * sum(all_reads[all_reads$Identity_Chiron == 0.0,]$Length) / total_length

aligned_proportion <- data.frame(Basecaller = factor(), Aligned = numeric(), Unaligned = numeric())
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Nanonet v2.0.0"), Aligned = nanonet_aligned, Unaligned = nanonet_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Albacore v0.8.4"), Aligned = albacore_v0.8.4_aligned, Unaligned = albacore_v0.8.4_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Albacore v0.9.1"), Aligned = albacore_v0.9.1_aligned, Unaligned = albacore_v0.9.1_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Albacore v1.0.4"), Aligned = albacore_v1.0.4_aligned, Unaligned = albacore_v1.0.4_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Albacore v1.1.2"), Aligned = albacore_v1.1.2_aligned, Unaligned = albacore_v1.1.2_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Albacore v1.2.6"), Aligned = albacore_v1.2.6_aligned, Unaligned = albacore_v1.2.6_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Albacore v2.0.2"), Aligned = albacore_v2.0.2_aligned, Unaligned = albacore_v2.0.2_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Scrappie events v1.0.0"), Aligned = scrappie_events_v1.0.0_aligned, Unaligned = scrappie_events_v1.0.0_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Scrappie raw v1.0.0"), Aligned = scrappie_raw_v1.0.0_aligned, Unaligned = scrappie_raw_v1.0.0_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Scrappie events v1.1.0"), Aligned = scrappie_events_v1.1.0_aligned, Unaligned = scrappie_events_v1.1.0_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Scrappie raw v1.1.0"), Aligned = scrappie_raw_rgrgr_r94_v1.1.0_aligned, Unaligned = scrappie_raw_rgrgr_r94_v1.1.0_unaligned))
aligned_proportion <- rbind(aligned_proportion, data.frame(Basecaller = factor("Chiron (847ad10)"), Aligned = chiron_aligned, Unaligned = chiron_unaligned))



# Make a read vs assembly identity data frame
nanonet_read_id <- matrixStats::weightedMedian(all_reads$Identity_Nanonet, all_reads$Length, na.rm = TRUE)
nanonet_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_Nanonet, all_assemblies$Length, na.rm = TRUE)

albacore_v0.8.4_read_id <- matrixStats::weightedMedian(all_reads$Identity_0.8.4, all_reads$Length, na.rm = TRUE)
albacore_v0.8.4_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_0.8.4, all_assemblies$Length, na.rm = TRUE)

albacore_v0.9.1_read_id <- matrixStats::weightedMedian(all_reads$Identity_0.9.1, all_reads$Length, na.rm = TRUE)
albacore_v0.9.1_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_0.9.1, all_assemblies$Length, na.rm = TRUE)

albacore_v1.0.4_read_id <- matrixStats::weightedMedian(all_reads$Identity_1.0.4, all_reads$Length, na.rm = TRUE)
albacore_v1.0.4_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_1.0.4, all_assemblies$Length, na.rm = TRUE)

albacore_v1.1.2_read_id <- matrixStats::weightedMedian(all_reads$Identity_1.1.2, all_reads$Length, na.rm = TRUE)
albacore_v1.1.2_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_1.1.2, all_assemblies$Length, na.rm = TRUE)

albacore_v1.2.6_read_id <- matrixStats::weightedMedian(all_reads$Identity_1.2.6, all_reads$Length, na.rm = TRUE)
albacore_v1.2.6_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_1.2.6, all_assemblies$Length, na.rm = TRUE)

albacore_v2.0.2_read_id <- matrixStats::weightedMedian(all_reads$Identity_2.0.2, all_reads$Length, na.rm = TRUE)
albacore_v2.0.2_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_2.0.2, all_assemblies$Length, na.rm = TRUE)

scrappie_events_v1.0.0_read_id <- matrixStats::weightedMedian(all_reads$Identity_Scrappie_events_1.0.0, all_reads$Length, na.rm = TRUE)
scrappie_events_v1.0.0_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_Scrappie_events_1.0.0, all_assemblies$Length, na.rm = TRUE)

scrappie_raw_v1.0.0_read_id <- matrixStats::weightedMedian(all_reads$Identity_Scrappie_raw_1.0.0, all_reads$Length, na.rm = TRUE)
scrappie_raw_v1.0.0_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_Scrappie_raw_1.0.0, all_assemblies$Length, na.rm = TRUE)

scrappie_events_v1.1.0_read_id <- matrixStats::weightedMedian(all_reads$Identity_Scrappie_events_1.1.0, all_reads$Length, na.rm = TRUE)
scrappie_events_v1.1.0_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_Scrappie_events_1.1.0, all_assemblies$Length, na.rm = TRUE)

scrappie_raw_rgrgr_r94_v1.1.0_read_id <- matrixStats::weightedMedian(all_reads$Identity_Scrappie_raw_rgrgr_r94_1.1.0, all_reads$Length, na.rm = TRUE)
scrappie_raw_rgrgr_r94_v1.1.0_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_Scrappie_raw_rgrgr_r94_1.1.0, all_assemblies$Length, na.rm = TRUE)

chiron_read_id <- matrixStats::weightedMedian(all_reads$Identity_Chiron, all_reads$Length, na.rm = TRUE)
chiron_assembly_id <- matrixStats::weightedMedian(all_assemblies$Identity_Chiron, all_assemblies$Length, na.rm = TRUE)

read_vs_assembly_identity <- data.frame(Basecaller = factor(), Read_identity = numeric(), Assembly_identity = numeric())
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Nanonet v2.0.0"), Read_identity = nanonet_read_id, Assembly_identity = nanonet_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Albacore v0.8.4"), Read_identity = albacore_v0.8.4_read_id, Assembly_identity = albacore_v0.8.4_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Albacore v0.9.1"), Read_identity = albacore_v0.9.1_read_id, Assembly_identity = albacore_v0.9.1_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Albacore v1.0.4"), Read_identity = albacore_v1.0.4_read_id, Assembly_identity = albacore_v1.0.4_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Albacore v1.1.2"), Read_identity = albacore_v1.1.2_read_id, Assembly_identity = albacore_v1.1.2_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Albacore v1.2.6"), Read_identity = albacore_v1.2.6_read_id, Assembly_identity = albacore_v1.2.6_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Albacore v2.0.2"), Read_identity = albacore_v2.0.2_read_id, Assembly_identity = albacore_v2.0.2_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Scrappie events v1.0.0"), Read_identity = scrappie_events_v1.0.0_read_id, Assembly_identity = scrappie_events_v1.0.0_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Scrappie raw v1.0.0"), Read_identity = scrappie_raw_v1.0.0_read_id, Assembly_identity = scrappie_raw_v1.0.0_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Scrappie events v1.1.0"), Read_identity = scrappie_events_v1.1.0_read_id, Assembly_identity = scrappie_events_v1.1.0_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Scrappie raw v1.1.0"), Read_identity = scrappie_raw_rgrgr_r94_v1.1.0_read_id, Assembly_identity = scrappie_raw_rgrgr_r94_v1.1.0_assembly_id))
read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = factor("Chiron (847ad10)"), Read_identity = chiron_read_id, Assembly_identity = chiron_assembly_id))



















# Finally, time for plots!
my_theme <- theme_bw() + theme(panel.grid.major.x = element_blank())



# Bar plot of unaligned fraction
ggplot(aligned_proportion, aes(x = Basecaller, y = Unaligned, fill = Basecaller)) +
  geom_bar(stat="identity", colour="black") +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 2), minor_breaks = seq(0, 100, 0.5), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x,fixed=TRUE)) +
  coord_cartesian(ylim=c(0, 16.0)) +
  labs(title = "Unaligned reads", x = "", y = "")



# Violin plots
ggplot(read_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) + 
  geom_violin(draw_quantiles = c(0.5)) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x,fixed=TRUE)) +
  coord_cartesian(ylim=c(65, 100)) +
  labs(title = "Read identities", x = "", y = "")

ggplot(read_rel_lengths, aes(x = Basecaller, y = Relative_length, weight = Length, fill = Basecaller)) + 
  geom_hline(yintercept = 100) + 
  geom_violin(draw_quantiles = c(0.5)) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 200, 2), minor_breaks = seq(0, 200, 1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x,fixed=TRUE)) +
  coord_cartesian(ylim=c(86, 114)) +
  labs(title = "Relative read lengths", x = "", y = "")

ggplot(assembly_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) + 
  geom_violin(draw_quantiles = c(0.5)) +
  fill_scale + my_theme + guides(fill=FALSE) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.5), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x,fixed=TRUE)) +
  coord_cartesian(ylim=c(98.5, 100)) +
  labs(title = "Assembly identities (pre-Nanopolish)", x = "", y = "")

ggplot(assembly_rel_lengths, aes(x = Basecaller, y = Relative_length, weight = Length, fill = Basecaller)) + 
  geom_hline(yintercept = 100) + 
  geom_violin(draw_quantiles = c(0.5)) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 200, 0.5), minor_breaks = seq(0, 200, 0.1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x,fixed=TRUE)) +
  coord_cartesian(ylim=c(98.5, 101.5)) +
  labs(title = "Relative assembly lengths", x = "", y = "")

ggplot(nanopolish_identities, aes(x = factor(Basecaller), y = Identity, weight = Length, fill = Basecaller)) + 
  geom_violin(draw_quantiles = c(0.5)) +
  fill_scale + my_theme + guides(fill=FALSE) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.5), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x,fixed=TRUE)) +
  coord_cartesian(ylim=c(98.5, 100)) +
  labs(title = "Assembly identities (post-Nanopolish)", x = "", y = "")



# Scatter plot
ggplot(read_vs_assembly_identity, aes(x = Read_identity, y = Assembly_identity, fill = Basecaller)) + 
  geom_point(shape = 21, size = 5, stroke = 0.5, alpha = 0.8) +
  fill_scale + theme_bw() +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.5), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
  coord_cartesian(xlim=c(80, 100), ylim=c(98.5, 100)) +
  labs(title = "Read and assembly identities", x = "", y = "")



# This code produces a single plot made of two violin plots:
# * one for the majority of the read identity at the top of the range
# * one for the unaligned reads at the bottom of the range
p1 <- ggplot(read_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) + 
  geom_violin(draw_quantiles = c(0.5)) +
  fill_scale + my_theme + guides(fill=FALSE) +
  theme(axis.ticks.x = element_blank()) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=NULL) +
  coord_cartesian(ylim=c(65, 100)) +
  labs(title = "Read identities", x = "", y = "")
p2 <- ggplot(read_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) + 
  geom_violin(draw_quantiles = c(0.5)) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x,fixed=TRUE)) +
  coord_cartesian(ylim=c(0, 5)) +
  labs(x = "", y = "")
gA <- ggplot_gtable(ggplot_build(p1))
gB <- ggplot_gtable(ggplot_build(p2))
maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
grid.arrange(gA, gB, ncol=1, heights=c(3.5, 1))



# # Joyplots
#
# library(ggjoy)
#
# ggplot(read_identities, aes(x = Identity, y = Basecaller, weight = Length, fill = Basecaller)) +
#   geom_joy(scale = 0.9, draw_quantiles = c(0.5)) +
#   fill_scale + theme_bw() + guides(fill=FALSE) + theme(axis.text.y = element_text(vjust = 0)) +
#   scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
#   scale_y_discrete(expand = c(0.05, 0)) +
#   coord_cartesian(xlim=c(70, 100)) +
#   labs(title = "Read identities", x = "", y = "")
#
# ggplot(assembly_identities, aes(x = Identity, y = Basecaller, weight = Length, fill = Basecaller)) +
#   geom_joy(scale = 1.0, draw_quantiles = c(0.5)) +
#   fill_scale + theme_bw() + guides(fill=FALSE) + theme(axis.text.y = element_text(vjust = 0)) +
#   scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.5), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
#   scale_y_discrete(expand = c(0.05, 0)) +
#   coord_cartesian(xlim=c(98, 100)) +
#   labs(title = "Assembly identities (pre-Nanopolish)", x = "", y = "")
#
# ggplot(nanopolish_identities, aes(x = Identity, y = Basecaller, weight = Length, fill = Basecaller)) +
#   geom_joy(scale = 0.9, draw_quantiles = c(0.5)) +
#   fill_scale + theme_bw() + guides(fill=FALSE) + theme(axis.text.y = element_text(vjust = 0)) +
#   scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.5), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
#   scale_y_discrete(expand = c(0.05, 0)) +
#   coord_cartesian(xlim=c(98, 100)) +
#   labs(title = "Assembly identities (post-Nanopolish)", x = "", y = "")
