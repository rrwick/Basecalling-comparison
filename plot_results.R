library(readr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)


basecaller_names <- c("Nanonet v2.0.0",
                      "Albacore v0.8.4", "Albacore v0.9.1", "Albacore v1.0.4", "Albacore v1.1.2", "Albacore v1.2.6", "Albacore v2.0.2",
                      "Scrappie events v1.0.0", "Scrappie events v1.1.0",
                      "Scrappie raw v1.0.0", "Scrappie raw v1.1.0 raw_r94", "Scrappie raw v1.1.0 rgr_r94", "Scrappie raw v1.1.0 rgrgr_r94",
                      "Chiron 847ad10")
basecaller_colours <- c("#F0E57F",                                   # Nanonet colour
                        brewer.pal(9, "Reds")[2:7],                  # Albacore colours
                        "#788CC8", "#6175B1",                        # Scrappie event colours
                        "#C09AC8", "#A56BB1", "#A56BB1", "#A56BB1",  # Scrappie raw colours
                        "#43B164")                                   # Chiron colour
names(basecaller_colours) <- basecaller_names
fill_scale <- scale_fill_manual(name = "Basecaller", values = basecaller_colours)
my_theme <- theme_bw() + theme(panel.grid.major.x = element_blank())



load_tsv_data <- function(filename, column_names) {
  if(file.exists(filename)) {
    data <- read_tsv(filename, skip = 1, col_names = column_names)
  }
  else {
    data <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(data) <- column_names
  }
  return(data)
}



all_reads <- data.frame(Name = character())
all_assemblies <- data.frame(Name = numeric())
all_nanopolish <- data.frame(Name = numeric())
basecaller_identities <- c()
basecaller_rel_lengths <- c()

for (basecaller in basecaller_names) {
  no_spaces <- gsub(" ", "_", basecaller)
  
  read_data_filename <- paste("results/", tolower(no_spaces), "_reads.tsv", sep="")
  assembly_data_filename <- paste("results/", tolower(no_spaces), "_assembly.tsv", sep="")
  nanopolish_data_filename <- paste("results/", tolower(no_spaces), "_nanopolished_assembly.tsv", sep="")
  
  length_column <- paste("Length_", no_spaces, sep="")
  identity_column <- paste("Identity_", no_spaces, sep="")
  rel_length_column <- paste("Rel_len_", no_spaces, sep="")
  
  basecaller_identities <- c(basecaller_identities, identity_column)
  basecaller_rel_lengths <- c(basecaller_rel_lengths, rel_length_column)
  
  column_names = c("Name", length_column, identity_column, rel_length_column)
  read_data <- load_tsv_data(read_data_filename, column_names)
  assembly_data <- load_tsv_data(assembly_data_filename, column_names)
  nanopolish_data <- load_tsv_data(nanopolish_data_filename, column_names)

  all_reads <- merge(all_reads, read_data, by=1, all=TRUE)
  all_assemblies <- merge(all_assemblies, assembly_data, by=1, all=TRUE)
  all_nanopolish <- merge(all_nanopolish, nanopolish_data, by=1, all=TRUE)
}


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


# Load the total bases called data frame.
total_bases <- read_tsv("results/read_counts_and_yields.tsv", skip = 1, col_names = c("Basecaller", "Read_count", "Total_bases"))


# Prepare the data frame for the read vs assembly identity scatterplot.
read_vs_assembly_identity <- data.frame(Basecaller = factor(), Read_identity = numeric(), Assembly_identity = numeric())
for (basecaller in basecaller_names) {
  no_spaces <- gsub(" ", "_", basecaller)
  identity_column <- paste("Identity_", no_spaces, sep="")
  read_ids <- all_reads[,identity_column]
  assembly_ids <- all_assemblies[,identity_column]
  if (!all(is.na(read_ids)) && !all(is.na(assembly_ids)) > 0) {
    median_read_id <- matrixStats::weightedMedian(read_ids, all_reads$Length, na.rm = TRUE)
    median_assembly_id <- matrixStats::weightedMedian(assembly_ids, all_assemblies$Length, na.rm = TRUE)
    read_vs_assembly_identity <- rbind(read_vs_assembly_identity, data.frame(Basecaller = basecaller, Read_identity = median_read_id, Assembly_identity = median_assembly_id))
  }
}











# Bar plots
ggplot(total_bases, aes(x = Basecaller, y = Total_bases, fill = Basecaller)) + 
  geom_bar(stat="identity", colour="black", width = 0.8) +
  fill_scale + my_theme + guides(fill=FALSE) +
  labs(title = "Total basecalling yield", x = "", y = "") +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x,fixed=TRUE), limits = basecaller_names) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2000000000, 500000000), minor_breaks = seq(0, 2000000000, 100000000), labels = scales::unit_format("M", 0.000001)) +
  coord_cartesian(ylim=c(0, 1500000000))



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
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 200, 4), minor_breaks = seq(0, 200, 1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=function(x) gsub(" ","\n",x,fixed=TRUE)) +
  coord_cartesian(ylim=c(84, 116)) +
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





# Scatter plot
ggplot(read_vs_assembly_identity, aes(x = Read_identity, y = Assembly_identity, fill = Basecaller)) + 
  geom_point(shape = 21, size = 4, stroke = 0.5, alpha = 0.85) +
  fill_scale + theme_bw() +
  # scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 1), minor_breaks = seq(0, 100, 0.5), labels = scales::unit_format("%")) +
  # scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.2), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
  # coord_cartesian(xlim=c(80, 90), ylim=c(98.5, 100)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 2), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.2), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
  coord_cartesian(xlim=c(80, 100), ylim=c(98.5, 100)) +
  labs(title = "Read and assembly identities", x = "Read identity", y = "Assembly identity")


# This one is square and has the diagonal line.
poly <- data.frame(x=c(0, 100, 100), y=c(0, 0, 100))
ggplot(read_vs_assembly_identity, aes(x = Read_identity, y = Assembly_identity, fill = Basecaller)) + 
  geom_polygon(data=poly, aes(x=x,y=y),alpha=0.5,fill="black") +
  geom_point(shape = 21, size = 2, stroke = 0.5, alpha = 0.85) +
  fill_scale + theme_bw() + theme(aspect.ratio=1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  coord_cartesian(xlim=c(80, 100), ylim=c(80, 100)) +
  labs(title = "Read and assembly identities", x = "Read identity", y = "Assembly identity")


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
