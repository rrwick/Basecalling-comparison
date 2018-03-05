library(readr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

basecaller_names <- c()
basecaller_colours <- c()

basecaller_names <- c(basecaller_names, "Nanonet v2.0.0")
basecaller_colours <- c(basecaller_colours, "#F0DF46")

basecaller_names <- c(basecaller_names, "Albacore v0.8.4", "Albacore v0.9.1", "Albacore v1.1.2", "Albacore v1.2.6", "Albacore v2.0.2", "Albacore v2.1.10")
basecaller_colours <- c(basecaller_colours, "#FCBBA1", "#F29A87", "#E87A6C", "#DF5952", "#D53937", "#CB181D")

basecaller_names <- c(basecaller_names, "Guppy v0.3.0", "Guppy v0.5.1")
basecaller_colours <- c(basecaller_colours, "#E39C54", "#CF7923")

basecaller_names <- c(basecaller_names, "Scrappie events v1.0.0", "Scrappie events v1.3.0")
basecaller_colours <- c(basecaller_colours, "#788CC8", "#6175B1")

basecaller_names <- c(basecaller_names, "Scrappie raw v1.0.0", "Scrappie raw v1.3.0 raw_r94", "Scrappie raw v1.3.0 rgr_r94", "Scrappie raw v1.3.0 rgrgr_r94", "Scrappie raw v1.3.0 rnnrf_r94")
basecaller_colours <- c(basecaller_colours, "#C4B2C8", "#BA9AC0", "#B38ABB", "#AC7BB6", "#A56BB1")

basecaller_names <- c(basecaller_names, "DeepNano e8a621e")
basecaller_colours <- c(basecaller_colours, "#6BB275")

basecaller_names <- c(basecaller_names, "Chiron v0.2", "Chiron v0.3")
basecaller_colours <- c(basecaller_colours, "#7BB8B8", "#1DB2B2")


names(basecaller_colours) <- basecaller_names
fill_scale <- scale_fill_manual(name = "Basecaller", values = basecaller_colours)
my_theme <- theme_bw() + theme(panel.grid.major.x = element_blank())

basecaller_labels <- gsub(" ", "\n", basecaller_names, fixed=TRUE)
basecaller_labels <- gsub("basecRAWller", "base-\ncRAWller", basecaller_labels, fixed=TRUE)
basecaller_labels <- gsub("DeepNano", "Deep-\nNano", basecaller_labels, fixed=TRUE)


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


all_reads <- load_tsv_data("results/read_data.tsv", column_names=c("Name", "Fast5_name", "Run_name", "Signal_length", "Start_time"))
all_assemblies <- data.frame(Name = numeric())
all_nanopolish <- data.frame(Name = numeric())
all_nanopolish_meth <- data.frame(Name = numeric())
all_medaka <- data.frame(Name = numeric())
basecaller_identities <- c()
basecaller_rel_lengths <- c()

for (basecaller in basecaller_names) {
  no_spaces <- gsub(" ", "_", basecaller)

  read_data_filename <- paste("results/", tolower(no_spaces), "_reads.tsv", sep="")
  assembly_data_filename <- paste("results/", tolower(no_spaces), "_assembly.tsv", sep="")
  nanopolish_data_filename <- paste("results/", tolower(no_spaces), "_nanopolish.tsv", sep="")
  nanopolish_meth_data_filename <- paste("results/", tolower(no_spaces), "_nanopolish_meth.tsv", sep="")
  medaka_data_filename <- paste("results/", tolower(no_spaces), "_medaka.tsv", sep="")

  length_column <- paste("Length_", no_spaces, sep="")
  identity_column <- paste("Identity_", no_spaces, sep="")
  rel_length_column <- paste("Rel_len_", no_spaces, sep="")

  basecaller_identities <- c(basecaller_identities, identity_column)
  basecaller_rel_lengths <- c(basecaller_rel_lengths, rel_length_column)

  column_names = c("Name", length_column, identity_column, rel_length_column)
  read_data <- load_tsv_data(read_data_filename, column_names)
  assembly_data <- load_tsv_data(assembly_data_filename, column_names)
  nanopolish_data <- load_tsv_data(nanopolish_data_filename, column_names)
  nanopolish_meth_data <- load_tsv_data(nanopolish_meth_data_filename, column_names)
  medaka_data <- load_tsv_data(medaka_data_filename, column_names)

  all_reads <- merge(all_reads, read_data, by=1, all=TRUE)
  all_assemblies <- merge(all_assemblies, assembly_data, by=1, all=TRUE)
  all_nanopolish <- merge(all_nanopolish, nanopolish_data, by=1, all=TRUE)
  all_nanopolish_meth <- merge(all_nanopolish_meth, nanopolish_meth_data, by=1, all=TRUE)
  all_medaka <- merge(all_medaka, medaka_data, by=1, all=TRUE)
}


# Each read's length is the median value of the different basecallers' lengths.
read_lengths <- all_reads[grepl("Length_", names(all_reads))]
all_reads["Length"] <- round(apply(read_lengths, 1, median, na.rm = TRUE))

assembly_lengths <- all_assemblies[grepl("Length_", names(all_assemblies))]
all_assemblies["Length"] <- round(apply(assembly_lengths, 1, median, na.rm = TRUE))

nanopolish_lengths <- all_nanopolish[grepl("Length_", names(all_nanopolish))]
all_nanopolish["Length"] <- round(apply(nanopolish_lengths, 1, median, na.rm = TRUE))

nanopolish_meth_lengths <- all_nanopolish_meth[grepl("Length_", names(all_nanopolish_meth))]
all_nanopolish_meth["Length"] <- round(apply(nanopolish_meth_lengths, 1, median, na.rm = TRUE))

medaka_lengths <- all_medaka[grepl("Length_", names(all_medaka))]
all_medaka["Length"] <- round(apply(medaka_lengths, 1, median, na.rm = TRUE))


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

nanopolish_meth_identities <- all_nanopolish_meth[,c("Name", "Length", basecaller_identities)]
colnames(nanopolish_meth_identities) <- c("Name", "Length", basecaller_names)
nanopolish_meth_identities <- melt(nanopolish_meth_identities, id=c("Name", "Length"))
colnames(nanopolish_meth_identities) <- c("Read_name", "Length", "Basecaller", "Identity")

medaka_identities <- all_medaka[,c("Name", "Length", basecaller_identities)]
colnames(medaka_identities) <- c("Name", "Length", basecaller_names)
medaka_identities <- melt(medaka_identities, id=c("Name", "Length"))
colnames(medaka_identities) <- c("Read_name", "Length", "Basecaller", "Identity")


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
read_vs_assembly_identity$Basecaller_with_newlines <- sapply(read_vs_assembly_identity$Basecaller, function(x) gsub(" v","\nv",x,fixed=TRUE))
















total_yield_plot <- ggplot(total_bases, aes(x = Basecaller, y = Total_bases, fill = Basecaller)) +
  geom_bar(stat="identity", colour="black", width = 0.8) +
  fill_scale + my_theme + guides(fill=FALSE) +
  labs(title = "", x = "", y = "total basecalling yield") +
  scale_x_discrete(labels=basecaller_labels, limits = basecaller_names) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2000000000, 500000000), minor_breaks = seq(0, 2000000000, 100000000), labels = scales::unit_format("Mbp", 0.000001)) +
  coord_cartesian(ylim=c(0, 1500000000))
total_yield_plot
ggsave(total_yield_plot, file='plots/total_yield.pdf', width = 12, height = 3)

rel_read_length_plot <- ggplot(read_rel_lengths, aes(x = Basecaller, y = Relative_length, weight = Length, fill = Basecaller)) +
  geom_hline(yintercept = 100) +
  geom_violin(draw_quantiles = c(0.5), width=1.1, bw=0.25) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 200, 4), minor_breaks = seq(0, 200, 1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=basecaller_labels) +
  coord_cartesian(ylim=c(88, 112)) +
  labs(title = "", x = "", y = "read length / reference length")
rel_read_length_plot
ggsave(rel_read_length_plot, file='plots/rel_read_length.pdf', width = 12, height = 4)

assembly_identity_plot <- ggplot(assembly_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) +
  geom_violin(draw_quantiles = c(0.5), bw=0.06) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.5), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=basecaller_labels) +
  coord_cartesian(ylim=c(98.5, 100)) +
  labs(title = "", x = "", y = "assembly identity")
assembly_identity_plot
ggsave(assembly_identity_plot, file='plots/assembly_identity.pdf', width = 12, height = 5)

rel_assembly_length_plot <- ggplot(assembly_rel_lengths, aes(x = Basecaller, y = Relative_length, weight = Length, fill = Basecaller)) +
  geom_hline(yintercept = 100) +
  geom_violin(draw_quantiles = c(0.5), bw=0.06) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 200, 0.5), minor_breaks = seq(0, 200, 0.1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=basecaller_labels) +
  coord_cartesian(ylim=c(98.75, 101.25)) +
  labs(title = "", x = "", y = "assembly length / reference length")
rel_assembly_length_plot
ggsave(rel_assembly_length_plot, file='plots/rel_assembly_length.pdf', width = 12, height = 4)

nanopolish_identity_plot <- ggplot(nanopolish_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) +
  geom_violin(data = assembly_identities, draw_quantiles = c(0.5), bw=0.06, alpha=0.2, colour=NA) +
  geom_violin(draw_quantiles = c(0.5), bw=0.06) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.5), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=basecaller_labels) +
  coord_cartesian(ylim=c(98.5, 100)) +
  labs(title = "", x = "", y = "assembly identity")
nanopolish_identity_plot
ggsave(nanopolish_identity_plot, file='plots/nanopolish_identity.pdf', width = 12, height = 5)

nanopolish_meth_identity_plot <- ggplot(nanopolish_meth_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) +
  geom_violin(data = assembly_identities, draw_quantiles = c(0.5), bw=0.06, alpha=0.2, colour=NA) +
  geom_violin(draw_quantiles = c(0.5), bw=0.06) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.5), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=basecaller_labels) +
  coord_cartesian(ylim=c(98.5, 100)) +
  labs(title = "", x = "", y = "assembly identity")
nanopolish_meth_identity_plot
ggsave(nanopolish_meth_identity_plot, file='plots/nanopolish_meth_identity.pdf', width = 12, height = 5)

medaka_identity_plot <- ggplot(medaka_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) +
  geom_violin(data = assembly_identities, draw_quantiles = c(0.5), bw=0.06, alpha=0.2, colour=NA) +
  geom_violin(draw_quantiles = c(0.5), bw=0.06) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.5), minor_breaks = seq(0, 100, 0.1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=basecaller_labels) +
  coord_cartesian(ylim=c(98.5, 100)) +
  labs(title = "", x = "", y = "assembly identity")
medaka_identity_plot
ggsave(medaka_identity_plot, file='plots/medaka_identity.pdf', width = 12, height = 5)


# This code produces a single plot made of two violin plots:
# * one for the majority of the read identity at the top of the range
# * one for the unaligned reads at the bottom of the range
p1 <- ggplot(read_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) +
  geom_violin(draw_quantiles = c(0.5), width=1.2, bw=0.6) +
  fill_scale + my_theme + guides(fill=FALSE) +
  theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=NULL) +
  coord_cartesian(ylim=c(65, 100)) +
  labs(title = "", x = "", y = "read identity")
p2 <- ggplot(read_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) +
  geom_violin(draw_quantiles = c(0.5), width=1.1, bw=0.6) +
  fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=basecaller_labels) +
  coord_cartesian(ylim=c(0, 5)) +
  labs(x = "", y = "")
gA <- ggplot_gtable(ggplot_build(p1))
gB <- ggplot_gtable(ggplot_build(p2))
maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
read_identity_plot <- grid.arrange(gA, gB, ncol=1, heights=c(3, 1))
ggsave(read_identity_plot, file='plots/read_identity.pdf', width = 12, height = 5.75)



# Two-part scatter plot
poly <- data.frame(x=c(0, 100, 100), y=c(0, 0, 100))
p1 <- ggplot(read_vs_assembly_identity, aes(x = Read_identity, y = Assembly_identity, fill = Basecaller)) +
  geom_polygon(data=poly, aes(x=x,y=y),alpha=0.3,fill="black") +
  geom_point(shape = 21, size = 2, stroke = 0.5, alpha = 0.85) +
  fill_scale + theme_bw() + theme(aspect.ratio=1) + guides(fill=FALSE) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1)) +
  coord_cartesian(xlim=c(70, 100), ylim=c(70, 100)) +
  labs(title = "", x = "read identity (%)", y = "assembly identity (%)")
p2 <- ggplot(read_vs_assembly_identity, aes(x = Read_identity, y = Assembly_identity, fill = Basecaller)) +
  geom_point(shape = 21, size = 3, stroke = 0.5, alpha = 0.85) +
  fill_scale + theme_bw() + theme(aspect.ratio=1) +
  guides(fill=guide_legend(title="")) + theme(legend.key.size = unit(0.9, 'lines')) +
  scale_x_continuous(expand = c(0.0, 0.0), breaks = seq(0, 100, 2), minor_breaks = seq(0, 100, 0.5)) +
  scale_y_continuous(expand = c(0.0, 0.0), breaks = seq(0, 100, 0.4), minor_breaks = seq(0, 100, 0.1)) +
  coord_cartesian(xlim=c(82, 89), ylim=c(98.7, 100.0)) +
  labs(title = "", x = "read identity (%)", y = "assembly identity (%)")
read_assembly_scatter_plot <- grid.arrange(p1, p2, ncol=2, widths=c(2,3))
ggsave(read_assembly_scatter_plot, file='plots/read_assembly_scatter.pdf', width = 9, height = 3.25)





# # Joyplots
# library(ggjoy)
#
# ggplot(read_identities, aes(x = Identity, y = Basecaller, weight = Length, fill = Basecaller)) +
#   geom_joy(scale = 2.0, draw_quantiles = c(0.5)) +
#   fill_scale + theme_bw() + guides(fill=FALSE) + theme(axis.text.y = element_text(vjust = 0)) +
#   scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
#   scale_y_discrete(expand = c(0.05, 0)) +
#   coord_cartesian(xlim=c(60, 100)) +
#   labs(title = "Read identities", x = "", y = "")
#
# ggplot(assembly_identities, aes(x = Identity, y = Basecaller, weight = Length, fill = Basecaller)) +
#   geom_joy(scale = 1.8, draw_quantiles = c(0.5)) +
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




# Signal length vs read length plots
for (basecaller in basecaller_names) {
  no_spaces <- gsub(" ", "_", basecaller)
  length_column <- paste("Length_", no_spaces, sep="")
  read_vs_signal_plot <- ggplot(all_reads, aes_string(x = "Signal_length", y = length_column, color="Run_name")) +
    geom_point(size = 0.25, alpha = 1, stroke = 0) +
    scale_color_manual(values=c("#5268AB", "#AD5151", "#AD5151", "#5268AB"),
                       breaks=c("klebs_033_sequencing_run", "klebs_033_restart_sequencing_run"),
                       labels=c("Initial run", "Restart")) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0), labels=scales::unit_format("k", 1e-3, sep="")) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim=c(0, 700000), ylim=c(0, 70000)) +
    labs(title = "", x = "signal length (samples)", y = paste(basecaller, "read length (bp)")) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
    theme(legend.title=element_blank())
  print(read_vs_signal_plot)
  png_filename =  paste("plots/read_vs_signal_", tolower(no_spaces), ".png", sep="")
  ggsave(read_vs_signal_plot, file=png_filename, width = 6, height = 4)
}







# Scrappie v1.1.0 vs v1.1.1
# Note: this section only works if both "Scrappie raw v1.1.0 rgrgr_r94" and
# "Scrappie raw v1.1.1 rgrgr_r94" are included in basecaller_names at the
# top of this script.
scrappie_read_identities <- read_identities[read_identities$Basecaller == "Scrappie raw v1.1.0 rgrgr_r94" | read_identities$Basecaller == "Scrappie raw v1.1.1 rgrgr_r94",]
scrappie_assembly_identities <- assembly_identities[assembly_identities$Basecaller == "Scrappie raw v1.1.0 rgrgr_r94" | assembly_identities$Basecaller == "Scrappie raw v1.1.1 rgrgr_r94",]

scrappie_names <- c("Scrappie raw v1.1.0 rgrgr_r94", "Scrappie raw v1.1.1 rgrgr_r94")
scrappie_labels <- gsub(" ", "\n", scrappie_names, fixed=TRUE)
scrappie_colours <- c("#C4B2C8", "#A56BB1")
names(scrappie_colours) <- scrappie_names
scrappie_fill_scale <- scale_fill_manual(name = "Basecaller", values = scrappie_colours)

p1 <- ggplot(scrappie_read_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) +
  geom_violin(draw_quantiles = c(0.5), bw=0.6) +
  scrappie_fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 5), minor_breaks = seq(0, 100, 1), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=scrappie_labels) +
  coord_cartesian(ylim=c(65.0, 100)) +
  labs(title = "", x = "", y = "read identity")

p2 <- ggplot(scrappie_assembly_identities, aes(x = Basecaller, y = Identity, weight = Length, fill = Basecaller)) +
  geom_violin(draw_quantiles = c(0.5), bw=0.06) +
  scrappie_fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.2), minor_breaks = seq(0, 100, 0.05), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=scrappie_labels) +
  coord_cartesian(ylim=c(99.0, 100)) +
  labs(title = "", x = "", y = "assembly identity")

blank <- rectGrob(gp=gpar(col="white"))
scrappie_comparison_plot <- grid.arrange(p1, blank, p2, ncol=3, widths=c(0.425, 0.15, 0.425))
# ggsave(scrappie_comparison_plot, file='plots/scrappie_comparison.pdf', width = 6, height = 4)











# The following code is for a violin plot comparing different polishing strategies.
basecaller <- "Albacore v2.1.10"
no_spaces <- gsub(" ", "_", basecaller)

polishing_names <- c(paste(basecaller, "(nopolishing)"),
                     paste(basecaller, "+Medaka"),
                     paste(basecaller, "+Nanopolish"),
                     paste(basecaller, "+Medaka +Nanopolish"),
                     paste(basecaller, "+Nanopolish +Medaka"))
polishing_labels <- gsub(" ", "\n", polishing_names, fixed=TRUE)
polishing_labels <- gsub("(nopolishing)", "(no polishing)", polishing_labels, fixed=TRUE)
polishing_colours <- c("#CCCCCC", "#4F84CC", "#DD616C", "#846DC2", "#AB61BA")
names(polishing_colours) <- polishing_names
polishing_fill_scale <- scale_fill_manual(name = "Polishing", values = polishing_colours)

no_polish_data_filename <- paste("results/", tolower(no_spaces), "_assembly.tsv", sep="")
medaka_data_filename <- paste("results/", tolower(no_spaces), "_medaka.tsv", sep="")
nanopolish_data_filename <- paste("results/", tolower(no_spaces), "_nanopolish_meth.tsv", sep="")
medaka_nanopolish_data_filename <- paste("results/", tolower(no_spaces), "_medaka_nanopolish_meth.tsv", sep="")
nanopolish_medaka_data_filename <- paste("results/", tolower(no_spaces), "_nanopolish_meth_medaka.tsv", sep="")

length_column <- paste("Length_", no_spaces, sep="")
rel_length_column <- paste("Rel_len_", no_spaces, sep="")

basecaller_identities <- c()

identity_column <- paste("Identity_", no_spaces, "_no_polish", sep="")
no_polish_data <- load_tsv_data(no_polish_data_filename, c("Name", length_column, identity_column, rel_length_column))
basecaller_identities <- c(basecaller_identities, identity_column)

identity_column <- paste("Identity_", no_spaces, "_medaka", sep="")
medaka_data <- load_tsv_data(medaka_data_filename, c("Name", length_column, identity_column, rel_length_column))
basecaller_identities <- c(basecaller_identities, identity_column)

identity_column <- paste("Identity_", no_spaces, "_nanopolish", sep="")
nanopolish_data <- load_tsv_data(nanopolish_data_filename, c("Name", length_column, identity_column, rel_length_column))
basecaller_identities <- c(basecaller_identities, identity_column)

identity_column <- paste("Identity_", no_spaces, "_medaka_nanopolish", sep="")
medaka_nanopolish_data <- load_tsv_data(medaka_nanopolish_data_filename, c("Name", length_column, identity_column, rel_length_column))
basecaller_identities <- c(basecaller_identities, identity_column)

identity_column <- paste("Identity_", no_spaces, "_nanopolish_medaka", sep="")
nanopolish_medaka_data <- load_tsv_data(nanopolish_medaka_data_filename, c("Name", length_column, identity_column, rel_length_column))
basecaller_identities <- c(basecaller_identities, identity_column)

polish_data <- data.frame(Name = numeric())
polish_data <- merge(polish_data, no_polish_data, by=1, all=TRUE)
polish_data <- merge(polish_data, medaka_data, by=1, all=TRUE)
polish_data <- merge(polish_data, nanopolish_data, by=1, all=TRUE)
polish_data <- merge(polish_data, medaka_nanopolish_data, by=1, all=TRUE)
polish_data <- merge(polish_data, nanopolish_medaka_data, by=1, all=TRUE)

assembly_lengths <- polish_data[grepl("Length_", names(polish_data))]
polish_data["Length"] <- round(apply(assembly_lengths, 1, median, na.rm = TRUE))

polish_identities <- polish_data[,c("Name", "Length", basecaller_identities)]
colnames(polish_identities) <- c("Name", "Length", polishing_names)
polish_identities <- melt(polish_identities, id=c("Name", "Length"))
colnames(polish_identities) <- c("Read_name", "Length", "Polishing", "Identity")

polishing_plot <- ggplot(polish_identities, aes(x = Polishing, y = Identity, weight = Length, fill = Polishing)) +
  geom_violin(draw_quantiles = c(0.5), bw=0.06, width=1.1) +
  polishing_fill_scale + my_theme + guides(fill=FALSE) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.1), minor_breaks = seq(0, 100, 0.05), labels = scales::unit_format("%")) +
  scale_x_discrete(labels=polishing_labels) +
  coord_cartesian(ylim=c(99.2, 100)) +
  labs(title = "", x = "", y = "assembly identity")
polishing_plot
ggsave(polishing_plot, file='plots/polishing_methods.pdf', width = 4.75, height = 4)

