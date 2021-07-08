# Plot gene cluster profiles across dataset, with relevant metadata

# These plots are a recreation of the profile plots included in Graphia (graphia.app)

# Author: Josh Harling-Lee

# Part of the GraPPLE package (github.com/JDHarlingLee/GraPPLE)

# load required packages
library(dplyr, quietly = T)
library(ggplot2, quietly = T)
library(cowplot, quietly = T)
library(optparse, quietly = T)

# OPTIONS
option_list <- list(
  make_option(c("-g", "--gene_presc_absc"), type = "character", metavar = "FILE",
              help = "Required: binary gene presence absence file, .tab/.tsv file"),
  
  make_option(c("-c", "--clusters"), type = "character", metavar = "FILE",
              help = "Required: list of genes and clusters, .csv file"),
  
  make_option(c("-m", "--metadata"), type = "character", metavar = "FILE",
              help = "Required: metadata for genomes (e.g. lineage, location, cluster), .csv file"),
  
  make_option(c("-o", "--output_dir"), type = "character", metavar = "DIR",
              help = "Directory for printing gene cluster plots too"),
  
  make_option(c("-N", "--N_metadata_groups"), type = "character", default = "10",
              help = "number or list specifying number of isolate metadata groupings - e.g. \"5,5,10\""),
  
  make_option(c("-p", "--plot_colours"), type = "character", default = "#7187a8",
              help = "colour or list of colours for gene cluster profiles"),
  
  make_option(c("-l", "--list_of_clusters"), type = "character",
              help = "Optional: list of clusters (number only) to output - e.g. \"1,2,5,8,10\". Will override "),
  
  make_option(c("-s", "--start"), type = "integer", default = 1,
              help = "starting cluster number"),
  
  make_option(c("-e", "--end"), type = "integer", default = 20,
              help = "end cluster number")
)

args <- parse_args(OptionParser(option_list = option_list))

# for debugging
#print(args)
#stop(0)

# input parameters
gene_p_a_file <- args$gene_presc_absc
clusters_file <- args$clusters
metadata_file <- args$metadata
out_dir <- args$output_dir

if(is.null(gene_p_a_file) || is.null(clusters_file) || is.null(metadata_file) || is.null(out_dir)){
  stop("Required variables not provided. Check gene pres_absc, gene clusters, metadata files and output directory options") 
}

N_meta_vars <- unlist(strsplit(args$N_metadata_groups, split = ","))
plot_colour <- unlist(strsplit(args$plot_colours, split = ","))
clusters_oi <- args$list_of_clusters
if (!is.null(clusters_oi)) {clusters_oi <- unlist(strsplit(clusters_oi, split = ","))}
start_cluster <- args$start
end_cluster <- args$end
if (end_cluster < start_cluster) {stop("End cluster must be a larger number than start cluster!")}

# read in data
gene_p_a <- read.csv(gene_p_a_file, header = F, stringsAsFactors = F, sep = "\t")
gene_clu <- read.csv(clusters_file)
isol_meta <- read.csv(metadata_file, sep = ",", stringsAsFactors = F)
out_dir <- args$output_dir

# file checks
if(sum(isol_meta[,1]%in%gene_p_a[1,])<=2){
  stop("genomes don't match... Check names in isolate metadata match those in gene_p_a file, and are in the first column of both files")
}

if(dim(gene_clu)[2]!=2){
  stop("wrong number of columns in gene cluster file")
}

if(dim(isol_meta)[2]>4){
  print("Warning: more than 3 metadata columns detected, this may result in cramped plots!")
}

if(dir.exists(out_dir)){
  print("Warning: output directory already exists")
} else {
  dir.create(out_dir)
}

# determine which clusters to investigate from options (default first 20 clusters)
if(is.null(clusters_oi)){
  clusters_oi <- start_cluster:end_cluster
}

# force N_meta_vars to numerical vector of correct length
N_meta_vars <- rep(as.integer(N_meta_vars), ncol(isol_meta))[1:ncol(isol_meta)]

# for plot colours to correct length
plot_colour <- rep(plot_colour, length(clusters_oi))[1:length(clusters_oi)]

# for debug
#str(as.list(.GlobalEnv))
#stop(0)



#### PROCESS METADATA ####
# Set col names
meta_vars <- colnames(isol_meta)[2:ncol(isol_meta)]
colnames(isol_meta)[1] <- "Isolates"

# Remove blanks/NAs
isol_meta[is.na(isol_meta)|isol_meta==""] <- "na" # replace <NA> and blanks with "na" - this is very rudimentary!

# Simplify metadata variables to only top N variables
for (i in seq_along(meta_vars)){
  j <- i+1 # as first column in genome names
  N <- min(N_meta_vars[i], length(unique(isol_meta[,j]))) # ensures no error if fewer than N levels exist in metadata
  
  isol_meta %>% 
    group_by(isol_meta[,j]) %>% 
    summarise(count=n()) %>% 
    arrange(-count) -> tmp
  
  colnames(tmp) <- c("variable", "count")
  topN_meta <- tmp$variable[1:N]
  isol_meta[,j][!isol_meta[,j]%in%topN_meta] <- "Other" # set all other levels to Other
}

isol_meta <- data.frame(isol_meta[1], apply(isol_meta[2:ncol(isol_meta)], 2, function(x) as.factor(x)))
#isol_meta$ST <- factor(isol_meta$ST, levels(isol_meta$ST)[c(1,2,3,4,7,8,9,10,11,12,5,6)])

#### PROCESS GENE INFORMATION ####
# Process gene_p_a file
colnames(gene_p_a) <- gene_p_a[1,]; gene_p_a <- gene_p_a[-1,]
rownames(gene_p_a) <- gene_p_a[,1]; gene_p_a <- gene_p_a[,-1]
r2 <- rownames(gene_p_a)
gene_p_a <- as.data.frame(apply(gene_p_a, 2, function(x) as.numeric(as.character(x))))
gene_p_a$genes <- r2

# Process gene clusters file
colnames(gene_clu) <- c("Gene", "Cluster")
gene_clu$Gene <- as.character(unlist(gene_clu$Gene))

for (i in seq_along(clusters_oi)){
  cluster <- paste("Cluster", clusters_oi[i])
  if(!cluster%in%gene_clu$Cluster){
    stop(paste(cluster[i], ": Invalid cluster # provided", sep = ""))
  }
}

# calculate cluster means per genome
genes_clusters <- inner_join(gene_p_a, gene_clu, by = c("genes" = "Gene"))

cluster_means <- genes_clusters %>% 
  select(-c(genes)) %>% 
  group_by(Cluster) %>% 
  summarise_all(mean) 

cluster_means <- as.data.frame(t(cluster_means))
colnames(cluster_means) <- as.character(unlist(cluster_means[1,])); cluster_means <- cluster_means[-1,] # move clusters to colnames

r3 <- rownames(cluster_means)
cluster_means <- as.data.frame(apply(cluster_means, 2, function(x) as.numeric(as.character(unlist(x)))))
cluster_means$Isolates <- r3 #add rownames as a col for inner_join

# combine metadata to cluster averages
cluster_means_meta <- inner_join(isol_meta, cluster_means, by = c("Isolates"="Isolates"))

# reorder dataframe by metadata - note this will only work on the first 3 columns
cluster_means_meta$Isolates <- factor(cluster_means_meta$Isolates,
                                      levels = cluster_means_meta$Isolates[order(cluster_means_meta[,2]
                                                                                 , cluster_means_meta[,3]
                                                                                 , cluster_means_meta[,4])])

#### PRODUCE GGPLOTS ####

## METADATA PLOTS
meta_plots <- list()
for (i in seq_along(meta_vars)){
  meta_plots[[i]] <- local({
    i <- i # assign local variable
    
    nb.cols <- length(unique(cluster_means_meta[,i+1]))
    # colSet <- paste("Set", i, sep = "")
    #colours <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
    colours <- grDevices::rainbow(N_meta_vars[i]+1, s = 0.6)
    meta_plot <- ggplot(cluster_means_meta) +
      geom_tile(aes(x=Isolates, y=meta_vars[i], fill=cluster_means_meta[,i+1])) +
      scale_fill_manual(values = colours) +
      scale_y_discrete(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0)) +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(color = "black", size = 10),
            axis.title.y = element_blank(),
            legend.position = "none",
            panel.border = element_rect(colour = "#4d4e4f", fill = NA, size = 0.5),
            plot.margin = unit(c(0, 0.1, 0, 0), "cm"))
  })
}

metadata_plots <- plot_grid(plotlist = meta_plots, ncol = 1, align = "v")


## GENE PROFILE PLOT FUNCTION
# define function to plot a full profile graph for any given cluster
quick_plot <- function(cluster, plot_colour){
  no_genes_in_cluster <- nrow(gene_clu[gene_clu$Cluster==cluster,])
  plt <- cluster_means_meta %>% select(Isolates,
                                       ends_with(cluster),
                                       any_of(meta_vars)
                                       )
  g1 <- ggplot(plt) +
    geom_col(aes(x=Isolates, y=plt[,2]), width = .8, color = plot_colour, fill = plot_colour) +
    ggtitle(paste(cluster), subtitle = paste(no_genes_in_cluster, "genes")) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0), labels = scales::percent) +
#    ylab("Av. proportion of genes in genome") +
    theme_classic() +
    theme(plot.title = element_text(size=12),
          plot.subtitle = element_text(size=10),
          axis.text.x = element_blank(),
          axis.text.y = element_text(vjust = 0.2, size=10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          #          axis.title.y = element_text(color = "black", vjust = 5, size=7),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.5, 0.1, 0.1, 0), "cm"))
  all_plots <- list()
  all_plots[[1]] <- g1
  all_plots <- append(all_plots, meta_plots)

  print(cluster)
  
  q <- plot_grid(plotlist = all_plots, ncol = 1, rel_heights = c(8,1,1,1,1,1,1,1), align = "v", axis = "l")
  return(q)
}

## LOOP FUNCTION FOR ALL SPECIFIED CLUSTERS
g <- list()

for (i in seq_along(clusters_oi)){
  cluster <- paste("Cluster", clusters_oi[i])
  g[[i]] <- quick_plot(cluster, plot_colour[i])
}

gene_cluster_profile_plots <- plot_grid(plotlist = g, ncol = 1)

## METADATA LEGENDS
# create metadata legends - duplicate of code above, but taking only legend
# kept separate for ease of formatting
meta_legends <- list()
for (i in seq_along(meta_vars)){
  meta_legends[[i]] <- local({
    i <- i # assign local variable
    
    nb.cols <- length(unique(cluster_means_meta[,i+1]))
    colours <- grDevices::rainbow(N_meta_vars[i]+1, s = 0.5)
    
    meta_plot <- get_legend(ggplot(cluster_means_meta) +
      geom_tile(aes(x=Isolates, y=meta_vars[i], fill=cluster_means_meta[,i+1])) +
      scale_fill_manual(values = colours, name = meta_vars[i]) +
      scale_y_discrete(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0)) +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(color = "black", size = 7),
            axis.title.y = element_blank(),
     #       legend.position = "none",
            panel.border = element_rect(colour = "#4d4e4f", fill = NA, size = 1))
    )
  })
}

metadata_legends <- plot_grid(plotlist = meta_legends, nrow = 1)


#### PLOT AND SAVE ALL PLOTS TO OUTPUT DIR ####
# attempt to roughly scale plots to data size
height = max(2, 2*length(g))
width = max(5, dim(cluster_means_meta)[1]/50)

# filename
if(length(clusters_oi)<5){
  # specify clusters where possible
  out_file <- paste(out_dir, "/gene_cluster_profiles_", paste(clusters_oi, collapse = "-"), sep = "")
} else {
  # avoid printing massive file name if lots of clusters!
  out_file <- paste(out_dir, "/gene_cluster_profiles_", clusters_oi[1], "-to-", tail(clusters_oi, 1), sep = "")
}

# print plots to output directory
print("Writing to file...")
# cluster plots
ggsave(paste(out_file, ".png", sep = ""), gene_cluster_profile_plots,
       height = height, width = width, units = "in", limitsize = FALSE)
# legends
ggsave(paste(out_file, "_legend.png", sep = ""), metadata_legends,
       height = 4, width = 4, units = "in", limitsize = FALSE)

print("Script Completed")




