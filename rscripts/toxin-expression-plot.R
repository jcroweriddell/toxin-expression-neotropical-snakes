##########
### Script to identify and rename toxin families based on blast hits & plot ###
##########
### Packages
library(tidyverse)
library(viridis)
library(scales)
library(RColorBrewer)

### Directories
dir_path <- "C:/Users/jmcr/Documents/toxin-expression-neotropical-snakes/"
dir_data <- paste0(dir_path,"data/")
dir_plot <- paste0(dir_path,"plots/")
dir_output <- paste0(dir_path,"outputs/")

### Colour palettes
# palette
plasma<-plasma(10)
show_col(plasma)

# label colours by toxin gene family
name_vector <- c("PLA2", "SVMP", "ThreeFTXs", "CRISP", "KUNITZ", "SVSP", "PLB", "Waprin", "Cystatin", "Other")
alpha_names <- sort(name_vector)
names(plasma) <- alpha_names

###########
### Functions

## Rename toxin families
rename_genes <- function(data){
  # Remove geneFamily column and punctuation from blastHit column
  #data$blastHit <- str_replace_all(data$blastHit, "[[:punct:]]", " ")
  #data$geneFamily <- data$blastHit
  
  # venom genes
  # venom gene family grep searches
  PLA2 <- ".*[pP]hospholipase A2.*|.*PLA2.*|[bB]eta[ -][bB]ungarotoxin.*"
  SVMP <- ".*[Mm]etall.*|.*[Dd]isintegrin.*|.*Poly-His-poly.*|.*(SVMP).*"
  ThreeFTXs <- ".*[Ff][Tt][Xx].*|.*finger.*|.*(?<=fron|neuro|three[ -]finger[ -]|neuro|elapi|cobra|pseudonaja|micruro|colubri|cobro|cyto|cardio|weak[ -])toxin.*|.*Muscarinic.*|.*[Aa]ctiflagelin.*|[Kk]appa[ -]\\d*\\w*-*[bB]ungarotoxin.*|^(3FTX).+e"
  CTL <- ".*C-type [Ll]ecti.*n.*|.*Snaclec.*|.*(CTL).*"
  CRISP <- ".*[Cc]ysteine.*protein.*|.*(CRISP).*"
  KUNITZ <- ".*[kK]unitz.*|.*[Mm]ambaquaretin.*|.*BPTI.*"
  SVSP <- ".*[sS]erine[ -]protease.*|.*[Pp]rotease.*|.*[Ff]ibrinogenase.*|.*[Tt]hrombin.*|.*[Pp]lasminogen.*|.*[pP]roteinase.*|.*(SVSP).*"
  BPP <- ".*[bB]radykinin.*|.*(BPP).*"
  NGF <- ".*[gG]rowth [Ff]actor.*"
  PLB <- ".*[Pp]hospholipase[ -]B.+|.*(PLB).*"
  Waprin <- ".*[Ww]aprin.*"
  SVNP <- ".*[Nn]atriuretic[ -][Pp]eptide.*|.*(SVNP).*"
  Myotoxin <- ".*[mM]yotoxin.*|.*[cC]rotamine.*"
  LAAO <- ".*L-amino.*|.*(LAAO).*"
  Cathelicidin <- ".*Cathelicidin.*"
  NUC <- ".*5'-nucleotidase.*"
  SARAFO <- ".*[bB]ibro.*|.*[sS]arafo.*"
  Cystatin <- ".*[Cc]ystatin.*"
  Hyaluronidase <- ".*[Hh]yaluronidase.*"
  NGF <- ".*(NGF).*"
  VGF <- '.*(VEGF).*'
  svMMP <- '.*(svMMP).*'
  seMMP <- ".*(seMMP).*"
  
  # vectors
  name_vector <- c("PLA2", "SVMP", "ThreeFTXs", "CTL", "CRISP", "KUNITZ", "SVSP", "BPP", "NGF", "PLB", "Waprin", "SVNP", "Myotoxin", "LAAO", "Cathelicidin", "NUC", "SARAFO", "Cystatin", "Hyaluronidase", "svMMP", "seMMP")
  sub_vector <- c(PLA2, SVMP, ThreeFTXs, CTL, CRISP, KUNITZ, SVSP, BPP, NGF, PLB, Waprin, SVNP, Myotoxin, LAAO, Cathelicidin, NUC, SARAFO, Cystatin, Hyaluronidase, svMMP, seMMP)
  tempdf <- data.frame(name_vector, sub_vector)

  # rename in a loop
  for(i in 1:nrow(tempdf)){
    
    data$geneFamily <- gsub(tempdf[i,2], tempdf[i,1], ignore.case = TRUE, perl = TRUE, data$geneFamily)
  }
  
  ## Naming the "others"
  #name_vector <- name_vector[-length(name_vector)]   # remove last 'TBD' column from the name_vector
  name_all_toxins <- paste0(name_vector, collapse = "|") # create a single grep search for all toxins
  name_vector_topToxins <- "PLA2|SVMP|ThreeFTXs|CRISP|KUNITZ|SVSP|PLB|Waprin|Myotoxin|Crystatin"
  name_vector_otherToxins <- "CTL|BPP|NGF|Myotoxin|LAAO|Cathelicidin|NUC|SARAFO|VGF|svMMP|seMMP|Hyaluronidase|SVNP"
  
  ## 
  
  # filter df by all toxin genes that were renamed
  # Non-toxin
  data_non_toxin <- data %>% 
    filter(!grepl(name_all_toxins, geneFamily)) %>%
    mutate(geneFamily = "Non-toxin")
  
  # top 20 toxins
  data_toxin <- data %>% 
    filter(grepl(name_vector_topToxins , geneFamily)) 
  
  # all other toxins
  data_toxin_other <- data %>%
    filter(grepl(name_vector_otherToxins, geneFamily)) %>%
    mutate(geneFamily = "Other")
  
  # Join to create full renamed data frame
  renamed_data <- data_non_toxin %>%
    full_join(data_toxin) %>% 
    full_join(data_toxin_other) %>%
    mutate(genes_type = case_when(geneFamily %in% c(name_vector, "Other") ~ "Toxin",
                                  !geneFamily %in% name_vector ~ "Non-toxin")) # for plot 1

  
  return(renamed_data)
  
}

## Plot 1: total FPKM for toxin and other genes
plot_1 <- function(data, taxon_ID = NULL, save_plot = FALSE, image_type = 'pdf', path_plot = NULL){
  data_FPKM <- data %>%
    filter(eValue < 1e-4, TPM > 1) %>%
    select(taxonID, seqID, FPKM) %>% 
    group_by(taxonID, seqID) %>%  # seqID lets us grab one for each sequence
    summarize(`log2FPKM+1`=log2(FPKM+1)) 
  
  plot_dat <- data_FPKM %>%
    filter(taxonID %in% taxon_ID) %>%
    arrange(-`log2FPKM+1`) %>%
    mutate(transcript_rank= 1:n()) %>% 
    full_join(data)
  
  plot_1 <- plot_dat %>%
    ggplot(aes(x= transcript_rank, y=`log2FPKM+1`, fill = genes_type)) +
    geom_bar(stat="identity", width = 1) +
    labs(x = '', y = expression(log[2] (FPKM + 1)), title = taxon_ID) +
    scale_x_continuous(limits = c(0, 1000), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
    scale_fill_manual(values=c("#D8D8D8", "#84201D")) +
    theme(strip.text = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          strip.placement = "inside",
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          legend.position = "none",
          #legend.title = element_blank(),
          #legend.text = element_text(size = 42),
          axis.text.y =  element_text(size = 50),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title = element_text(size = 64))
  
  ## Save plot or just return it
  if(isTRUE(save_plot)){
    
    ## Create output directory if it doesn't exist already
    dir.create(path = normalizePath(path_plot), showWarnings = FALSE, recursive = TRUE)
    
    ## Build file name
    n <- paste(taxon_ID, "plot1", image_type, sep = ".")
    
    ggsave(filename = n,
           device = image_type,
           path = normalizePath(path_plot),
           plot = plot_1,
           units = "mm",
           width = 420,
           height = 240)
  } else {
    return(plot_1)
  }
}

## Plot 2: total FPKM for top 50 toxin genes####
plot_2 <- function(data, taxon_ID = NULL, save_plot = FALSE, image_type = 'pdf', path_plot = NULL, path_table = NULL){
  data_FPKM <- data %>%
    filter(eValue < 1e-4, TPM > 1) %>%
    select(taxonID, seqID, FPKM) %>% 
    group_by(taxonID, seqID) %>%  # seqID lets us grab one for each sequence
    summarize(`log2FPKM+1`=log2(FPKM+1)) 
  
  # plot
  plot_dat <- data_FPKM %>%
    full_join(data) %>%
    filter(taxonID %in% taxon_ID) %>%
    filter(genes_type %in% "Toxin") %>%
    arrange(-`log2FPKM+1`) %>%
    mutate(transcript_rank= 1:n())
  
  plot_2 <- plot_dat %>%
    ggplot(aes(x= transcript_rank, y=`log2FPKM+1`, fill = geneFamily)) +
    geom_bar(stat="identity") +
    #coord_cartesian(xlim = c(0,50), ylim = c(0,20)) +
    scale_x_continuous(limits = c(0,50), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,20), expand = c(0,0)) +
    labs(x = 'Transcript rank', y =  expression(log[2] (FPKM + 1))) +
    theme(strip.text = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          strip.placement = "inside",
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          axis.text = element_text(size = 50),
          axis.title = element_text(size =64)) +
    aes(fill= geneFamily) +
    scale_fill_manual(values=plasma)  +
    facet_wrap(~taxonID, ncol = 1)
    #scale_fill_manual(values = brewcols)
  #could also use values=viridis or plasma
  
  ## Save plot or just return it
  if(isTRUE(save_plot)){
    
    ## Create output directory if it doesn't exist already
    dir.create(path = normalizePath(path_plot), showWarnings = FALSE, recursive = TRUE)
    
    ## Build file name
    n <- paste(taxon_ID, "plot2", image_type, sep = ".")
    
    ggsave(filename = n,
           device = image_type,
           path = normalizePath(path_plot),
           plot = plot_2,
           units = "mm",
           width = 420,
           height = 240)
    
    nt <- paste(taxon_ID, "plot2", "csv", sep = ".")
    
    write_csv(x = plot_dat, file = paste0(path_table, nt))
    
  } else {
    return(plot_2)
  }
}

## Plot 3: percentage expression of all toxin genes
plot_3 <- function(data, taxon_ID = NULL, save_plot = FALSE, image_type = 'pdf', path_plot = NULL, path_table = NULL){
  
  taxon_filter <- data %>%
    filter(eValue < 1e-4, TPM > 1,
           taxonID %in% taxon_ID, genes_type %in% "Toxin") 
 
  taxon_TPM <- select(taxon_filter, TPM) 
  
  sum_taxon_TPM <- sum(taxon_TPM)
  
  # plot
  plot_dat <- taxon_filter %>%
    group_by(geneFamily) %>%
    summarise(total_TPM = sum(TPM)) %>%
    mutate(percent_TPM = (total_TPM/sum_taxon_TPM)*100)
  
  plot_3 <- plot_dat %>%
    ggplot(aes(x = "", y = percent_TPM, fill = geneFamily)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    ggtitle(taxon_ID) +
    theme_void() +
    theme(plot.title = element_text(vjust = 0, hjust=0.5))+
    aes(fill= geneFamily)+
    scale_fill_manual(values=plasma) 
  #could also use values=viridis or plasma
  ## Save plot or just return it
  if(isTRUE(save_plot)){
    
    ## Create output directory if it doesn't exist already
    dir.create(path = normalizePath(path_plot), showWarnings = FALSE, recursive = TRUE)
    
    ## Build file name
    n <- paste(taxon_ID, "plot3", image_type, sep = ".")
    
    ggsave(filename = n,
           device = image_type,
           path = normalizePath(path_plot),
           plot = plot_3)
    
    nt <- paste(taxon_ID, "plot3", "csv", sep = ".")
    
    write_csv(x = plot_dat, file = paste0(path_table, nt))
    
  } else {
    return(plot_3)
  }
}


###########
### Data prep

### Load data 

### I used the following code to load in each datasheet as a single file
### I then saved as an R data object, so now you can just load in that below
### If you need to load in new data, just unhash this section
d <- list.files(path = dir_data, pattern = ".mergedExpression.csv", full.names = TRUE) %>%
  set_names(str_remove(string = basename(.), pattern = ".mergedExpression.csv")) %>%
  map(read.csv, stringsAsFactors = FALSE) %>%
  bind_rows(.id = "taxonID") 

# Prep data for renaming script
d$blastHit <- str_replace_all(d$blastHit, "[[:punct:]]", " ")
d$geneFamily <- d$blastHit

## rename gene families for the combined data
dat <- rename_genes(d)


################
### OPTIONAL
### Save and load as an rds
# save as R data object
#saveRDS(object = dat, file = paste0(dir_path,"../genecounts_renamed"))

# Load in the following r data object
#dat <- readRDS(file =  paste0(dir_path,"../genecounts_renamed") )
#################

## Explore new df with renamed toxin gene families
head(dat)
unique(dat$geneFamily)

###
## Gene count summary
# Number of toxin transcripts per sample
gene_count <- dat %>%
  filter(genes_type == "Toxin") %>%
  group_by(taxonID, geneFamily) %>%
  summarize(toxin_count = n())

print(as_tibble(gene_count), n = 50)

# save table in outputs folder
write.table(gene_count, file = paste0(dir_output, "number_toxin_transcripts.csv"), sep = ",", col.names = TRUE, row.names = FALSE)
##


### Plotting
################################################################################

##
## Save each plot for each sample
# for saving all samples
samples <- unique(dat$taxonID)

# save plots Total transcript expression (FPKM) of venom transcriptomes
for(i in samples){
  print(paste("Saving plot for ", i))
  
  plot_1(data = dat,  
         taxon_ID = i, # name of sample you want to plot
         save_plot = TRUE, # do you want to save plot to a dir?
         image_type = 'png', # file type for saving
         path_plot = dir_plot) # path to plots dir
  
}

## Save plots for plot_2 and plot_3

for(i in samples){
  print(paste("Saving plot for ", i))
  
  plot_2(data = dat,  
         taxon_ID = i, 
         save_plot = TRUE, 
         image_type = 'png', 
         path_plot = dir_plot,
         path_table = dir_output)
  
  plot_3(data = dat,  
         taxon_ID = i, 
         save_plot = TRUE, 
         image_type = 'png', 
         path_plot = dir_plot,
         path_table = dir_output)
  
}

##
sessionInfo()
##
