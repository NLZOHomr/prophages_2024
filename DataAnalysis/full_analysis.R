## Script for the analysis of prophage data
library(tidyverse)
library(micropan)
library(stringr)
library(ggsci)
library(ggnested)
library(ComplexHeatmap)
library(gplots)
library(reshape2)
library(ggvenn)

# https://github.com/wanyuac/handyR/blob/master/R/phylo__makeFastANIDistMatrix.R
makeFastANIDistMatrix <- function(f, keep_asym = FALSE, frac = FALSE, suffix = "fasta") {
  # Initiation
  ani <- read.delim(file = f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1 : 3]
  names(ani) <- c("Query", "Reference", "ANI")
  ani$Query <- sapply(ani$Query, .extractSampleName, suffix)
  ani$Reference <- sapply(ani$Reference, .extractSampleName, suffix)
  ani$D <- 100 - ani$ANI  # Calculate distances from ANIs
  ani <- ani[, -3]  # Remove the column "ANI"
  
  if (frac) {
    ani$D <- ani$D / 100  # Convert percentages to decimal fractions
    precision <- 6  # Number of decimals to keep
  } else {
    precision <- 4  # The same as FastANI
  }
  
  ids <- sort(union(ani$Query, ani$Reference), decreasing = FALSE)
  n <- length(ids)
  M <- matrix(data = NA, nrow = n, ncol = n, dimnames = list(ids, ids))
  diag(M) <- 0
  
  # Stage one: copy values from the data frame to matrix M
  for (i in 1 : nrow(ani)) {
    rw <- ani[i, ]  # Extract one row from the data frame to increase the speed
    q <- rw$Query
    r <- rw$Reference
    if (r != q) {
      M[q, r] <- rw$D
    }
  }
  
  # Stage two: convert M into a symmetric matrix by taking the mean distance between every pair of genomes
  # This is the same method that FastANI uses for generating the PHYLIP-formatted lower triangular matrix.
  # See https://github.com/ParBLiSS/FastANI/issues/36
  if (keep_asym) {
    M_asym <- M
  }
  
  for (i in 1 : (n - 1)) {
    for (j in (i + 1) : n) {
      val_up <- M[i, j]  # The value in the upper triangle
      val_lo <- M[j, i]  # The value in the lower triangle
      v <- round((val_up + val_lo) / 2, digits = precision)  # The same precision as FastANI (after dividing values by 100)
      M[i, j] <- v
      M[j, i] <- v
    }
  }
  
  # Return the result
  if (keep_asym) {
    out <- list("D" = M, "D_asym" = M_asym)
  } else {
    out <- M
  }
  
  return(out)
}

.extractSampleName <- function(fasta_path, suffix) {
  fields <- unlist(strsplit(x = fasta_path, split = "/", fixed = TRUE))
  f <- fields[length(fields)]
  f <- gsub(pattern = paste0(".", suffix), replacement = "", x = f, fixed = TRUE)
  return(f)
}

####### DATA IMPORT ########

#' Mapping files:
#bacterial isolates
bacteriaMap=read_csv('~/Desktop/HQ_prophage_collection/bacteria_metadata.csv')

#prophages from different tools
#phaster
phasterMap=read_csv('~/Desktop/HQ_prophage_collection/phaster_metadata.csv')

#vibrant
vibrantMap=read_csv('~/Desktop/HQ_prophage_collection/vibrant_metadata.csv')

#cenote-taker3
ct3Map=read_csv('~/Desktop/HQ_prophage_collection/ct3_metadata.csv')

#vOTU data
vOTU=read_delim('~/Desktop/HQ_prophage_collection/mmseqs2output/mmseqs95_cluster.tsv') %>%
  mutate(header_id=gsub("MB20_", "MB20-", header_id),
         header_id=gsub("MB18_", "MB18-", header_id),
         header_id=gsub("Btheta_", "", header_id))

#' mapping files for fecal water


#' detection of vOTUs in other databases
mgv_match=read.delim("~/Desktop/HQ_prophage_collection/MGV_besthit.m8.tsv", header=F)
gpd_match=read.delim("~/Desktop/HQ_prophage_collection/GPD_besthit.m8.tsv", header=F)
gvd_match=read.delim("~/Desktop/HQ_prophage_collection/GVD_besthit.m8.tsv", header=F)

mgv_map=read.delim("~/Desktop/HQ_prophage_collection/mgv_matched_metadata.tsv", header=F)
gpd_map=read.delim("~/Desktop/HQ_prophage_collection/gpd_matched_metadata.tsv", header=F)
gvd_map=read.delim("~/Desktop/HQ_prophage_collection/gvp_matched_metadata.tsv", header=F)




#### GENERAL OVERVIEW OF BACTERIAL ISOLATES AND THEIR GENOMES #####
nrow(bacteriaMap)
bacteriaMap %>% 
  group_by(Identification) %>%
  summarise(NoGenomes=length(unique(Strain_label)))

bacteriaMap %>% 
  group_by(genus) %>%
  summarise(NoGenomes=length(unique(Strain_label)))

Fig1=ggnested(bacteriaMap, 
         aes(genus, 
             main_group = genus, 
             sub_group = Identification)) + 
  geom_bar() +
  theme_bw() +
  facet_grid(~person+time, scales = 'free_x', space = 'free') +
  labs(x=NULL, y='Number of sequenced isolates') +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,17),
                     breaks = seq(0,17,2))

####' fastANI with targeted bacterial genomes
x <- read.table("~/Desktop/HQ_prophage_collection/fastani_wgs_output.txt")
matrix <- acast(x, V1~V2, value.var="V3")
matrix[is.na(matrix)] <- 70

### define the colors within 2 zones
breaks = seq(min(matrix), max(100), length.out=100)
gradient1 = colorpanel( sum( breaks[-1]<=95), "red", "white" )
gradient2 = colorpanel( sum( breaks[-1]>95 & breaks[-1]<=100), "white", "darkblue" )
hm.colors = c(gradient1, gradient2)

#hm.colors = c(gradient1, gradient2)
Fig2=heatmap.2(matrix, scale = "none", trace = "none", 
          #col = gradient2,
          col = hm.colors,
          cexRow=.30, cexCol=.30)



##### PROPHAGES: DISCOVERY #########
bact_col <- c(
  'Bacteroides faecis' = "#00355B",
  'Bacteroides ovatus' = "#20557B",
  'Bacteroides thetaiotamicron' = "#6297BD",
  'Bacteroides uniformis' = "#83B8DE",
  'Phocaeicola massiliensis' = "#B24A00",
  'Phocaeicola vulgatus' = "#FF974D"
)

vibrantMap_update=vibrantMap %>%
  mutate(scaffold_new=str_split(pattern = '_fragment', scaffold, simplify = T)[,1])

phasterMap_update= phasterMap %>% 
  mutate(scaffold_new=str_split(pattern = ':',contig, simplify = T)[,1])

#' prophage quality
Pmap=phasterMap %>%
  left_join(bacteriaMap) %>%
  group_by(Identification, contig, quality, newName) %>%
  summarise(Identification=unique(Identification),
            prophageID=unique(sample),
            quality=unique(quality),
            name=unique(newName),
            scaffold_new=str_split(pattern = ':',contig, simplify = T)[,1]) %>%
  mutate(tool='PHASTER')

Vmap=vibrantMap %>%
  left_join(bacteriaMap) %>%
  group_by(Identification, Quality, prophageID,scaffold) %>%
  summarise(Identification=unique(Identification),
            quality=unique(Quality),
            prophageID=unique(prophageID),
            name=unique(newName)) %>%
  ungroup() %>%
  select(-Quality) %>%
  mutate(tool='VIBRANT',
         scaffold_new=str_split(pattern = '_fragment', scaffold, simplify = T)[,1])

Cmap=ct3Map %>%
  group_by(Identification,contig,input_name) %>%
  summarise(Identification=unique(Identification),
            prophageID=unique(contig),
            scaffold_new=unique(input_name),
            name=unique(newName)
            ) %>%
  mutate(quality='complete',
         tool='Cenote-Taker3') 

allMap=rbind(Pmap[c(1,3,8)], Vmap[c(1,4,6)], Cmap[c(1,7,8)])
allMap %>% group_by(Identification) %>%
  summarise(NoGenomes=length(quality))
  
Fig3=allMap %>% mutate(
  quality=factor(quality, levels = c('complete','complete circular','intact',
                                     'high quality draft','medium quality draft',
                                     'questionable','incomplete', 'low quality draft')
  )) %>%
  ggplot(aes(x=quality,fill=Identification)) +
  theme_bw() +
  geom_bar() +
  scale_fill_manual(values = bact_col) +

  facet_grid(~tool, scales = 'free_x', space = 'free')+
#  theme(strip.background=element_rect(fill=facet_colors))+
  labs(x=NULL, y='Number of predicted prophages', fill=NULL)


#' VENN diagram 

x <- list(
  PHASTER = phasterMap_update$scaffold_new, 
  VIBRANT = vibrantMap_update$scaffold_new, 
  CENOTETAKER3 = ct3Map$input_name
)

Fig4=ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)



##### PROPHAGES: DYNAMICS #########

names(vOTU)[2]='prophageID' 
Cmap=Cmap %>% mutate(prophageID=gsub("MB20_", "MB20-", prophageID),
                     prophageID=gsub("MB18_", "MB18-", prophageID),
                     prophageID=gsub("Btheta_", "", prophageID))

vOTUs_v=vOTU %>% filter(tool =='vibrant') %>% left_join(Vmap[-6], by=c('prophageID'='scaffold'))
vOTUs_c=vOTU %>% filter(tool =='ct3') %>% left_join(Cmap[-c(2,8)], by='prophageID')
vOTUs_p=vOTU %>% filter(tool =='phaster') %>% left_join(Pmap[-c(4,8)], by='prophageID')

# VENN diagram with only HQ as determined by checkM
x_2 <- list(
  PHASTER = phasterMap_update$scaffold_new[phasterMap_update$scaffold_new %in% vOTUs_p$scaffold_new], 
  VIBRANT = unique(vibrantMap_update$scaffold_new[vibrantMap_update$scaffold_new %in% vOTUs_v$scaffold_new]), 
  CENOTETAKER3 = ct3Map$input_name[ct3Map$input_name %in% vOTUs_c$scaffold_new]
)

Fig4_2=ggvenn(
  x_2, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)


head(vOTUs_p)
vOTUs_p=vOTUs_p[-c(6,7)]

head(vOTUs_v)
vOTUs_v=vOTUs_v[-c(2,7)]
names(vOTUs_v)[5]="prophageID"

vOTUs_c=vOTUs_c[-c(6,9)]

vOTU_all=rbind(vOTUs_p, vOTUs_c, vOTUs_v) %>% left_join(bacteriaMap[-2], by=c('name'='newName'))

vOTU_all %>%
  group_by(person, time, genus, vOTU) %>%
  summarise(nProphages=length(prophageID)) %>%
  group_by(vOTU, genus, person, time) %>%
  summarise(nProphages=sum(nProphages)) %>%
  filter(!is.na(vOTU)) %>%
  arrange(desc(nProphages)) %>%
  ggplot(aes(x=factor(vOTU), 
             size=nProphages, 
             col=genus, y=genus)) +
  theme_bw()+
  geom_point() +
  scale_color_uchicago()+
  labs(#title = 'Change in prophage diversity over time, subject and host range',
    #subtitle='vOTUs defined by mmseqs2 (95% similarity)',
    x=NULL, y=NULL,
    size='# of prophages\nwithin a cluster') +
  facet_grid(time~person, scales = 'free',space = 'free_x')+
  guides(color='none') +
  theme(axis.text.x = element_text(angle=45, hjust = 1))


df_mtx=vOTU_all%>% select("prophageID","vOTU") %>% mutate(No=1) %>%
  pivot_wider(names_from = prophageID, values_from = No, values_fill = 0, values_fn=sum)

make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}

mtx_votus_v2=make_matrix(select(df_mtx,-vOTU),
                         pull(df_mtx,vOTU))

clust_votus <- hclust(dist(mtx_votus_v2))

tmp4_v2=vOTU_all %>%
  dplyr::group_by(vOTU, Strain_label, Identification) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

tmp4_order_v2=tmp4_v2[order(tmp4_v2$Identification),]
tmp4_order_v2%>%group_by(vOTU) %>%
  mutate(vOTU_no=length(unique(Identification))) %>%
  arrange(desc(vOTU_no)) %>%
  filter(vOTU_no>1)
# there is 1 vOTU that is found in more than one genus. 
# There are no vOTUs that would be detected in different species within a genus.

# (Fig5=ggplot(tmp4_order_v2, 
#                    aes(vOTU, factor(Strain_label, levels=unique(Strain_label)),
#                        fill = Identification)) +
#     geom_tile() +
#     scale_fill_manual(values = bact_col) +
#     coord_fixed()+
#     #coord_flip()+
#     theme_bw() +
#     geom_vline(xintercept=2, col="grey60",linetype='dashed')+
#     scale_x_discrete(limits = rownames(mtx_votus_v2)[clust_votus$order])+
#     labs(y='Bacterial genomes', x='vOTUs', fill=NULL, alpha="No. of vOTUs")
#   +
#     theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8))
#   #    axis.text.x=element_blank()))
# )

# creating a fastANI figure of WGS where prophages are derived from 
x <- read.table("~/Desktop/HQ_prophage_collection/fastani_wgs_output.txt")

mergedTmp=tmp4_order_v2 %>% left_join(bacteriaMap)
includedWGS=unique(mergedTmp$newName)

subset_x <- x[x$V1  %in% includedWGS, ]
subset_x_2 <- subset_x[subset_x$V2 %in% includedWGS, ]

matrix_filt <- acast(subset_x_2, V1~V2, value.var="V3")
matrix_filt[is.na(matrix_filt)] <- 70

### define the colors within 2 zones
breaks = seq(min(matrix_filt), max(100), length.out=100)
gradient1 = colorpanel( sum( breaks[-1]<=95), "red", "white" )
gradient2 = colorpanel( sum( breaks[-1]>95 & breaks[-1]<=100), "white", "darkblue" )
hm.colors = c(gradient1, gradient2)

Fig6=heatmap.2(matrix_filt, scale = "none", trace = "none", dendrogram="row",
               col = hm.colors,
               cexRow=.30, cexCol=.30,)

plot(Fig6$rowDendrogram)

#Use the ordering of the fastANI plot for the prophage host figure:
ordered_axis=rownames(matrix_filt)[Fig6$rowInd]

(Fig5_v2=ggplot(mergedTmp, 
             aes(vOTU, factor(newName, levels=ordered_axis),
                 fill = Identification)) +
    geom_tile() +
    scale_fill_manual(values = bact_col) +
    coord_fixed()+
    #coord_flip()+
    theme_bw() +
    geom_vline(xintercept=2, col="grey60",linetype='dashed')+
    scale_x_discrete(limits = rownames(mtx_votus_v2)[clust_votus$order])+
    labs(y='Bacterial genomes', x='vOTUs', fill=NULL, alpha="No. of vOTUs")
  +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8))
)


#### Prophages over time

vOTU_all$name <- factor(vOTU_all$name, levels = unique(vOTU_all$name[order(vOTU_all$Identification)]))

Fig7=vOTU_all %>%
    ggplot(aes(x = factor(vOTU), col = Identification, y = reorder(name, Identification))) +
    geom_point() +
    scale_color_manual(values = bact_col) +  # Assuming bact_col is defined elsewhere
    labs(x = NULL, y = NULL) +
    facet_grid(time ~ person, scales = 'free', space = 'free_x') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  


###### PROPHAGES IN FECAL WATER SAMPLES ######
#' coverage calculated using coverM 
coverage_df=read_delim("~/Desktop/HQ_prophage_collection/prophage_metag.cov") 

v_map=vOTU_all %>% filter(tool == 'vibrant') %>% 
  left_join(Vmap[-c(4,6)], by=c('scaffold_new', "Identification", "name", "prophageID"))

names(vOTU_all)
names(v_map)
v_map=v_map[-2]
names(v_map)[12]='prophageID'

tmp_map=vOTU_all %>% filter(tool!='vibrant') %>% rbind(v_map)

coverage_df$Contig %in% tmp_map$prophageID

coverage_full=left_join(coverage_df,tmp_map, by=c('Contig'='prophageID'))%>%
  pivot_longer(cols = starts_with("D"),
               names_to = "sample",
               values_to = "TPM")

Fig5=coverage_full %>%
  filter(TPM>0) %>%
  mutate(TPM=if_else(TPM==0, NA, TPM)) %>%
  ggplot(aes(vOTU, sample, col=sample, size=TPM)) +
  geom_point() +
  theme_bw() +
  facet_grid(~person, scales = 'free') +
  scale_color_futurama() +
  labs(x=NULL, y=NULL) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  guides(color="none")


###### PROPHAGES IN OTHER DATASETS ######
gvd_match=gvd_match[c(1,2)]
names(gvd_match)=c('prophageID', 'ref')
gvd=left_join(gvd_match, gvd_map, by=c('ref' = 'V1'))

gpd_match=gpd_match[c(1,2)]
names(gpd_match)=c('prophageID', 'ref')
gpd=left_join(gpd_match, gpd_map, by=c('ref' = 'V1'))

mgv_match=mgv_match[c(1,2)]
names(mgv_match)=c('prophageID', 'ref')
mgv=left_join(mgv_match, mgv_map, by=c('ref' = 'V1'))

######

# Summarizing genomic information of identified prophages

checkV_summ=read_csv("OneDrive - Ministrstvo za javno upravo - zdravstvo-65868816-Nacionalni laboratorij za zdravje, okolje in hrano/HQ_prophage_collection/checkV_summary_prophages.csv") %>%
  mutate(contig_id=gsub("Btheta_", "", contig_id),
         contig_id=str_split(contig_id, '_intact', simplify = T)[,1],
         contig_id=gsub("MB20-", "MB20_", contig_id),
         contig_id=gsub("MB18-", "MB18_", contig_id))

vOTU_v2= vOTU %>%
  mutate(prophageID=gsub("MB20-", "MB20_", prophageID),
         prophageID=gsub("MB18-", "MB18_", prophageID))

prophage_stats=left_join(checkV_summ,vOTU_v2, by=c("contig_id"="prophageID"))

summary_prophage_quality=prophage_stats %>%
  group_by(vOTU, tool.x) %>%
  summarise(nOTUs=length(contig_id),
            meanLength=mean(contig_length),
            maxLength=max(contig_length),
            minLenght=min(contig_length),
            meanORF=mean(gene_count),
            maxORF=max(gene_count),
            minORF=min(gene_count),
            #vibrant=if_else(unique(tool.x) == 'vibrant', 1, 0),
            #ct3=if_else(unique(tool.x) == 'ct3', 1, 0),
            #phaster=if_else(unique(tool.x) == 'phaster', 1, 0)
            )

#write.csv(summary_prophage_quality, 'Desktop/HQ_prophage_collection/summary_stats_high_quality_prophages.csv')


#########

# How many prophage genes are annotated?

df <- readxl::read_xlsx('OneDrive - Ministrstvo za javno upravo - zdravstvo-65868816-Nacionalni laboratorij za zdravje, okolje in hrano/HQ_prophage_collection/manuscript/BMCmicrobiology_submissionn/Supplemental Table S2.xlsx')

df=df%>% column_to_rownames('ORF')
df=df[8:30]

# Create a new column 'is_annotated' that checks if there is any annotation in any column
df$is_annotated <- apply(df[, -1], 1, function(x) any(!is.na(x)))

# Count annotated and non-annotated genes
annotated_count <- sum(df$is_annotated)
non_annotated_count <- sum(!df$is_annotated)

cat("Number of annotated genes:", annotated_count, "\n")
cat("Number of non-annotated genes:", non_annotated_count, "\n")

anno_per=annotated_count/(annotated_count+non_annotated_count)
non_anno_per=non_annotated_count/(annotated_count+non_annotated_count)




