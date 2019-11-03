library(readxl)
library(readr)
library(purrr)
library(stringr)
library(tidyverse)
library(ggrepel)
MS110.con <- read_excel('../20190725 YLS SILAC/20190725YLS proteinGroups_110.xlsx')
MS211.siRNA <- read_excel('../20190725 YLS SILAC/20190725 YLSproteinGroups_211.xlsx')
MS312.siRNA <- read_excel('../20190725 YLS SILAC/20190725 YLS proteinGroups_312.xlsx')
MS74.con <- read_excel('../20190725 YLS SILAC/20190725YLS proteinGroups 74.xlsx')
MS85.siRNA <- read_excel('../20190725 YLS SILAC/20190725YLS proteinGroup 85.xlsx')
MS96.siRNA <- read_excel('../20190725 YLS SILAC/20190725YLS ProteinGroup_96.xlsx')

# modify the protein id ---------------------------------------------------

extract.protein.name<- function(x){
  p <- parent.frame()
  if(is.data.frame(x)){
    name <- deparse(substitute(x))
    x$`Protein IDs` <- str_extract(x$`Protein IDs`, pattern = '(?<=\\|)\\w+(?=\\_HUMAN)')
    if(exists(name, envir = p, inherits = F)){
      assign(name, x, envir = p)
    } else {
      assign(name, x, envir = where(name))
      }
  }else if(is.character(x)){
    name <- x
    obj <- get(x, envir = p)
    obj$`Protein IDs` <- str_extract(obj$`Protein IDs`, pattern = '(?<=\\|)\\w+(?=\\_HUMAN)')
    if (exists(name, envir = p, inherits = F)) {
      assign(name, obj, envir = p)
    }else{
      assign(name, obj, envir = pryr::where(name))
    }
  } else {
    errorCondition("The type of argument should be 'character' or 'dataframe'")
  }
  invisible()
}
data.obj <- ls(pattern = 'MS*') #get the obj names
lapply(data.obj, extract.protein.name)

# normalized --------------------------------------------------------------

z_score <- function(x){
  (x - mean(x))/sd(x)
}

# remove na names ---------------------------------------------------------

library(dplyr)
row.na.rm <- function(x){
  name <- x
  x <- get(x)
  x <- filter(x, !is.na(`Protein IDs`))
  assign(name, x, envir = pryr::where(name))
  invisible()
}
lapply(data.obj, row.na.rm)

# select rows -------------------------------------------------------------

select.intensity <- function(x){
  name <- x
  x <- get(x)
  x <- select(x, `Protein IDs`,starts_with('Intensity'))
  assign(name, x, envir = pryr::where(name))
  invisible()
}

lapply(data.obj, select.intensity)

#names(MS110.con) <- str_replace_all(names(MS110.con), ' ', '.')

combinename <- function(dataframe) {
  names(dataframe) <- str_replace_all(names(dataframe), ' ', '.')
  return(dataframe)
}

for (i in seq_len(length(data.obj))) {
  name <- data.obj[i]
  temp <- combinename(get(data.obj[i]))
  temp <- select(temp, "Protein.IDs","Intensity.L", "Intensity.H")
  temp[2:3] <- map(temp[2:3], z_score)
  assign(name, temp, envir = pryr::where(name))
  rm(temp)
}

temp <- intersect(get(data.obj[1])$Protein.IDs, get(data.obj[2])$Protein.IDs)
temp <- intersect(temp, get(data.obj[3])$Protein.IDs)
temp <- intersect(temp, get(data.obj[4])$Protein.IDs)
temp <- intersect(temp, get(data.obj[5])$Protein.IDs)
temp <- intersect(temp, get(data.obj[6])$Protein.IDs)

filterbytmp <- function(obj) {
  obj <- filter(obj, Protein.IDs %in% temp)
}

for (i in seq_len(length(data.obj))) {
  name <- data.obj[i]
  obj <- get(data.obj[i])
  obj <- filterbytmp(obj) %>% 
    group_by(Protein.IDs) %>% 
    summarise(Intensity.L = mean(Intensity.L),
              Intensity.H = mean(Intensity.H))
  assign(name, obj, envir = pryr::where(name))
  rm(obj)
}

con <- inner_join(MS110.con, MS74.con, by = "Protein.IDs")
siRNA <- inner_join(MS211.siRNA, MS312.siRNA, by = "Protein.IDs")
siRNA <- inner_join(siRNA, MS85.siRNA, by = "Protein.IDs")
siRNA <- inner_join(siRNA, MS96.siRNA, by = "Protein.IDs")

files <- con$Protein.IDs
U3Cko <- con[1]
U3Cko <- add_column(U3Cko,
                    p=(1:length(files)),
                    con = (1:length(files)),
                    siRNA = (1:length(files)))
ttest <- vector('list',length(files))

for (i in seq_along(files)) {
  con.tmp <- filter(con, Protein.IDs == files[i])
  siRNA.tmp <- filter(siRNA, Protein.IDs == files[i])
  con.t <- pull(as.data.frame(t(con.tmp), stringsAsFactors = F), var = -1)[2:5]
  siRNA.t <- pull(as.data.frame(t(siRNA.tmp), stringsAsFactors = F), var = -1)[2:9]
  ttest[[i]] <- t.test(as.numeric(con.t), as.numeric(siRNA.t))
  U3Cko$p[[i]] <- ttest[[i]]$p.value
  U3Cko$con[[i]] <- mean(as.numeric(con.t))
  U3Cko$siRNA[[i]] <- mean(as.numeric(siRNA.t))
}

# plot --------------------------------------------------------------------
TUBE3CKO_ms <- U3Cko
TUBE3CKO_ms <- mutate(TUBE3CKO_ms, log2KOdividecon = log2(siRNA / con))
TUBE3CKO_ms <- mutate(TUBE3CKO_ms, log10p = -log10(p))
TUBE3CKO_ms <- mutate(TUBE3CKO_ms, KOdividecon = siRNA / con)

TUBE3CKO_ms$label <- 'total'
upgene <- filter(TUBE3CKO_ms,log10p > (-log10(0.05)) & log2KOdividecon > 1 & log2KOdividecon != Inf)
upgene$label <- 'upgene'
downgene <- filter(TUBE3CKO_ms,log10p > (-log10(0.05)) & log2KOdividecon < -1 & log2KOdividecon != -Inf)
downgene$label <- 'downgene'

TBB <- tibble(Protein.IDs = str_subset(TUBE3CKO_ms$Protein.IDs, pattern = "^TBB\\w*"))
TBBjoin <- left_join(TBB,TUBE3CKO_ms, by = "Protein.IDs")
TBBjoin$label <- 'tubulin'
ACT <- tibble(Protein.IDs = str_subset(TUBE3CKO_ms$Protein.IDs, pattern = "^ACT\\w*"))
ACTjoin <- right_join(TUBE3CKO_ms,ACT, by = "Protein.IDs")
ACTjoin$label <- 'actin'

ggplot(TUBE3CKO_ms,aes(log2KOdividecon,log10p)) + 
  geom_vline(xintercept  = 1, color = 'white', size = 1.5) +
  geom_vline(xintercept  = -1, color = 'white', size = 1.5) +
  geom_hline(yintercept = -log10(0.05), color = 'white', size = 1.5) +
  
  geom_point(aes(color = label), alpha = 0.2) +
  labs(y = '-log10(P-value)',
       x = 'log2(Fold Change)',
       title = 'The differences of protein expression in the UBE3C siRNA KO 293T cells') +
  geom_point(data = upgene, aes(log2KOdividecon,log10p,color = label) ) +
  geom_point(data = downgene, aes(log2KOdividecon,log10p, color = label)) +
  geom_point(data = TBBjoin, aes(log2KOdividecon,log10p, color = label)) +
  geom_point(data = ACTjoin, aes(log2KOdividecon,log10p, color = label)) +
  geom_text_repel(data = upgene, aes(label = Protein.IDs),color = 'red') +
  geom_text_repel(data = downgene, aes(label = Protein.IDs),color = 'blue') +
  geom_text_repel(data = ACTjoin, aes(label = Protein.IDs),color = 'green') +
  geom_text_repel(data = TBBjoin, aes(label = Protein.IDs),color = 'yellow') +
  
  scale_color_manual(name = 'Significant', 
                     limits = c('upgene','total','downgene','actin','tubulin'),
                     label = c('Up','No','Down','actin','tubulin'), 
                     values = c('red','gray20','blue','green','yellow')) +
  theme(axis.title = element_text(family = "Helvetica"), 
        axis.text = element_text(family = "Helvetica"), 
        axis.text.x = element_text(family = "Helvetica"), 
        axis.text.y = element_text(family = "Helvetica"), 
        plot.title = element_text(family = "Helvetica", hjust = 0.5), 
        legend.text = element_text(family = "Helvetica"), 
        legend.title = element_text(family = "Helvetica")) +
  guides(color = guide_legend(override.aes = list(alpha = 1,size = 3))) +
  scale_x_continuous(breaks = c(-4,-2, 0, 2, 4),
                     labels = c('-4','-2','0','2','4')) +
  scale_y_continuous(breaks = c(0,1,-log10(0.05),2),
                     labels = c('0','1','p = 0.05','2')) +
  
  ggsave(filename = 'UBEKO total protein MS.tiff', width = 15, height = 9.27)

write_csv(upgene,path = 'upgene.csv')
write_csv(downgene,path = 'downgene.csv')
write_csv(TUBE3CKO_ms,path = 'Total.csv')