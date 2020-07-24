library(RMySQL)
library(tidyverse)
library(gt23)
library(GenomicRanges)
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
source('lib.R')

compPatients <- c('p03712-11', 'p04409-17', 'p03712-40')

if(! file.exists('intSites.rds')){
  
  tab <- read.table('JCI130144_tableS1.txt', header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  
  # Retrieve intSites for subjects in JCI paper as well as additional subjects for comparison. 
  p <- unique(c(compPatients, tab$Patient))
  dbConn <- dbConnect(MySQL(), group='specimen_management')
  sampleData <- dbGetQuery(dbConn, paste0('select * from gtsp where Patient in (', 
                                           paste0(sQuote(p), collapse = ', '), ')'))
  dbDisconnect(dbConn)

  frags <- getDBgenomicFragments(sampleData$SpecimenAccNum, 'specimen_management', 'intsites_miseq') 
  intSites <- unlist(GRangesList(lapply(split(frags, frags$patient), 
                                        function(x){ 
                                          stdIntSiteFragments(x) %>%
                                          collapseReplicatesCalcAbunds() %>%
                                          annotateIntSites()
                                        })))
  intSites$VCN <- sampleData[match(intSites$GTSP, sampleData$SpecimenAccNum),]$VCN
  intSites$disease <- tab[match(intSites$patient, tab$Patient),]$Disease
  intSites$trial <- tab[match(intSites$patient, tab$Patient),]$Trial
  intSites$response <- tab[match(intSites$patient, tab$Patient),]$Response
  intSites$responseClass <- tab[match(intSites$patient, tab$Patient),]$RespClass
  saveRDS(intSites, file = 'intSites.rds')
} else {
  intSites <- readRDS('intSites.rds')
}

intSites$cellType <- sub('T\\-Cells', 'T cells', intSites$cellType)
compSites <- subset(intSites, patient %in% compPatients)

intSites <- subset(intSites, disease == 'CLL')
intSites$response <- sub('CompleteWithRelapse', 'CompleteWithRelapse', intSites$response)
intSites <- subset(intSites, response %in% c('Complete', 'None'))



JCI_d0 <- subset(intSites, ! patient %in% compPatients & timePoint == 'D0' & cellType == 'T cells')


# Loop through the D28 JCI paper samples and return PBMC samples where available.
# Return Whole blood if PBMC is not available. 
o <- subset(intSites, ! patient %in% compPatients & timePoint == 'D28')

JCI_d28 <- unlist(GRangesList(lapply(split(o, o$patient), function(x){
  if('PBMC' %in% unique(x$cellType)){
    return(subset(x, cellType == 'PBMC'))
  } else {
    return(subset(x, cellType == 'Whole Blood'))
  }
})))

# Only include subjects that have both a d0 and d28 timepoint
i <- base::intersect(unique(JCI_d0$patient), unique(JCI_d28$patient))
JCI_d0 <- subset(JCI_d0, patient %in% i)
JCI_d28 <- subset(JCI_d28, patient %in% i)



comp_d0 <- subset(compSites, timePoint == 'D0' & cellType == 'T cells')

o <- subset(compSites, timePoint == 'D28')

comp_d28 <- unlist(GRangesList(lapply(split(o, o$patient), function(x){
  if('PBMC' %in% unique(x$cellType)){
    return(subset(x, cellType == 'PBMC'))
  } else {
    return(subset(x, cellType == 'Whole Blood'))
  }
})))

# Group sites by timepoint and response using the patient field which 
# is the default separator for the heatmap software.
JCI_d0$patient.org  <- JCI_d0$patient
JCI_d0$patient      <- paste('d0', JCI_d0$response)
JCI_d28$patient.org <- JCI_d28$patient
JCI_d28$patient     <- paste('d28', JCI_d28$response)

comp_d0$patient.org <- comp_d0$patient
comp_d0$patient <- paste('d0', comp_d0$patient)
comp_d28$patient.org <- comp_d28$patient
comp_d28$patient <- paste('d28', comp_d28$patient)

compSiteTable <- 
  group_by(data.frame(c(comp_d0, comp_d28)), patient, cellType) %>% 
  summarise(nSites = n_distinct(posid))


a <- bind_rows(data.frame(comp_d0), data.frame(JCI_d0))
b <- bind_rows(data.frame(comp_d28), data.frame(JCI_d28))

if(! dir.exists('genomicHeatMapData_d0')) createGenomicHeatMapData(a, outputDir = 'genomicHeatMapData_d0')
if(! dir.exists('epiGenomicHeatMapData_d0')) createEpiGenomicHeatMapData(a, outputDir = 'epiGenomicHeatMapData_d0')

genomicHeatMap_d0 <- createGenomicHeatMapPlot('genomicHeatMapData_d0', sampleOrder = unique(a$patient))
epiGenomicHeatMap_d0 <- createGenomicHeatMapPlot('epiGenomicHeatMapData_d0', sampleOrder = unique(a$patient))
ggsave(genomicHeatMap_d0$gen_plot, file = 'genomicHeatMap_d0.pdf', units = 'in', height = 8, width = 6)
ggsave(epiGenomicHeatMap_d0$gen_plot, file = 'epiGenomicHeatMap_d0.pdf', units = 'in', height = 12, width = 6)

if(! dir.exists('genomicHeatMapData_d28')) createGenomicHeatMapData(b, outputDir = 'genomicHeatMapData_d28')
if(! dir.exists('epiGenomicHeatMapData_d28')) createEpiGenomicHeatMapData(b, outputDir = 'epiGenomicHeatMapData_d28')

genomicHeatMap_d28 <- createGenomicHeatMapPlot('genomicHeatMapData_d28', sampleOrder = unique(b$patient))
epiGenomicHeatMap_d28 <- createGenomicHeatMapPlot('epiGenomicHeatMapData_d28', sampleOrder = unique(b$patient))
ggsave(genomicHeatMap_d28$gen_plot, file = 'genomicHeatMap_d28.pdf', units = 'in', height = 8, width = 6)
ggsave(epiGenomicHeatMap_d28$gen_plot, file = 'epiGenomicHeatMap_d28.pdf', units = 'in', height = 12, width = 6)


# Restore patient IDs.
comp_d0$patient  <- comp_d0$patient.org
comp_d28$patient <- comp_d28$patient.org
JCI_d0$patient   <- paste('JCI', JCI_d0$patient.org)
JCI_d28$patient  <- paste('JCI', JCI_d28$patient.org)

d <- c(comp_d0, comp_d28, JCI_d0, JCI_d28)

promoters <- bind_rows(lapply(list.files(pattern = '^Brd'), 
                    function(x){ 
                      load(x)
                      bind_rows(lapply(split(d, paste(d$patient, d$timePoint)), function(p){
                        n <- distanceToNearest(p, epigenData)
                        tibble(patient = p$patient[1], timePoint = p$timePoint[1], response = p$response[1],
                               promoter = sub('_promoters.RData', '', x),
                               p50KB = n_distinct(queryHits(n[mcols(n)$distance <= 50000,]))/n_distinct(p$posid),
                               p25KB = n_distinct(queryHits(n[mcols(n)$distance <= 25000,]))/n_distinct(p$posid),
                               p10KB = n_distinct(queryHits(n[mcols(n)$distance <= 10000,]))/n_distinct(p$posid))
                       }))
                    }))



d <- bind_rows(lapply(split(promoters, paste(promoters$promoter, promoters$patient)), function(x){
       tibble(promoter = x$promoter[1], patient = x$patient[1], response = x$response[1],
              x50 = x[x$timePoint == 'D28',]$p50KB, y50 = x[x$timePoint == 'D0',]$p50KB,
              x25 = x[x$timePoint == 'D28',]$p50KB, y25 = x[x$timePoint == 'D0',]$p25KB,
              x10 = x[x$timePoint == 'D28',]$p50KB, y10 = x[x$timePoint == 'D0',]$p10KB)
    }))
         
d$p <- ifelse(grepl('^JCI', d$patient), 'JCI patient', d$patient)
d[is.na(d$response),]$response <- 'None'

write.table(d, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, file = 'BrdPromoterRelProxData.tsv')


p <- ggplot(d, aes(x50, y50, color = p, shape = response)) +
     theme_bw()+
     labs(x = 'day 28', y = 'day 0') +
     scale_x_continuous(labels = scales::percent) +
     scale_y_continuous(labels = scales::percent) +
     scale_color_manual(name = 'Subject', values = c('gray75', 'red2', 'blue', 'green3')) +
     geom_point(size = 3) +
     facet_grid(~promoter) +
  theme(text = element_text(size=18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(p, file = 'BrdPromoterRelProx50.pdf')


