## Plot circular ideogram

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(circlize)
library(dplyr)

x = rnorm(1600)
sectors = sample(c(1:22), 1600, replace = TRUE)
sectors = paste0("chr",sectors)

circos.initializeWithIdeogram(x=x)
circos.trackHist(sectors, x = x, col = "#999999", 
                 border = "#999999")
circos.trackHist(sectors, x = x, bin.size = 0.1, 
                 col = "#999999", border = "#999999")

data = as.data.frame(x)
circos.genomicDensity(data, col = c("#FF000080"), track.height = 0.1)

load(system.file(package = "circlize", "extdata", "DMR.RData"))

circos.initializeWithIdeogram()
circos.genomicDensity(DMR_hyper, col = c("#FF000080"), track.height = 0.25,border="#FFFFFFFF")

# This histogram example works
test.data = data.frame(chr = paste0("chr",sample(c(1:22), 1600, replace = TRUE)),
                       start = rnorm(1600)*1e7)
test.data$stop = test.data$start+20
circos.initializeWithIdeogram()
circos.genomicDensity(test.data, col = c("#FF000080"), track.height = 0.25,border="#FFFFFFFF")

## Import and parse Brunello guides

# Import Brunello guides
Brnlo.df = read.table("BrunelloGuidesInfo.tsv",sep="\t",header=T,stringsAsFactors=F)
Brnlo.df = Brnlo.df[!is.na(Brnlo.df$Target.Gene.ID),] # Remove NTCs

# Record the chromosome each guide aligns to

Brnlo.df$chr = gsub("NC_0*|\\.[0-9]*","",Brnlo.df$Genomic.Sequence)

# Reorder in increasing chromosome number
Brnlo.df = Brnlo.df %>%
  arrange(chr)

# Parse chromosome names
# chr23 = X
# chr24 = Y
Brnlo.df$chr[Brnlo.df$chr==23] = "X"
Brnlo.df$chr[Brnlo.df$chr==24] = "Y"

Brnlo.df$chr.name = paste0("chr",Brnlo.df$chr)


# This histogram example works
plot.hist = data.frame(chr = Brnlo.df$chr.name,
                       start = Brnlo.df$Position.of.Base.After.Cut..1.based.)
plot.hist$stop = plot.hist$start
# plot.hist$color = ifelse(plot.hist$chr%in%paste0("chr",c(1:10)),"orange","green")
plot.hist$color = viridis::viridis(24)[as.numeric(gsub("chr","",plot.hist$chr))]
plot.hist$color[is.na(plot.hist$color)] = "orange"

cols = rainbow(24)

bed = data.frame( # TO DO: change start values to middle of chromosome
  chr = c("chr1","chr5"),
  start = 0,
  end = 20,
  value = runif(2)
)

# tiff(file="BrunelloHistogram.tiff",
#      width=2000, height=2000,res=100)
circos.clear()
circos.par("start.degree"=90)
circos.initializeWithIdeogram(plotType = NULL)
circos.trackHist(plot.hist$chr, x = plot.hist$start, col = cols, 
                 border = cols,bin.size = 1e6,bg.border = "#FFFFFF",
                 track.height = 0.4) 
circos.genomicIdeogram(track.height=0.05) # put ideogram as the inside track
circos.genomicTrackPlotRegion(bed,bg.border = "#FFFFFF", panel.fun = function(region, value, ...) {
  circos.genomicText(region, value,y=1 ,labels="chr1")
})
# dev.off()
circos.clear()


circos.par("track.height" = 0.1, cell.padding = c(0, 0, 0, 0))
circos.initializeWithIdeogram(plotType = NULL)

bed = generateRandomBed(nr = 20)

circos.genomicTrack(bed, ylim = c(0, 1), panel.fun = function(region, value, ...) {
  circos.genomicText(region, value, y = 0.5, labels = "text", ...)
})

## Use "rainbow" as a color scale so that it reflects the EGFR library schematic

circos.genomicIdeogram() # put ideogram as the third track


circos.trackHist(plot.hist$chr, x = plot.hist$start, col = col_fun, 
                 border = "#FF000080",bin.size = 1e6) 

circos.genomicDensity(plot.hist, col = c("#FF000080"), track.height = 0.25,border="#FFFFFFFF")



circos.genomicLines(plot.hist$chr, plot.hist$start)



# Note chromosomes for each gene
chr.df = data.frame(gene = unique(Brnlo.df$Target.Gene.Symbol))

chr.df$chr = NA
x = getBM(attributes = "chromosome_name",
      filters = "hgnc_symbol", values = chr.df$gene,
      mart = human)
