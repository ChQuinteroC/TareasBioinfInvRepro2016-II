#Cargamos librerias

library(SimRAD)
library(seqinr)


#Importar genoma de referian en archvo .fasta

Laccaria <- ref.DNAseq("../data/Lacbi2_AssemblyScaffolds.fasta", subselect.contigs = FALSE)

#Exploramos el objeto, número de bases y %de GC

width(Laccaria)

GC(s2c(Laccaria))

#Enzimas, metemos los sitios de restricción de cada enzima de restricción

cut_site_5prime1 <- "CTGCA"
cut_site_3prime1 <- "G"

cut_site_5prime3 <- "C"
cut_site_3prime3 <- "CGG"

cut_site_5prime4 <- "G"
cut_site_3prime4 <- "CAGC"

#Simulación de la digestión utilizando la función insilico.digest, con la combinacion de enzimas

PstI <- insilico.digest(Laccaria, cut_site_5prime1, cut_site_3prime1, verbose = TRUE)
MspI <- insilico.digest(Laccaria, cut_site_5prime3, cut_site_3prime3, verbose = TRUE)
ApeKIA <- insilico.digest(Laccaria, cut_site_5prime4, cut_site_3prime4, verbose = TRUE)

PstI.ApKIA<- insilico.digest(Laccaria, cut_site_5prime1, cut_site_3prime1, cut_site_5prime4, cut_site_3prime4, verbose = TRUE)
PstI.MspI<- insilico.digest(Laccaria, cut_site_5prime1, cut_site_3prime1, cut_site_5prime3, cut_site_3prime3, verbose = TRUE)
MspI.ApkIA<- insilico.digest(Laccaria, cut_site_5prime3, cut_site_3prime3, cut_site_5prime4, cut_site_3prime4, verbose = TRUE)

#Adapt.select solo se utiliza cuando se usaron dos enzimas de restricción para la digestión

PstI.ApKIAadap<-adapt.select(PstI.ApKIA, type = "AB+BA",
                             cut_site_5prime1, cut_site_3prime1, cut_site_5prime4,
                             cut_site_3prime4)
PstI.MspIadap<- adapt.select(PstI.MspI, type = "AB+BA",
                             cut_site_5prime1, cut_site_3prime1, cut_site_5prime3,
                             cut_site_3prime3)
MspI.ApkIAadap<- adapt.select(MspI.ApkIA, type = "AB+BA",
                              cut_site_5prime3, cut_site_3prime3, cut_site_5prime4,
                              cut_site_3prime4)

#Selección de fragmentos


MspI.select <-size.select(MspI, min.size= 200, max.size= 400, graph = TRUE, verbose = TRUE)
PstI.select <- size.select(PstI, min.size= 200, max.size= 400, graph = TRUE, verbose = TRUE)
ApeKIA.select <- size.select(ApeKIA, min.size= 200, max.size= 400, graph = TRUE, verbose = TRUE)

PstI_ApKIAselect <- size.select(PstI.ApKIAadap, min.size= 200, max.size= 400, graph = TRUE, verbose = TRUE)
PstI_MspIselect <- size.select(PstI.MspIadap, min.size= 200, max.size= 400, graph = TRUE, verbose = TRUE)
MspI_ApkIAselect <- size.select(MspI.ApkIAadap, min.size= 200, max.size= 400, graph = TRUE, verbose = TRUE)

#Agregramos nombres a cada secuencia. A cada objeto del output de size.select se le agrega nombre a sus secuencias

names(MspI.select)<- 1:length(MspI.select)
names(PstI.select)<- 1:length(PstI.select)
names(ApeKIA.select)<- 1:length(ApeKIA.select)
names(PstI_ApKIAselect)<- 1:length(PstI_ApKIAselect)
names(PstI_MspIselect)<- 1:length(PstI_MspIselect)
names(MspI_ApkIAselect)<- 1:length(MspI_ApkIAselect)

#Los outputs de size.select son de clase caracter, por lo que debemos transformarlos a clase
#DNAstrinset para poder gurdarlos a formato FASTA

MspI.DNA <- DNAStringSet(MspI.select, start = NA, end = NA, width = NA, use.names = TRUE)
PstI.DNA <- DNAStringSet(PstI.select, start = NA, end = NA, width = NA, use.names = TRUE)
ApeKIA.DNA <- DNAStringSet(ApeKIA.select, start = NA, end = NA, width = NA, use.names = TRUE)


#Guardamos los archivos DNAStringset en formato FASTA

writeXStringSet(MspI.DNA, "MspI_sec_lb", append = FALSE, compress = FALSE, compression_level = NA, format = "fasta")
writeXStringSet(PstI.DNA, "PstI_sec_lb", append = FALSE, compress = FALSE, compression_level = NA, format = "fasta")
writeXStringSet(ApeKIA.DNA, "ApeKIA_sec_lb", append = FALSE, compress = FALSE, compression_level = NA, format = "fasta")
writeXStringSet(PstI_ApKIAselect, "PstI_ApkIA_sec_lb", append = FALSE, compress = FALSE, compression_level = NA, format = "fasta")
writeXStringSet(PstI_MspIselect, "PstI_MspI_sec_lb", append = FALSE, compress = FALSE, compression_level = NA, format = "fasta")
writeXStringSet(MspI_ApkIAselect, "MspI_ApkIA_sec_lb", append = FALSE, compress = FALSE, compression_level = NA, format = "fasta")
