# bioinformatica
Alvos
#importar os ficheiros das interacções para o R

omni=read.delim("interactions_omni.txt",header=T, stringsAsFactors = F)

dorothea=read.delim("dorothea_AB.txt",header=T, stringsAsFactors = F)

#Extrair a lista de proteínas que interagem fisicamente com uma das subunidades do NFkB
# NFKBI (P105 -> P50) - P19838
# RELA (P65) - Q04206

#Extrair a lista de genes que são regulados pelo NFkB (alvos)
#Da base de dados OmniPath
linhas_alvos_omni_N=which(omni$source=="P19838")
linhas_alvos_omni_R=which(omni$source=="Q04206")

alvos_omni=union(omni$target[linhas_alvos_omni_N],omni$target[linhas_alvos_omni_R])

#Da base de dados DoRothEA
linhas_alvos_dorothea_N=which(dorothea$source=="P19838")
linhas_alvos_dorothea_R=which(dorothea$source=="Q04206")

alvos_dorothea=union(dorothea$target[linhas_alvos_dorothea_N],dorothea$target[linhas_alvos_dorothea_R])

#união entre as 2 bases de dados
alvos=union(alvos_omni,alvos_dorothea)

#Determinar quantos alvos ao todo são extraidos das bases de dados
a=length(alvos)

#Instalar as packages necessárias
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)

#Defenir o universo como estando adaptado a células HeLa
helaprot=read.delim("E-PROT-19-query-results.tsv",header=T,skip=4,stringsAsFactors = F)
helaprot=helaprot[!is.na(helaprot$HeLa),]
hela_uniprot=mapIds(org.Hs.eg.db,keys=helaprot$Gene.ID,column="UNIPROT",keytype="ENSEMBL")

#Encontrar os alvos do NF-kB nas HeLa
alvos_hela=intersect(alvos,hela_uniprot)
#Fazer o enrriquecimento
bp_ciclos_alvos_hela=enrichGO(alvos_hela,OrgDb=org.Hs.eg.db, keyType="UNIPROT",ont="BP",pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=hela_uniprot, minGSSize = 50, maxGSSize=300, readable = T)

#converter em dataframe
bp_ciclos_tab=as.data.frame(bp_ciclos_alvos_hela)
to_drop=bp_ciclos_tab$ID[bp_ciclos_tab$Count<10]
bp_ciclos_alvos_hela=dropGO(bp_ciclos_alvos_hela,term=to_drop)
bp_ciclos_tab=as.data.frame(bp_ciclos_alvos_hela)

#remover termos semelhantes
bp_ciclos_simp=simplify(bp_ciclos_alvos_hela,cutoff=0.5)
bp_ciclos_simp_tab=as.data.frame(bp_ciclos_simp)

#Instalar packages extra
BiocManager::install("DOSE")
install.packages("dplyr")
library(DOSE)
library(dplyr)

#cálculo do fold enrichment
bp_ciclos_simp =mutate(bp_ciclos_simp,FoldEnrichment = parse_ratio(bp_ciclos_simp_tab$GeneRatio) / parse_ratio(bp_ciclos_simp_tab$BgRatio))

bp_ciclos_simp_tab=as.data.frame(bp_ciclos_simp)
