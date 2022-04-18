
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file

R_pipeline <- function (input.file="", output.file="", population="", distance.max="", include.QTL="", include.RGB="", network="", significance.threshold="", restart="", highlight.top="", placeholder="", ...){
	
	sT <- Sys.time()
	
	# read input file
	data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
	
	if(significance.threshold == "NULL"){
		significance.threshold <- NULL
	}else{
		significance.threshold <- as.numeric(significance.threshold)
	}
	
	distance.max <- as.numeric(distance.max)
	restart <- as.numeric(restart)
	highlight.top <- as.numeric(highlight.top)
	
	if(population=="NA"){
		LD.customised <- NULL
	}else{
		LD.customised <- oRDS(str_c("GWAS_LD.", population), placeholder=placeholder) %>% as.data.frame()
	}
	
	if(include.QTL=="QTL_all"){
		include.QTL <- c("eQTL_eQTLGen","pQTL_Plasma")
	}else if(include.QTL=="NA"){
		include.QTL <- NA
	}
	
	if(include.RGB=="NA"){
		include.RGB <- NA
	}
	
	#GR.SNP <- oRDS("dbSNP_GWAS", placeholder=placeholder)
	GR.SNP <- oRDS("dbSNP_Common", placeholder=placeholder)
	GR.Gene <- oRDS("UCSC_knownGene", placeholder=placeholder)

	network.customised <- oDefineNet(network=network, STRING.only=c("experimental_score","database_score"), placeholder=placeholder) %>% oNetInduce(nodes_query=names(GR.Gene), knn=0, largest.comp=F)
	
	message(sprintf("Customised the network with %d nodes and %d edges (%s) ...", gorder(network.customised), gsize(network.customised), as.character(Sys.time())), appendLF=TRUE)
	
	message(sprintf("Preparing predictors (%s) ...", as.character(Sys.time())), appendLF=TRUE)
	
	ls_pNode_genomic <- oPierSNPsAdv(data, score.cap=100, LD.customised=LD.customised, significance.threshold=significance.threshold, distance.max=distance.max, decay.kernel="constant", GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.TAD=NULL, include.QTL=include.QTL, include.RGB=include.RGB, network.customised=network.customised, placeholder=placeholder, verbose.details=F)
	
	message(sprintf("Doing prioritisation (%s) ...", as.character(Sys.time())), appendLF=TRUE)
	
	dTarget <- oPierMatrix(ls_pNode_genomic, displayBy="pvalue", aggregateBy="fishers", placeholder=placeholder)

	if(class(dTarget)=="dTarget"){
		# *_priority.txt
		dTarget$priority %>% as_tibble() %>% transmute(Gene=name, Rank=rank, Rating=signif(rating,digits=4), Type=ifelse(seed=="Y","Core","Peripheral"), Proximity=nGene, QTL=eGene, PCHiC=cGene, Description=description) %>% write_delim(output.file, delim="\t")
		
		# *_priority.xlsx
		output.file.priority <- gsub(".txt$", "_priority.xlsx", output.file, perl=T)
		dTarget$priority %>% as_tibble() %>% transmute(Gene=name, Rank=rank, Rating=signif(rating,digits=4), Type=ifelse(seed=="Y","Core","Peripheral"), Proximity=nGene, QTL=eGene, PCHiC=cGene, Description=description) %>% openxlsx::write.xlsx(output.file.priority)
		
		##########
		# *_evidence.xlsx
		df_Gene2SNP <- oPierEvidence(ls_pNode_genomic) %>% as_tibble() %>% arrange(Flag,-Pval) %>% transmute(Gene,SNP,SNP_type=ifelse(Flag=="Lead","Input","LD"),Evidence=Context) %>% distinct() %>% mutate(Evidence=str_replace_all(Evidence,"nGene_","Proximity_")) %>% mutate(Evidence=str_replace_all(Evidence,"_constant","bp")) %>% mutate(Evidence=str_replace_all(Evidence,"eGene_","QTL_")) %>% mutate(Evidence=str_replace_all(Evidence,"cGene_",""))
		df_evidence <- dTarget$priority %>% as_tibble() %>% filter(seed=="Y") %>% transmute(Gene=name, Rating=signif(rating,digits=4), Proximity=nGene, QTL=eGene, PCHiC=cGene) %>% left_join(df_Gene2SNP, by="Gene") %>% arrange(-Rating, SNP_type, SNP, Evidence)
		output.file.evidence <- gsub(".txt$", "_evidence.xlsx", output.file, perl=T)
		df_evidence %>% openxlsx::write.xlsx(output.file.evidence)
		##########
		
		# Manhattan plot
		message(sprintf("Drawing manhattan plot (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		top.label.query <- dTarget$priority %>% top_frac(1,rating) %>% pull(name)
		label.query.only <- T
		gp_manhattan_dTarget <- oPierManhattan(dTarget, color=c("darkred", "#2071b8"), point.size=0.2, top=highlight.top, top.label.type="text", top.label.size=2.5, top.label.col="black", top.label.query=top.label.query, label.query.only=label.query.only, y.lab="Priority rating (scored 0-5)", signature=F, placeholder=placeholder)
		output.file.manhattan.pdf <- gsub(".txt$", "_manhattan.pdf", output.file, perl=T)
		ggsave(output.file.manhattan.pdf, gp_manhattan_dTarget, device=cairo_pdf, width=9, height=3)
		output.file.manhattan.png <- gsub(".txt$", "_manhattan.png", output.file, perl=T)
		ggsave(output.file.manhattan.png, gp_manhattan_dTarget, type="cairo", width=9, height=3)
		
		######################################
		# RMD
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$runtime <- str_c(runtime," seconds")
		ls_rmd$data_input <- data %>% set_names(c("SNP","PVALUE"))
		ls_rmd$highlight_top <- highlight.top
		
		ls_rmd$xlsx_priority <- gsub("public/", "", output.file.priority, perl=T)
		ls_rmd$xlsx_evidence <- gsub("public/", "", output.file.evidence, perl=T)
		ls_rmd$pdf_manhattan <- gsub("public/", "", output.file.manhattan.pdf, perl=T)
		ls_rmd$png_manhattan <- gsub("public/", "", output.file.manhattan.png, perl=T)
		
		output_dir <- gsub("cTGene.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_cTGene.Rmd", bookdown::html_document2(number_sections=F,theme="readable", hightlight="default"), output_dir=output_dir)
	}
	
	##########################################
}

startT <- Sys.time()

library(tidyverse)
library(igraph)
library(GenomicRanges)

# galahad
vec <- list.files(path='/home/hfang/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# huawei
vec <- list.files(path='/root/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# mac
vec <- list.files(path='/Users/hfang/Sites/XGR/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))

R_pipeline(input.file="public/tmp/data.SNPs.STRING_high.74228649.txt", output.file="public/tmp/cTGene.SNPs.STRING_high.74228649.txt", population="EUR", distance.max="20000", include.QTL="pQTL_Plasma", include.RGB="PCHiC_PMID27863249_Monocytes", network="STRING_high", significance.threshold="5e-8", restart="0.7", highlight.top="30", placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('
--- cTGene: ',runTime,' secs ---
'), appendLF=TRUE)
