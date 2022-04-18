
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
#/home/hfang/R-3.6.2/lib/R/library
# rm -rf /home/hfang/R-3.6.2/lib/R/library/00*
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file

R_pipeline <- function (input.file="", output.file="", population="", crosslink="", network="", significance.threshold="", restart="", highlight.top="", placeholder="", ...){
	
	sT <- Sys.time()
	
	# read input file
	data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
	
	if(significance.threshold == "NULL"){
		significance.threshold <- NULL
	}else{
		significance.threshold <- as.numeric(significance.threshold)
	}
	
	restart <- as.numeric(restart)
	highlight.top <- as.numeric(highlight.top)
	
	if(population=="NA"){
		LD.customised <- NULL
	}else{
		LD.customised <- oRDS(str_c("GWAS_LD.", population), placeholder=placeholder) %>% as.data.frame()
	}
	
	relative.importance <- c(0,0,0)
	include.QTL <- NULL
	include.RGB <- NULL
	if(str_detect(crosslink, "proximity")){
		relative.importance <- c(1,0,0)
		distance.max <- str_replace_all(crosslink, "proximity_", "") %>% as.numeric()
	}else if(str_detect(crosslink, "QTL")){
		relative.importance <- c(0,1,0)
		include.QTL <- crosslink
	}else if(str_detect(crosslink, "PCHiC")){
		relative.importance <- c(0,0,1)
		include.RGB <- crosslink
	}
	
	#GR.SNP <- oRDS("dbSNP_GWAS", placeholder=placeholder)
	GR.SNP <- oRDS("dbSNP_Common", placeholder=placeholder)
	GR.Gene <- oRDS("UCSC_knownGene", placeholder=placeholder)

	network.customised <- oDefineNet(network=network, STRING.only=c("experimental_score","database_score"), placeholder=placeholder) %>% oNetInduce(nodes_query=names(GR.Gene), knn=0, largest.comp=F)
	
	message(sprintf("Customised the network with %d nodes and %d edges (%s) ...", gorder(network.customised), gsize(network.customised), as.character(Sys.time())), appendLF=TRUE)
	
	pNode <- oPierSNPs(data, score.cap=100, LD.customised=LD.customised, significance.threshold=significance.threshold, distance.max=distance.max, decay.kernel="constant", GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.QTL=include.QTL, include.RGB=include.RGB, relative.importance=relative.importance, network.customised=network.customised, placeholder=placeholder)

	if(class(pNode)=="pNode"){
		# *_CG.txt
		df_scores <- pNode$priority %>% as_tibble() %>% filter(seed==1) %>% transmute(Gene=name, Scores=weight, Description=description, Evidence=crosslink) %>% mutate(Evidence=ifelse(str_detect(Evidence,"proximity"), str_c(Evidence,"bp"), Evidence)) %>% mutate(Evidence=ifelse(str_detect(Evidence,"QTL"), str_c("QTL_",Evidence), Evidence)) %>% mutate(Evidence=str_replace_all(Evidence,"proximity","Proximity"))
		df_scores %>% write_delim(output.file, delim="\t")
		
		# *_CG.xlsx
		output.file.CG <- gsub(".txt$", "_CG.xlsx", output.file, perl=T)
		df_scores %>% openxlsx::write.xlsx(output.file.CG)
		#df_scores %>% openxlsx::write.xlsx("/Users/hfang/Sites/XGR/PiER-site/app/examples/eV2CG.xlsx")
		
		##########
		# *_CG_evidence.xlsx
		df_Gene2SNP <- oPierEvidence(pNode) %>% as_tibble() %>% arrange(Flag,-Pval) %>% transmute(Gene,SNP,SNP_type=ifelse(Flag=="Lead","Input","LD"),Evidence=crosslink) %>% distinct() %>% mutate(Evidence=ifelse(str_detect(Evidence,"proximity"), str_c(Evidence,"bp"), Evidence)) %>% mutate(Evidence=ifelse(str_detect(Evidence,"QTL"), str_c("QTL_",Evidence), Evidence)) %>% mutate(Evidence=str_replace_all(Evidence,"proximity","Proximity"))
		df_evidence <- df_scores %>% select(-Evidence) %>% left_join(df_Gene2SNP, by="Gene") %>% arrange(-Scores, SNP_type, SNP, Evidence)
		output.file.CG_evidence <- gsub(".txt$", "_CG_evidence.xlsx", output.file, perl=T)
		df_evidence %>% openxlsx::write.xlsx(output.file.CG_evidence)
		#df_evidence %>% openxlsx::write.xlsx("/Users/hfang/Sites/XGR/PiER-site/app/examples/eV2CG_evidence.xlsx")
		##########
		
		# Manhattan plot
		message(sprintf("Drawing manhattan plot (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		pNode$priority$priority <- pNode$priority$weight
		top.label.query <- pNode$priority %>% filter(seed==1) %>% pull(name)
		label.query.only <- T
		gp_manhattan_pNode <- oPierManhattan(pNode, color=c("darkred", "#2071b8"), point.size=0.2, top=highlight.top, top.label.type="text", top.label.size=2.5, top.label.col="black", top.label.query=top.label.query, label.query.only=label.query.only, y.lab="Core gene scores\n(scored 0-100)", signature=F, placeholder=placeholder)
		output.file.manhattan.pdf <- gsub(".txt$", "_manhattan.pdf", output.file, perl=T)
		ggsave(output.file.manhattan.pdf, gp_manhattan_pNode, device=cairo_pdf, width=9, height=3)
		output.file.manhattan.png <- gsub(".txt$", "_manhattan.png", output.file, perl=T)
		ggsave(output.file.manhattan.png, gp_manhattan_pNode, type="cairo", width=9, height=3)
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/PiER-site/pier_app/public
		## but outputs at public/tmp/eV2CG.SNPs.STRING_high.72959383_priority.xlsx
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$runtime <- str_c(runtime," seconds")
		ls_rmd$data_input <- data %>% set_names(c("SNP","PVALUE"))
		ls_rmd$highlight_top <- highlight.top
		
		ls_rmd$xlsx_CG <- gsub("public/", "", output.file.CG, perl=T)
		ls_rmd$xlsx_CG_evidence <- gsub("public/", "", output.file.CG_evidence, perl=T)
		ls_rmd$pdf_manhattan <- gsub("public/", "", output.file.manhattan.pdf, perl=T)
		ls_rmd$png_manhattan <- gsub("public/", "", output.file.manhattan.png, perl=T)
		
		output_dir <- gsub("eV2CG.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_eV2CG.Rmd", bookdown::html_document2(number_sections=F,theme="readable", hightlight="default"), output_dir=output_dir)
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

R_pipeline(input.file="public/tmp/data.SNPs.STRING_highest.53083304.txt", output.file="public/tmp/eV2CG.SNPs.STRING_highest.53083304.txt", population="EUR", crosslink="proximity_20000", network="STRING_highest", significance.threshold="5e-8", restart="0.7", highlight.top="30", placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('
--- eV2CG: ',runTime,' secs ---
'), appendLF=TRUE)
