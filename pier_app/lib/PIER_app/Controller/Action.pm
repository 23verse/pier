package PIER_app::Controller::Action;
use PIER_app::Controller::Utils;
use Mojo::Base 'Mojolicious::Controller';
#use JSON;
use JSON::Parse;
use LWP::Simple;

# Render template "index.html.ep"
# Render template "demo.html.ep"
sub index {
	my $c = shift;
  	$c->render();
}



# Render template "cTCrosstalk.html.ep"
sub PiER_cTCrosstalk_default {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $snplist = $c->param('snplist');
  	my $population = $c->param('pop') || 'NA'; # by default: NA
	my $distance_max = $c->param('distance');
	my $include_QTL = $c->param('QTL') || 'NA';
	my $include_RGB = $c->param('RGB') || 'NA';
  	my $network = $c->param('network') || 'STRING_highest'; # by default: STRING_highest
  	
	my $significance_threshold = $c->param('significance_threshold');
  	
	my $restart = $c->param('restart');
	my $highlight_top = $c->param('highlight_top');
	my $pathway_top = $c->param('pathway_top');
	my $subnet_size = $c->param('subnet_size');
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($snplist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $network.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.SNPs.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'cTCrosstalk.SNPs.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'cTCrosstalk.SNPs.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n/, $snplist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
		}
		PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}elsif(-e '/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS'){
			$placeholder="/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS";
		}elsif(-e '/var/www/bigdata_dev'){
			$placeholder="/var/www/bigdata_dev";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", population="", distance.max="", include.QTL="", include.RGB="", network="", significance.threshold="", restart="", highlight.top="", pathway.top="", subnet.size="", placeholder="", ...){
	
	#############################
	## an example for coding test
	if(FALSE){
		
		if(FALSE){
			## copy Fang package R functions to galahad and source
			#rsync -avzhe ssh --progress --delete /Users/hfang/Sites/XGR/Fang galahad.well.ox.ac.uk:./
			#cd ~
			#/home/hfang/R-3.6.2/bin/R
			vec <- list.files(path="/home/hfang/Fang/R", pattern=".r", full.names=T)
			ls_tmp <- lapply(vec, function(x) source(x))
			
			## copy Fang package R functions to huawei and source
			#rsync -avzhe ssh --progress --delete /Users/hfang/Sites/XGR/Fang root@94.74.100.58:./
			#cd ~
			#/usr/local/bin/R
			vec <- list.files(path="/root/Fang/R", pattern=".r", full.names=T)
			ls_tmp <- lapply(vec, function(x) source(x))
			
		}else{
			#cd /Users/hfang/Sites/XGR/PiER-site/app/examples
			#R
			vec <- list.files(path="/Users/hfang/Sites/XGR/Fang/R", pattern=".r", full.names=T)
			ls_tmp <- lapply(vec, function(x) source(x))	
		}
		
		#library(Fang)
		#library(patchwork)
		library(tidyverse)
		library(igraph)
		library(GenomicRanges)
		guid <- NULL
		
		# mac
		placeholder <- "http://galahad.well.ox.ac.uk/bigdata_dev"
		placeholder <- "~/Sites/SVN/github/bigdata_dev"
		tb <- read_delim("~/Sites/XGR/PiER-site/app/examples/GWAS.cross.SNP.txt", delim="\t")
		tb <- read_delim("./GWAS.cross.SNP.txt", delim="\t")
		## data_psych
		data_psych <- tb %>% filter(pubmedid==31835028) %>% distinct(snp_id_current,pvalue) %>% as.data.frame()
		data_psych %>% filter(pvalue<5e-8) %>% count(snp_id_current)
		## data_inflam
		data_inflam <- tb %>% filter(pubmedid==26974007) %>% distinct(snp_id_current,pvalue) %>% as.data.frame()
		data_inflam %>% filter(pvalue<5e-8) %>% count(snp_id_current)
		data_inflam %>% openxlsx::write.xlsx("~/Sites/XGR/PiER-site/app/examples/data_inflam.xlsx")
		data <- openxlsx::read.xlsx("~/Sites/XGR/PiER-site/app/examples/data_inflam.xlsx")
		
		# galahad
		placeholder <- "/var/www/bigdata_dev"
		placeholder <- "/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS"
		tb <- read_delim("~/pier_app/public/app/examples/GWAS.cross.SNP.txt", delim="\t")
		data <- openxlsx::read.xlsx("~/pier_app/public/app/examples/data_inflam.xlsx")
		
		# huawei
		placeholder <- "/var/www/bigdata_dev"
		tb <- read_delim("~/pier_app/public/app/examples/GWAS.cross.SNP.txt", delim="\t")
		data <- openxlsx::read.xlsx("~/pier_app/public/app/examples/data_inflam.xlsx")
		
		#dbSNP_Common <- oRDS("dbSNP_Common", placeholder="http://galahad.well.ox.ac.uk/bigdata_dev")
		#dbSNP_GWAS <- Pi::xRDataLoader("dbSNP_GWAS", RData.location="http://galahad.well.ox.ac.uk/bigdata")
		
		############
		# Parameters
		############
		LD.customised <- oRDS("GWAS_LD.EUR", placeholder=placeholder, guid=guid) %>% as.data.frame()
		significance.threshold <- 5e-8
		distance.max <- 20000
		decay.kernel <- "constant"
		restart <- 0.7
		subnet.size <- 30
		#GR.SNP <- oRDS("dbSNP_GWAS", placeholder=placeholder, guid=guid)
		GR.SNP <- oRDS("dbSNP_Common", placeholder=placeholder, guid=guid)
		GR.Gene <- oRDS("UCSC_knownGene", placeholder=placeholder, guid=guid)
		include.TAD <- NULL
		include.RGB <- c("PCHiC_PMID27863249_Combined", "PCHiC_PMID27863249_Activated_total_CD4_T_cells","PCHiC_PMID27863249_Endothelial_precursors","PCHiC_PMID27863249_Erythroblasts","PCHiC_PMID27863249_Fetal_thymus","PCHiC_PMID27863249_Macrophages_M0","PCHiC_PMID27863249_Macrophages_M1","PCHiC_PMID27863249_Macrophages_M2","PCHiC_PMID27863249_Megakaryocytes","PCHiC_PMID27863249_Monocytes","PCHiC_PMID27863249_Naive_B_cells","PCHiC_PMID27863249_Naive_CD4_T_cells","PCHiC_PMID27863249_Naive_CD8_T_cells","PCHiC_PMID27863249_Neutrophils","PCHiC_PMID27863249_Nonactivated_total_CD4_T_cells","PCHiC_PMID27863249_Total_B_cells","PCHiC_PMID27863249_Total_CD4_T_cells","PCHiC_PMID27863249_Total_CD8_T_cells")[10]
		include.QTL <- c("eQTL_eQTLGen", "pQTL_Plasma")[2]
		
		network.customised <- oDefineNet(network="STRING_high", STRING.only=c("experimental_score","database_score"), placeholder=placeholder) %>% oNetInduce(nodes_query=names(GR.Gene), knn=0, largest.comp=F)
		
		message(sprintf("Customised the network with %d nodes and %d edges (%s) ...", gorder(network.customised), gsize(network.customised), as.character(Sys.time())), appendLF=TRUE)
		
		###################
		# genomic predictors
		###################
		ls_pNode_genomic <- oPierSNPsAdv(data, score.cap=100, LD.customised=LD.customised, significance.threshold=significance.threshold, distance.max=distance.max, decay.kernel=decay.kernel, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.TAD=include.TAD, include.QTL=include.QTL, include.RGB=include.RGB, network.customised=network.customised, placeholder=placeholder, guid=guid, verbose.details=T)
		
		###################
		# prioritisation
		###################
		dTarget <- oPierMatrix(ls_pNode_genomic, displayBy="pvalue", aggregateBy="fishers", placeholder=placeholder, guid=guid)
		dTarget$priority %>% as_tibble() %>% transmute(Gene=name, Rank=rank, Rating=signif(rating,digits=4), Type=ifelse(seed=="Y","Core","Peripheral"), Proximity=nGene, QTL=eGene, PCHiC=cGene, Description=description) %>% openxlsx::write.xlsx(str_c("dTarget_priority.xlsx"))
		
		## Evidence
		df_Gene2SNP <- oPierEvidence(ls_pNode_genomic) %>% as_tibble() %>% arrange(Flag,-Pval) %>% transmute(Gene,SNP,SNP_type=ifelse(Flag=="Lead","Input","LD"),Evidence=Context) %>% distinct() %>% mutate(Evidence=str_replace_all(Evidence,"nGene_","Proximity_")) %>% mutate(Evidence=str_replace_all(Evidence,"_constant","bp")) %>% mutate(Evidence=str_replace_all(Evidence,"eGene_","QTL_")) %>% mutate(Evidence=str_replace_all(Evidence,"cGene_",""))
		#df_Gene2SNP %>% count(Evidence)
		df_evidence <- dTarget$priority %>% as_tibble() %>% filter(seed=="Y") %>% transmute(Gene=name, Rating=signif(rating,digits=4), Proximity=nGene, QTL=eGene, PCHiC=cGene) %>% left_join(df_Gene2SNP, by="Gene") %>% arrange(-Rating, SNP_type, SNP, Evidence)
		#df_evidence %>% filter(Gene=="IL23R") %>% print(n=Inf)
		#df_evidence %>% filter(Gene=="IL6R") %>% print(n=Inf)
		df_evidence %>% openxlsx::write.xlsx(str_c("dTarget_evidence.xlsx"))
		
		###################
		# Manhattan plot
		###################
		top.label.query <- dTarget$priority %>% top_frac(1,rating) %>% pull(name)
		label.query.only <- T
		gp_manhattan_dTarget <- oPierManhattan(dTarget, color=c("darkred", "#2071b8"), point.size=0.2, top=30, top.label.type="text", top.label.size=2.5, top.label.col="black", top.label.query=top.label.query, label.query.only=label.query.only, y.lab="Priority rating (scored 0-5)", signature=F, placeholder=placeholder, guid=guid)
		ggsave("dTarget_manhattan.pdf", gp_manhattan_dTarget, width=9, height=3)
		ggsave("dTarget_manhattan.png", gp_manhattan_dTarget, type="cairo", width=9, height=3)
		
		###################
		# KEGG analysis
		###################
		vec <- dTarget$priority %>% top_frac(0.01,rating) %>% pull(name)
		background <- dTarget$priority %>% pull(name)
		set <- oRDS("org.Hs.egKEGG", placeholder=placeholder, guid=guid)
		eset <- oSEA(vec, set, background=background, size.range=c(15,1500), test="fisher", min.overlap=10)
		df_eTerm <- eset %>% oSEAextract() %>% filter(namespace %in% c("Environmental Process","Organismal System","Cellular Process","Human Disease","Genetic Process","Metabolism")[c(1,2,3,5,6)])
		#df_eTerm <- eset %>% oSEAextract()
		gp_kegg <- df_eTerm %>% filter(distance==3) %>% mutate(onto=namespace) %>% oSEAdotplot(label.top=5, size.title="Number of member genes")
		ggsave("dTarget_KEGG.pdf", gp_kegg, width=5, height=4)
		ggsave("dTarget_KEGG.png", gp_kegg, type="cairo", width=5, height=4)	
		df_eTerm %>% filter(distance==3) %>% as_tibble() %>% openxlsx::write.xlsx(str_c("dTarget_KEGG.xlsx"))

		###################
		# Pathway Crosstalk
		###################
		ig.KEGG.category <- oRDS("ig.KEGG.category", placeholder=placeholder, guid=guid)
		ls_ig <- ig.KEGG.category %>% filter(category %in% c("Organismal System","Environmental Process")) %>% deframe()
		ig <- oCombineNet(ls_ig)
		ig2 <- oNetInduce(ig, nodes_query=V(ig)$name, largest.comp=T) %>% as.undirected()
		subg <- oPierSubnet(dTarget, network.customised=ig2, subnet.size=30)
		subg <- subg %>% oLayout(c("layout_with_kk","graphlayouts.layout_with_stress")[2])
		gp_rating_evidence <- oVisEvidenceAdv(dTarget, nodes=V(subg)$name, g=subg, node.info="smart", node.label.size=2.5, node.label.color="black", node.label.alpha=0.95, node.label.force=0.2, node.color.title="Priority\nrating", colormap="spectral", node.color.alpha=0.95, node.size.range=5, edge.color="steelblue4", edge.color.alpha=0.5, edge.size=0.3, edge.curve=0.05, pie.color.alpha=0, pie.radius=0)
		gp_rating <- gp_rating_evidence + guides(fill="none")
		gp_rating
		ggsave("dTarget_crosstalk.pdf", gp_rating, device=cairo_pdf, width=6, height=6)
		ggsave("dTarget_crosstalk.png", gp_rating, type="cairo", width=6, height=6)
		df_subg <- dTarget$priority %>% as_tibble() %>% transmute(Gene=name, Rank=rank, Rating=signif(rating,digits=4), Proximity=nGene, QTL=eGene, PCHiC=cGene, Description=description) %>% semi_join(subg %>% oIG2TB("node"), by=c("Gene"="name"))
		df_subg %>% openxlsx::write.xlsx(str_c("dTarget_crosstalk.xlsx"))
		
        #igraph::degree(subg) %>% sort()
        #subg.sig <- oPierSubnet(dTarget, network.customised=ig2, subnet.size=30, test.permutation=T, num.permutation=10, respect=c("none","degree")[2], aggregateBy="fishers")
		#subg.sig$combinedP = 6.61e-130
		
		df_evidence_crosstalk <- df_evidence %>% semi_join(df_subg, by="Gene")
		df_evidence_crosstalk %>% openxlsx::write.xlsx(str_c("dTarget_crosstalk_evidence.xlsx"))
		
		###################
		# Drug repurposing
		###################
		ChEMBL_v30_DTT <- oRDS("ChEMBL_v30_DTT", placeholder=placeholder)
		#ChEMBL_v30_DTT <- ChEMBL_v30_DTT %>% mutate(efo_term=str_replace_all(efo_term,"\\{http.*",""))
		ChEMBL_v30_DTT <- ChEMBL_v30_DTT %>% filter(!(efo_term %in% c("cancer","neoplasm","immune system disease"))) %>% filter(!str_detect(efo_term,"neoplasm")) %>% filter(!str_detect(efo_term,"cancer")) %>%  mutate(pref_name_drug=str_to_lower(pref_name_drug)) %>% mutate(pref_name_drug=str_c(pref_name_drug," [",mechanism_of_action,"]")) %>% mutate(efo_term=str_to_title(efo_term))
		
		dtt <- ChEMBL_v30_DTT %>% mutate(target_number=1) %>% as.data.frame()
		vec_data <- dtt %>% select(Symbol,phase) %>% distinct() %>% arrange(-phase,Symbol) %>% semi_join(tibble(Symbol=V(subg)$name), by="Symbol") %>% pull(Symbol)
		DR <- oRepurpose(vec_data, phase.min=4, target.max=5, reorder="none", DTT=dtt, restricted=NULL, colormap="spectral.top", zlim=c(1,4), na.color="transparent", label.size=2.5, label.color="white", x.rotate=70, size=4, legend.title="Approved", x.text.size=7, y.text.size=7)
		gp_DR <- DR$gp + theme(legend.position="none") + scale_y_discrete(position="right")
		n_targets <- DR$df %>% count(Target) %>% nrow()
		n_diseases <- DR$df %>% count(Disease) %>% nrow()
		ggsave("Crosstalk.approved.pdf", gp_DR, device=cairo_pdf, width=n_diseases*0.16, height=2.5+n_targets*0.16)
		ggsave("Crosstalk.approved.png", gp_DR, type="cairo", width=n_diseases*0.16, height=2.5+n_targets*0.16)
		DR$df %>% openxlsx::write.xlsx("Crosstalk.approved.xlsx")
		#DR$index %>% openxlsx::write.xlsx("Crosstalk.approved_index.xlsx")
		
		if(0){
			library(gridExtra)
			tt <- ttheme_default(base_size=6, padding=unit(c(2,2),"mm"))
			tt_left <- ttheme_default(base_size=6, padding=unit(c(2,2),"mm"), core=list(fg_params=list(hjust=0,x=0.01)), colhead=list(fg_params=list(hjust=0, x=0.01)))
			tt_right <- ttheme_default(base_size=6, padding=unit(c(2,2),"mm"), core=list(fg_params=list(hjust=1,x=0.98)), colhead=list(fg_params=list(hjust=1, x=0.98)))
			# one-column table
			gp_index <- DR$index %>% mutate(Drug=str_replace_all(Drug,",","\n")) %>% gridExtra::tableGrob(rows=NULL,cols=c("Index","Approved drugs [mechanisms of action]"), theme=tt_right)
			grid::grid.draw(gp_index)
			
			# two-column tables
			tmp <- DR$index %>% separate_rows(Drug, sep=", ") 
			n <- tmp %>% dplyr::slice(ceiling(tmp %>% nrow()/2)) %>% pull(Drug_index)
			df <- DR$index %>% mutate(Drug=str_replace_all(Drug,", ","\n"))
			m <- df %>% nrow()
			gt1 <- df %>% dplyr::slice(1:n) %>% gridExtra::tableGrob(rows=NULL,cols=c("Index","Approved drugs [mechanisms of action]"), theme=tt_right)
			gt2 <- df %>% dplyr::slice((n+1):m) %>% gridExtra::tableGrob(rows=NULL,cols=c("Index","Approved drugs [mechanisms of action]"), theme=tt_right)
			ls_gt <- list(gt1, gt2)
			gp_table <- ls_gt %>% gridExtra::marrangeGrob(ncol=2, nrow=1, as.table=FALSE)
			gp_table
		}
		
	}
	#############################
	
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
	pathway.top <- as.numeric(pathway.top)
	subnet.size <- as.numeric(subnet.size)
	
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
	
	GR.SNP <- oRDS("dbSNP_GWAS", placeholder=placeholder)
	GR.Gene <- oRDS("UCSC_knownGene", placeholder=placeholder)
	network.customised <- oDefineNet(network=network, STRING.only=c("experimental_score","database_score"), placeholder=placeholder) %>% oNetInduce(nodes_query=names(GR.Gene), knn=0, largest.comp=F)
	
	message(sprintf("Preparing predictors (%s) ...", as.character(Sys.time())), appendLF=TRUE)
	
	ls_pNode_genomic <- oPierSNPsAdv(data, score.cap=100, LD.customised=LD.customised, significance.threshold=significance.threshold, distance.max=distance.max, decay.kernel="constant", GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.TAD=NULL, include.QTL=include.QTL, include.RGB=include.RGB, network.customised=network.customised, placeholder=placeholder, verbose.details=F)
	
	message(sprintf("Doing prioritisation (%s) ...", as.character(Sys.time())), appendLF=TRUE)
	
	dTarget <- oPierMatrix(ls_pNode_genomic, displayBy="pvalue", aggregateBy="fishers", placeholder=placeholder)

	if(class(dTarget)=="dTarget"){
		# *_priority.txt
		dTarget$priority %>% as_tibble() %>% write_delim(output.file, delim="\t")
		
		# *_priority.xlsx
		output.file.priority <- gsub(".txt$", "_priority.xlsx", output.file, perl=T)
		dTarget$priority %>% as_tibble() %>% openxlsx::write.xlsx(output.file.priority)
		
		# Manhattan plot
		message(sprintf("Drawing manhattan plot (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		top.label.query <- dTarget$priority %>% top_frac(1,rating) %>% pull(name)
		label.query.only <- T
		gp_manhattan_dTarget <- oPierManhattan(dTarget, color=c("darkred", "darkorange"), point.size=0.2, top=30, top.label.type="text", top.label.size=2, top.label.col="black", top.label.query=top.label.query, label.query.only=label.query.only, y.lab="Pi rating", signature=F, placeholder=placeholder)
		output.file.manhattan.pdf <- gsub(".txt$", "_manhattan.pdf", output.file, perl=T)
		ggsave(output.file.manhattan.pdf, gp_manhattan_dTarget, width=9, height=3)
		output.file.manhattan.png <- gsub(".txt$", "_manhattan.png", output.file, perl=T)
		#ggsave(output.file.manhattan.png, gp_manhattan_dTarget, device=png(), width=9, height=3)
		
		################
		# KEGG analysis
		################
		message(sprintf("Performing KEGG analysis based on %e genes (%s) ...", pathway.top, as.character(Sys.time())), appendLF=TRUE)
		vec <- dTarget$priority %>% top_frac(pathway.top, rating) %>% pull(name)
		set <- oRDS("org.Hs.egKEGG", placeholder=placeholder)
		eset <- oSEA(vec, set, size.range=c(15,1500), test="fisher", min.overlap=10)
		df_eTerm <- eset %>% oSEAextract() %>% filter(namespace %in% c("Environmental Process","Organismal System","Cellular Process","Human Disease","Genetic Process","Metabolism")[c(1,2,3,5,6)])
		
		# *_kegg.xlsx
		output.file.kegg <- gsub(".txt$", "_kegg.xlsx", output.file, perl=T)
		df_eTerm %>% filter(distance==3) %>% as_tibble() %>% openxlsx::write.xlsx(output.file.kegg)
		
		# *_kegg.pdf *_kegg.png
		gp_kegg <- df_eTerm %>% filter(distance==3) %>% mutate(onto=namespace) %>% oSEAdotplot(label.top=10, size.title="Number of member genes")
		output.file.kegg.pdf <- gsub(".txt$", "_kegg.pdf", output.file, perl=T)
		ggsave(output.file.kegg.pdf, gp_kegg, width=5, height=4)
		output.file.kegg.png <- gsub(".txt$", "_kegg.png", output.file, perl=T)
		#ggsave(output.file.kegg.png, gp_kegg, device=png(), width=5, height=4)
		
		
		###################
		# Pathway Crosstalk
		###################
		message(sprintf("Performing pathway crosstalk analysis restricted to %d network genes (%s) ...", subnet.size, as.character(Sys.time())), appendLF=TRUE)
		ig.KEGG.category <- oRDS("ig.KEGG.category", placeholder=placeholder)
		ls_ig <- ig.KEGG.category %>% filter(category %in% c("Organismal System","Environmental Process")) %>% deframe()
		ig <- oCombineNet(ls_ig)
		ig2 <- oNetInduce(ig, nodes_query=V(ig)$name, largest.comp=T) %>% as.undirected()
		subg <- oPierSubnet(dTarget, network.customised=ig2, subnet.size=subnet.size)
		subg <- subg %>% oLayout(c("layout_with_kk","graphlayouts.layout_with_stress")[2])
		gp_rating_evidence <- oVisEvidenceAdv(dTarget, nodes=V(subg)$name, g=subg, node.info="smart", node.label.size=2.5, node.label.color="steelblue4", node.label.force=0.4, colormap="spectral.top", node.color.alpha=0.7, node.size.range=4, edge.color="steelblue", edge.color.alpha=0.2, edge.curve=0.05, pie.color.alpha=0, pie.radius=0)
		gp_rating <- gp_rating_evidence + guides(fill="none")

		# *_crosstalk.xlsx
		output.file.crosstalk <- gsub(".txt$", "_crosstalk.xlsx", output.file, perl=T)
		df_subg <- dTarget$priority %>% as_tibble() %>% select(name,rank,rating,description) %>% semi_join(subg %>% oIG2TB("node"), by="name")
		df_subg %>% openxlsx::write.xlsx(output.file.crosstalk)
		
		# *_crosstalk.pdf *_crosstalk.png
		output.file.crosstalk.pdf <- gsub(".txt$", "_crosstalk.pdf", output.file, perl=T)
		ggsave(output.file.crosstalk.pdf, gp_rating, width=6, height=6)
		output.file.crosstalk.png <- gsub(".txt$", "_crosstalk.png", output.file, perl=T)
		#ggsave(output.file.crosstalk.png, gp_rating, device=png(), width=6, height=6)
		
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/PiER-site/pier_app/public
		## but outputs at public/tmp/cTCrosstalk.SNPs.STRING_high.72959383_priority.xlsx
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		ls_rmd <- list()
		ls_rmd$data_input <- data %>% set_names(c("SNP","PVALUE"))
		ls_rmd$highlight_top <- highlight.top
		
		ls_rmd$xlsx_priority <- gsub("public/", "", output.file.priority, perl=T)
		ls_rmd$pdf_manhattan <- gsub("public/", "", output.file.manhattan.pdf, perl=T)
		ls_rmd$png_manhattan <- gsub("public/", "", output.file.manhattan.png, perl=T)
		
		ls_rmd$xlsx_kegg <- gsub("public/", "", output.file.kegg, perl=T)
		ls_rmd$pdf_kegg <- gsub("public/", "", output.file.kegg.pdf, perl=T)
		ls_rmd$png_kegg <- gsub("public/", "", output.file.kegg.png, perl=T)
		
		ls_rmd$xlsx_crosstalk <- gsub("public/", "", output.file.crosstalk, perl=T)
		ls_rmd$pdf_crosstalk <- gsub("public/", "", output.file.crosstalk.pdf, perl=T)
		ls_rmd$png_crosstalk <- gsub("public/", "", output.file.crosstalk.png, perl=T)
		
		output_dir <- gsub("cTCrosstalk.*", "", output.file, perl=T)
		
		Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		rmarkdown::render("public/RMD_cTCrosstalk.Rmd", bookdown::html_document2(number_sections=F,theme="readable", hightlight="default"), output_dir=output_dir)
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
startT <- Sys.time()

library(tidyverse)
library(igraph)
vec <- list.files(path='/home/hfang/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", population=\"$population\", distance.max=\"$distance_max\", include.QTL=\"$include_QTL\", include.RGB=\"$include_RGB\", network=\"$network\", significance.threshold=\"$significance_threshold\", restart=\"$restart\", highlight.top=\"$highlight_top\", pathway.top=\"$pathway_top\", subnet.size=\"$subnet_size\", placeholder=\"$placeholder\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('cTCrosstalk_default: ',runTime,' secs\n'), appendLF=TRUE)
";

# for calling R function
PIER_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
if(-e $rscript_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			#############
			### for _priority.xlsx
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/_priority.xlsx/g;
			if(-e $tmp_file){
				$ajax_priority_xlsx_file=$tmp_file;
				$ajax_priority_xlsx_file=~s/^public//g;
				print STDERR "_priority.xlsx locates at $ajax_priority_xlsx_file\n";
			}
			
			#############
			### for _manhattan.pdf
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/_manhattan.pdf/g;
			if(-e $tmp_file){
				$ajax_manhattan_pdf_file=$tmp_file;
				$ajax_manhattan_pdf_file=~s/^public//g;
				print STDERR "_manhattan.pdf locates at $ajax_manhattan_pdf_file\n";
			}
			
			#############
			### for RMD_cTCrosstalk.html
			$tmp_file=$output_filename;
			$tmp_file=~s/cTCrosstalk.*//g;
			$tmp_file=$tmp_file."RMD_cTCrosstalk.html";
			if(-e $tmp_file){
				$ajax_rmd_html_file=$tmp_file;
				$ajax_rmd_html_file=~s/^public//g;
				print STDERR "RMD_cTCrosstalk locates at $ajax_rmd_html_file\n";
			}
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);
	
	# stash $xlsxfile;
	my $xlsxfile;
	$xlsxfile->{priority}=$ajax_priority_xlsx_file;
	$c->stash(xlsxfile => $xlsxfile);

	# stash $pdffile;
	my $pdffile;
	$pdffile->{manhattan}=$ajax_manhattan_pdf_file;
	$c->stash(pdffile => $pdffile);
	
  	$c->render();

}



# Render template "cTCrosstalk.html.ep"
sub PiER_cTCrosstalk {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $snplist = $c->param('snplist');
  	my $population = $c->param('pop') || 'NA'; # by default: NA
	my $distance_max = $c->param('distance');
	my $include_QTL = $c->param('QTL') || 'NA';
	my $include_RGB = $c->param('RGB') || 'NA';
  	my $network = $c->param('network') || 'STRING_highest'; # by default: STRING_highest
  	
	my $significance_threshold = $c->param('significance_threshold');
  	
	my $restart = $c->param('restart');
	my $highlight_top = $c->param('highlight_top');
	my $pathway_top = $c->param('pathway_top');
	my $subnet_size = $c->param('subnet_size');
  	
	my $crosstalk_sig = $c->param('crosstalk_sig');
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($snplist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $network.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.SNPs.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'cTCrosstalk.SNPs.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'cTCrosstalk.SNPs.'.$rand_file.'.r';
	
		my $my_input;
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $snplist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			# mac
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}elsif(-e '/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS'){
			# galahad
			$placeholder="/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS";
		}elsif(-e '/var/www/bigdata_dev'){
			# huawei
			$placeholder="/var/www/bigdata_dev";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", population="", distance.max="", include.QTL="", include.RGB="", network="", significance.threshold="", restart="", highlight.top="", pathway.top="", subnet.size="", crosstalk.sig="", placeholder="", ...){
	
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
	subnet.size <- as.numeric(subnet.size)
	highlight.top <- as.numeric(highlight.top)
	pathway.top <- as.numeric(pathway.top)
	
	#population <- "EUR"
	#highlight.top <- 30
	#pathway.top <- 0.01
	#network <- "STRING_high"
	
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
		
		###################
		# KEGG analysis
		###################
		vec <- dTarget$priority %>% top_frac(pathway.top,rating) %>% pull(name)
		background <- dTarget$priority %>% pull(name)
		set <- oRDS("org.Hs.egKEGG", placeholder=placeholder)
		eset <- oSEA(vec, set, background=background, size.range=c(15,1500), test="fisher", min.overlap=10)
		df_eTerm <- eset %>% oSEAextract() %>% filter(namespace %in% c("Environmental Process","Organismal System","Cellular Process","Human Disease","Genetic Process","Metabolism")[c(1,2,3,5,6)])
		gp_kegg <- df_eTerm %>% filter(distance==3) %>% mutate(onto=namespace) %>% oSEAdotplot(label.top=5, size.title="Number of member genes")
		output.file.kegg.pdf <- gsub(".txt$", "_kegg.pdf", output.file, perl=T)
		ggsave(output.file.kegg.pdf, gp_kegg, device=cairo_pdf, width=5, height=4)
		output.file.kegg.png <- gsub(".txt$", "_kegg.png", output.file, perl=T)
		ggsave(output.file.kegg.png, gp_kegg, type="cairo", width=5, height=4)
		# *_kegg.xlsx
		output.file.kegg <- gsub(".txt$", "_kegg.xlsx", output.file, perl=T)
		df_eTerm %>% filter(distance==3) %>% as_tibble() %>% openxlsx::write.xlsx(output.file.kegg)
		
		###################
		# Pathway Crosstalk
		###################
		message(sprintf("Performing pathway crosstalk analysis restricted to %d network genes (%s) ...", subnet.size, as.character(Sys.time())), appendLF=TRUE)
		ig.KEGG.category <- oRDS("ig.KEGG.category", placeholder=placeholder)
		ls_ig <- ig.KEGG.category %>% filter(category %in% c("Organismal System","Environmental Process")) %>% deframe()
		ig <- oCombineNet(ls_ig)
		ig2 <- oNetInduce(ig, nodes_query=V(ig)$name, largest.comp=T) %>% as.undirected()
		subg <- oPierSubnet(dTarget, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder)
		
		subg <- subg %>% oLayout(c("layout_with_kk","graphlayouts.layout_with_stress")[2])
		gp_rating_evidence <- oVisEvidenceAdv(dTarget, nodes=V(subg)$name, g=subg, node.info="smart", node.label.size=2.5, node.label.color="black", node.label.alpha=0.95, node.label.force=0.4, node.color.title="Priority\nrating", colormap="spectral", node.color.alpha=0.8, node.size.range=5, edge.color="steelblue4", edge.color.alpha=0.5, edge.size=0.3, edge.curve=0.05, pie.color.alpha=0, pie.radius=0)
		gp_rating <- gp_rating_evidence + guides(fill="none")
		
		# *_crosstalk.xlsx
		output.file.crosstalk <- gsub(".txt$", "_crosstalk.xlsx", output.file, perl=T)
		df_subg <- dTarget$priority %>% as_tibble() %>% transmute(Gene=name, Rank=rank, Rating=signif(rating,digits=4), Proximity=nGene, QTL=eGene, PCHiC=cGene, Description=description) %>% semi_join(subg %>% oIG2TB("node"), by=c("Gene"="name"))
		df_subg %>% openxlsx::write.xlsx(output.file.crosstalk)
		
		# *_crosstalk_evidence.xlsx
		output.file.crosstalk_evidence <- gsub(".txt$", "_crosstalk_evidence.xlsx", output.file, perl=T)
		df_evidence_crosstalk <- df_evidence %>% semi_join(df_subg, by="Gene")
		df_evidence_crosstalk %>% openxlsx::write.xlsx(output.file.crosstalk_evidence)
		
		# *_crosstalk.pdf *_crosstalk.png
		output.file.crosstalk.pdf <- gsub(".txt$", "_crosstalk.pdf", output.file, perl=T)
		ggsave(output.file.crosstalk.pdf, gp_rating, device=cairo_pdf, width=6, height=6)
		output.file.crosstalk.png <- gsub(".txt$", "_crosstalk.png", output.file, perl=T)
		ggsave(output.file.crosstalk.png, gp_rating, type="cairo", width=6, height=6)
		
		combinedP <- 1
		if(crosstalk.sig=="yes"){
			subg.sig <- oPierSubnet(dTarget, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder, test.permutation=T, num.permutation=10, respect=c("none","degree")[2], aggregateBy="fishers")
			combinedP <- signif(subg.sig$combinedP, digits=2)
		}
		
		###################
		# Drug repurposing
		###################
		message(sprintf("Performing drug repurposing analysis based on %d crosstalk genes (%s) ...", vcount(subg), as.character(Sys.time())), appendLF=TRUE)
		ChEMBL_v30_DTT <- oRDS("ChEMBL_v30_DTT", placeholder=placeholder)
		#ChEMBL_v30_DTT <- ChEMBL_v30_DTT %>% mutate(efo_term=str_replace_all(efo_term,"\\{http.*",""))
		ChEMBL_v30_DTT <- ChEMBL_v30_DTT %>% filter(!(efo_term %in% c("cancer","neoplasm","immune system disease"))) %>% filter(!str_detect(efo_term,"neoplasm")) %>% filter(!str_detect(efo_term,"cancer")) %>%  mutate(pref_name_drug=str_to_lower(pref_name_drug)) %>% mutate(pref_name_drug=str_c(pref_name_drug," [",mechanism_of_action,"]")) %>% mutate(efo_term=str_to_title(efo_term))
		## prepare ddt and vec_data for drug repurposing (approved)
		dtt <- ChEMBL_v30_DTT %>% mutate(target_number=1) %>% as.data.frame()
		vec_data <- dtt %>% select(Symbol,phase) %>% distinct() %>% arrange(-phase,Symbol) %>% semi_join(tibble(Symbol=V(subg)$name), by="Symbol") %>% pull(Symbol)
		DR <- oRepurpose(vec_data, phase.min=4, target.max=5, reorder="none", DTT=dtt, restricted=NULL, colormap="spectral.top", zlim=c(1,4), na.color="transparent", label.size=2.5, label.color="white", x.rotate=70, size=4, legend.title="Approved", x.text.size=7, y.text.size=7)
		gp_DR <- DR$gp + theme(legend.position="none") + scale_y_discrete(position="right")
		n_targets <- DR$df %>% count(Target) %>% nrow()
		n_diseases <- DR$df %>% count(Disease) %>% nrow()
		
		# *_repurposing.xlsx
		output.file.repurposing <- gsub(".txt$", "_repurposing.xlsx", output.file, perl=T)
		DR$df %>% openxlsx::write.xlsx(output.file.repurposing)
		
		# *_repurposing.pdf *_repurposing.png
		output.file.repurposing.pdf <- gsub(".txt$", "_repurposing.pdf", output.file, perl=T)
		ggsave(output.file.repurposing.pdf, gp_DR, device=cairo_pdf, width=n_diseases*0.16, height=2.5+n_targets*0.16)
		output.file.repurposing.png <- gsub(".txt$", "_repurposing.png", output.file, perl=T)
		ggsave(output.file.repurposing.png, gp_DR, type="cairo", width=n_diseases*0.16, height=2.5+n_targets*0.16)
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/PiER-site/pier_app/public
		## but outputs at public/tmp/cTCrosstalk.SNPs.STRING_high.72959383_priority.xlsx
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$runtime <- str_c(runtime," seconds")
		ls_rmd$data_input <- data %>% set_names(c("SNP","PVALUE"))
		ls_rmd$highlight_top <- highlight.top
		ls_rmd$pathway_top <- pathway.top
		ls_rmd$combinedP <- combinedP
		
		ls_rmd$xlsx_priority <- gsub("public/", "", output.file.priority, perl=T)
		ls_rmd$xlsx_evidence <- gsub("public/", "", output.file.evidence, perl=T)
		ls_rmd$pdf_manhattan <- gsub("public/", "", output.file.manhattan.pdf, perl=T)
		ls_rmd$png_manhattan <- gsub("public/", "", output.file.manhattan.png, perl=T)
		
		ls_rmd$xlsx_kegg <- gsub("public/", "", output.file.kegg, perl=T)
		ls_rmd$pdf_kegg <- gsub("public/", "", output.file.kegg.pdf, perl=T)
		ls_rmd$png_kegg <- gsub("public/", "", output.file.kegg.png, perl=T)
		
		ls_rmd$xlsx_crosstalk <- gsub("public/", "", output.file.crosstalk, perl=T)
		ls_rmd$xlsx_crosstalk_evidence <- gsub("public/", "", output.file.crosstalk_evidence, perl=T)
		ls_rmd$pdf_crosstalk <- gsub("public/", "", output.file.crosstalk.pdf, perl=T)
		ls_rmd$png_crosstalk <- gsub("public/", "", output.file.crosstalk.png, perl=T)
		
		ls_rmd$xlsx_repurposing <- gsub("public/", "", output.file.repurposing, perl=T)
		ls_rmd$pdf_repurposing <- gsub("public/", "", output.file.repurposing.pdf, perl=T)
		ls_rmd$png_repurposing <- gsub("public/", "", output.file.repurposing.png, perl=T)
		
		output_dir <- gsub("cTCrosstalk.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_cTCrosstalk.Rmd", bookdown::html_document2(number_sections=F,theme="readable", hightlight="default"), output_dir=output_dir)
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
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

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", population=\"$population\", distance.max=\"$distance_max\", include.QTL=\"$include_QTL\", include.RGB=\"$include_RGB\", network=\"$network\", significance.threshold=\"$significance_threshold\", restart=\"$restart\", highlight.top=\"$highlight_top\", pathway.top=\"$pathway_top\", subnet.size=\"$subnet_size\", crosstalk.sig=\"$crosstalk_sig\", placeholder=\"$placeholder\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- cTCrosstalk: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
PIER_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			#############
			### for _priority.xlsx
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/_priority.xlsx/g;
			if(-e $tmp_file){
				$ajax_priority_xlsx_file=$tmp_file;
				$ajax_priority_xlsx_file=~s/^public//g;
				print STDERR "_priority.xlsx locates at $ajax_priority_xlsx_file\n";
			}
			
			#############
			### for _manhattan.pdf
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/_manhattan.pdf/g;
			if(-e $tmp_file){
				$ajax_manhattan_pdf_file=$tmp_file;
				$ajax_manhattan_pdf_file=~s/^public//g;
				print STDERR "_manhattan.pdf locates at $ajax_manhattan_pdf_file\n";
			}
			
			#############
			### for RMD_cTCrosstalk.html
			$tmp_file=$output_filename;
			$tmp_file=~s/cTCrosstalk.*//g;
			$tmp_file=$tmp_file."RMD_cTCrosstalk.html";
			if(-e $tmp_file){
				$ajax_rmd_html_file=$tmp_file;
				$ajax_rmd_html_file=~s/^public//g;
				print STDERR "RMD_cTCrosstalk locates at $ajax_rmd_html_file\n";
			}
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);
	
	# stash $xlsxfile;
	my $xlsxfile;
	$xlsxfile->{priority}=$ajax_priority_xlsx_file;
	$c->stash(xlsxfile => $xlsxfile);

	# stash $pdffile;
	my $pdffile;
	$pdffile->{manhattan}=$ajax_manhattan_pdf_file;
	$c->stash(pdffile => $pdffile);
	
  	$c->render();

}


# Render template "cTGene.html.ep"
sub PiER_cTGene {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $snplist = $c->param('snplist');
  	my $population = $c->param('pop') || 'NA'; # by default: NA
	my $distance_max = $c->param('distance');
	my $include_QTL = $c->param('QTL') || 'NA';
	my $include_RGB = $c->param('RGB') || 'NA';
  	my $network = $c->param('network') || 'STRING_highest'; # by default: STRING_highest
  	
	my $significance_threshold = $c->param('significance_threshold');
  	
	my $restart = $c->param('restart');
	my $highlight_top = $c->param('highlight_top');
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($snplist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $network.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.SNPs.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'cTGene.SNPs.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'cTGene.SNPs.'.$rand_file.'.r';
	
		my $my_input;
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $snplist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			# mac
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}elsif(-e '/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS'){
			# galahad
			$placeholder="/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS";
		}elsif(-e '/var/www/bigdata_dev'){
			# huawei
			$placeholder="/var/www/bigdata_dev";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
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
';

# for calling R function
$my_rscript.="
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

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", population=\"$population\", distance.max=\"$distance_max\", include.QTL=\"$include_QTL\", include.RGB=\"$include_RGB\", network=\"$network\", significance.threshold=\"$significance_threshold\", restart=\"$restart\", highlight.top=\"$highlight_top\", placeholder=\"$placeholder\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- cTGene: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
PIER_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			#############
			### for _priority.xlsx
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/_priority.xlsx/g;
			if(-e $tmp_file){
				$ajax_priority_xlsx_file=$tmp_file;
				$ajax_priority_xlsx_file=~s/^public//g;
				print STDERR "_priority.xlsx locates at $ajax_priority_xlsx_file\n";
			}
			
			#############
			### for _manhattan.pdf
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/_manhattan.pdf/g;
			if(-e $tmp_file){
				$ajax_manhattan_pdf_file=$tmp_file;
				$ajax_manhattan_pdf_file=~s/^public//g;
				print STDERR "_manhattan.pdf locates at $ajax_manhattan_pdf_file\n";
			}
			
			#############
			### for RMD_cTGene.html
			$tmp_file=$output_filename;
			$tmp_file=~s/cTGene.*//g;
			$tmp_file=$tmp_file."RMD_cTGene.html";
			if(-e $tmp_file){
				$ajax_rmd_html_file=$tmp_file;
				$ajax_rmd_html_file=~s/^public//g;
				print STDERR "RMD_cTGene locates at $ajax_rmd_html_file\n";
			}
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);
	
	# stash $xlsxfile;
	my $xlsxfile;
	$xlsxfile->{priority}=$ajax_priority_xlsx_file;
	$c->stash(xlsxfile => $xlsxfile);

	# stash $pdffile;
	my $pdffile;
	$pdffile->{manhattan}=$ajax_manhattan_pdf_file;
	$c->stash(pdffile => $pdffile);
	
  	$c->render();

}


# Render template "eV2CG.html.ep"
sub PiER_eV2CG {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $snplist = $c->param('snplist');
  	my $population = $c->param('pop') || 'NA'; # by default: NA
	my $crosslink = $c->param('crosslink') || 'proximity_10000';
  	my $network = $c->param('network') || 'STRING_highest'; # by default: STRING_highest
  	
	my $significance_threshold = $c->param('significance_threshold');
  	
	my $restart = $c->param('restart') || 0.7;
	my $highlight_top = $c->param('highlight_top') || 30;
	my $subnet_size = $c->param('subnet_size') || 30;
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($snplist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $network.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.SNPs.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'eV2CG.SNPs.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'eV2CG.SNPs.'.$rand_file.'.r';
	
		my $my_input="";
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $snplist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			# mac
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}elsif(-e '/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS'){
			# galahad
			$placeholder="/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS";
		}elsif(-e '/var/www/bigdata_dev'){
			# huawei
			$placeholder="/var/www/bigdata_dev";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
#/home/hfang/R-3.6.2/lib/R/library
# rm -rf /home/hfang/R-3.6.2/lib/R/library/00*
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
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
';

# for calling R function
$my_rscript.="
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

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", population=\"$population\", crosslink=\"$crosslink\", network=\"$network\", significance.threshold=\"$significance_threshold\", restart=\"$restart\", highlight.top=\"$highlight_top\", placeholder=\"$placeholder\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- eV2CG: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
PIER_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			#############
			### for RMD_eV2CG.html
			$tmp_file=$output_filename;
			$tmp_file=~s/eV2CG.*//g;
			$tmp_file=$tmp_file."RMD_eV2CG.html";
			if(-e $tmp_file){
				$ajax_rmd_html_file=$tmp_file;
				$ajax_rmd_html_file=~s/^public//g;
				print STDERR "RMD_eV2CG locates at $ajax_rmd_html_file\n";
			}
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}



# Render template "eCG2PG.html.ep"
sub PiER_eCG2PG {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $genelist = $c->param('genelist');
  	my $network = $c->param('network') || 'STRING_high'; # by default: STRING_highest
  	
	my $restart = $c->param('restart') || 0.7;
	my $highlight_top = $c->param('highlight_top') || 30;
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $network.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'eCG2PG.Genes.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'eCG2PG.Genes.'.$rand_file.'.r';
	
		my $my_input;
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			# mac
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}elsif(-e '/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS'){
			# galahad
			$placeholder="/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS";
		}elsif(-e '/var/www/bigdata_dev'){
			# huawei
			$placeholder="/var/www/bigdata_dev";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", network="", restart="", highlight.top="", placeholder="", ...){
	
	sT <- Sys.time()
	
	# read input file
	data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
	#data <- openxlsx::read.xlsx("/Users/hfang/Sites/XGR/PiER-site/app/examples/eV2CG.xlsx")
	#data <- openxlsx::read.xlsx("~/pier_app/public/app/examples/eV2CG.xlsx")
	
	restart <- as.numeric(restart)
	highlight.top <- as.numeric(highlight.top)
	
	GR.Gene <- oRDS("UCSC_knownGene", placeholder=placeholder)

	network.customised <- oDefineNet(network=network, STRING.only=c("experimental_score","database_score"), placeholder=placeholder) %>% oNetInduce(nodes_query=names(GR.Gene), knn=0, largest.comp=F)
	
	message(sprintf("Customised the network with %d nodes and %d edges (%s) ...", gorder(network.customised), gsize(network.customised), as.character(Sys.time())), appendLF=TRUE)
	
	message(sprintf("Preparing predictors (%s) ...", as.character(Sys.time())), appendLF=TRUE)

	pNode <- oPierGenes(data=data, network.customised=network.customised, restart=restart, placeholder=placeholder)

	if(class(pNode)=="pNode"){
		
		# *_CG2PG.txt
		pNode$priority %>% filter(priority>0) %>% write_delim(output.file, delim="\t")
		
		# *_CG2PG.xlsx
		output.file.CG2PG <- gsub(".txt$", "_CG2PG.xlsx", output.file, perl=T)
		pNode$priority %>% filter(priority>0) %>% openxlsx::write.xlsx(output.file.CG2PG)
		#pNode$priority %>% filter(priority>0) %>% transmute(gene=name, score=priority) %>% openxlsx::write.xlsx("/Users/hfang/Sites/XGR/PiER-site/app/examples/eCG2PG.xlsx")
		
		# Manhattan plot
		message(sprintf("Drawing manhattan plot (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		top.label.query <- pNode$priority %>% filter(priority>0) %>% pull(name)
		label.query.only <- T
		gp_manhattan_pNode <- oPierManhattan(pNode, color=c("darkred", "#2071b8"), point.size=0.2, top=highlight.top, top.label.type="text", top.label.size=2.5, top.label.col="black", top.label.query=top.label.query, label.query.only=label.query.only, y.lab="Gene affinity scores\n(scored 0-1)", signature=F, placeholder=placeholder)
		output.file.manhattan.pdf <- gsub(".txt$", "_manhattan.pdf", output.file, perl=T)
		ggsave(output.file.manhattan.pdf, gp_manhattan_pNode, device=cairo_pdf, width=9, height=3)
		output.file.manhattan.png <- gsub(".txt$", "_manhattan.png", output.file, perl=T)
		ggsave(output.file.manhattan.png, gp_manhattan_pNode, type="cairo", width=9, height=3)
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/PiER-site/pier_app/public
		## but outputs at public/tmp/eCG2PG.SNPs.STRING_high.72959383_priority.xlsx
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$runtime <- str_c(runtime," seconds")
		ls_rmd$data_input <- data %>% set_names(c("Genes","Weights"))
		ls_rmd$highlight_top <- highlight.top
		
		ls_rmd$xlsx_CG2PG <- gsub("public/", "", output.file.CG2PG, perl=T)
		ls_rmd$pdf_manhattan <- gsub("public/", "", output.file.manhattan.pdf, perl=T)
		ls_rmd$png_manhattan <- gsub("public/", "", output.file.manhattan.png, perl=T)
		
		output_dir <- gsub("eCG2PG.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_eCG2PG.Rmd", bookdown::html_document2(number_sections=F,theme="readable", hightlight="default"), output_dir=output_dir)
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
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

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", network=\"$network\", restart=\"$restart\", highlight.top=\"$highlight_top\", placeholder=\"$placeholder\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- eCG2PG: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
PIER_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			#############
			### for RMD_eCG2PG.html
			$tmp_file=$output_filename;
			$tmp_file=~s/eCG2PG.*//g;
			$tmp_file=$tmp_file."RMD_eCG2PG.html";
			if(-e $tmp_file){
				$ajax_rmd_html_file=$tmp_file;
				$ajax_rmd_html_file=~s/^public//g;
				print STDERR "RMD_eCG2PG locates at $ajax_rmd_html_file\n";
			}
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}


# Render template "eCrosstalk.html.ep"
sub PiER_eCrosstalk {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $genelist = $c->param('genelist');
	my $subnet_size = $c->param('subnet_size');
  	
	my $crosstalk_sig = $c->param('crosstalk_sig');
  	
  	# The output json file (default: '')
	my $ajax_txt_file='';
  	# The output html file (default: '')
	my $ajax_rmd_html_file='';
	
	# The output _priority.xlsx file (default: '')
	my $ajax_priority_xlsx_file='';
  	
	# The output _manhattan.pdf file (default: '')
	my $ajax_manhattan_pdf_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'eCrosstalk.Genes.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'eCrosstalk.Genes.'.$rand_file.'.r';
	
		my $my_input;
		my $line_counts=0;
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
			
			$line_counts++;
		}
		# at least two lines otherwise no $input_filename written
		if($line_counts >=2){
			PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			# mac
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}elsif(-e '/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS'){
			# galahad
			$placeholder="/data/archive/ULTRADDR/create_ultraDDR_database/dir_output_RDS";
		}elsif(-e '/var/www/bigdata_dev'){
			# huawei
			$placeholder="/var/www/bigdata_dev";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/home/hfang/R-3.6.2/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", subnet.size="", crosstalk.sig="", placeholder="", ...){
	
	sT <- Sys.time()
	
	# read input file
	data <- read_delim(input.file, delim="\t") %>% as.data.frame() %>% select(1:2)
	#data <- openxlsx::read.xlsx("/Users/hfang/Sites/XGR/PiER-site/app/examples/eCG2PG.xlsx")
	#data <- openxlsx::read.xlsx("~/pier_app/public/app/examples/eCG2PG.xlsx")
	
	subnet.size <- as.numeric(subnet.size)
	
	message(sprintf("Performing pathway crosstalk analysis restricted to %d network genes (%s) ...", subnet.size, as.character(Sys.time())), appendLF=TRUE)
	ig.KEGG.category <- oRDS("ig.KEGG.category", placeholder=placeholder)
	ls_ig <- ig.KEGG.category %>% filter(category %in% c("Organismal System","Environmental Process")) %>% deframe()
	ig <- oCombineNet(ls_ig)
	ig2 <- oNetInduce(ig, nodes_query=V(ig)$name, largest.comp=T) %>% as.undirected()
	
	df_priority <- tibble(name=data[,1], seed=1, weight=data[,2], priority=data[,2]) %>% as.data.frame()
	rownames(df_priority) <- df_priority$name
	pNode <- list(priority=df_priority[, c("seed","weight","priority")])
	class(pNode) <- "pNode"
	subg <- oPierSubnet(pNode, network=NA, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder)
	
	if(vcount(subg)>0){
	
		ind <- match(V(subg)$name, df_priority$name)
		V(subg)$priority <- df_priority$priority[ind]
	
		subg <- subg %>% oLayout(c("layout_with_kk","graphlayouts.layout_with_stress")[2])
	
		gp_rating <- oGGnetwork(g=subg, node.label="name", node.label.size=3, node.label.color="black", node.label.alpha=0.95, node.label.padding=0.5, node.label.arrow=0, node.label.force=0.4, node.shape=19, node.xcoord="xcoord", node.ycoord="ycoord", node.color="priority", node.color.title="Gene\nscores\n(input)", colormap="spectral", node.size.range=5, title="", edge.color="steelblue4", edge.color.alpha=0.5, edge.size=0.3, edge.curve=0.05)
		df_subg <- oIG2TB(subg,"nodes") %>% transmute(name,score=as.numeric(priority),description) %>% arrange(desc(score))
		
		# *_crosstalk.txt
		df_subg %>% write_delim(output.file, delim="\t")
		# *_crosstalk.xlsx
		output.file.crosstalk <- gsub(".txt$", "_crosstalk.xlsx", output.file, perl=T)
		df_subg %>% openxlsx::write.xlsx(output.file.crosstalk)

		# *_crosstalk.pdf *_crosstalk.png
		output.file.crosstalk.pdf <- gsub(".txt$", "_crosstalk.pdf", output.file, perl=T)
		ggsave(output.file.crosstalk.pdf, gp_rating, device=cairo_pdf, width=6, height=6)
		output.file.crosstalk.png <- gsub(".txt$", "_crosstalk.png", output.file, perl=T)
		ggsave(output.file.crosstalk.png, gp_rating, type="cairo", width=6, height=6)
		
		combinedP <- 1
		if(crosstalk.sig=="yes"){
			subg.sig <- oPierSubnet(pNode, network=NA, network.customised=ig2, subnet.size=subnet.size, placeholder=placeholder, test.permutation=T, num.permutation=10, respect=c("none","degree")[2], aggregateBy="fishers")
			combinedP <- signif(subg.sig$combinedP, digits=2)
		}
		
		######################################
		# RMD
		## R at /Users/hfang/Sites/XGR/PiER-site/pier_app/public
		## but outputs at public/tmp/eCrosstalk.SNPs.STRING_high.72959383_priority.xlsx
		######################################
		message(sprintf("RMD (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		
		eT <- Sys.time()
		runtime <- as.numeric(difftime(strptime(eT, "%Y-%m-%d %H:%M:%S"), strptime(sT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		
		ls_rmd <- list()
		ls_rmd$runtime <- str_c(runtime," seconds")
		ls_rmd$data_input <- data %>% set_names(c("Genes","Scores"))
		ls_rmd$combinedP <- combinedP
		
		ls_rmd$xlsx_crosstalk <- gsub("public/", "", output.file.crosstalk, perl=T)
		ls_rmd$pdf_crosstalk <- gsub("public/", "", output.file.crosstalk.pdf, perl=T)
		ls_rmd$png_crosstalk <- gsub("public/", "", output.file.crosstalk.png, perl=T)
		
		output_dir <- gsub("eCrosstalk.*", "", output.file, perl=T)
		
		## rmarkdown
		if(file.exists("/usr/local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
		}else if(file.exists("/home/hfang/.local/bin/pandoc")){
			Sys.setenv(RSTUDIO_PANDOC="/home/hfang/.local/bin")
		}else{
			message(sprintf("PANDOC is NOT FOUND (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
		rmarkdown::render("public/RMD_eCrosstalk.Rmd", bookdown::html_document2(number_sections=F,theme="readable", hightlight="default"), output_dir=output_dir)
	}
	
	##########################################
}
';

# for calling R function
$my_rscript.="
startT <- Sys.time()

library(tidyverse)
library(igraph)

# galahad
vec <- list.files(path='/home/hfang/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# huawei
vec <- list.files(path='/root/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))
# mac
vec <- list.files(path='/Users/hfang/Sites/XGR/Fang/R', pattern='.r', full.names=T)
ls_tmp <- lapply(vec, function(x) source(x))

R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", subnet.size=\"$subnet_size\", crosstalk.sig=\"$crosstalk_sig\", placeholder=\"$placeholder\")

endT <- Sys.time()
runTime <- as.numeric(difftime(strptime(endT, '%Y-%m-%d %H:%M:%S'), strptime(startT, '%Y-%m-%d %H:%M:%S'), units='secs'))
message(str_c('\n--- eCrosstalk: ',runTime,' secs ---\n'), appendLF=TRUE)
";

# for calling R function
PIER_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
# $input_filename (and $rscript_filename) must exist
if(-e $rscript_filename and -e $input_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.6.2/bin/Rscript'){
    	# galahad
    	$command="/home/hfang/R-3.6.2/bin/Rscript $rscript_filename";
    }else{
    	# mac and huawei
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			my $tmp_file='';
			
			## notes: replace 'public/' with '/'
			$tmp_file=$output_filename;
			if(-e $tmp_file){
				$ajax_txt_file=$tmp_file;
				$ajax_txt_file=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file\n";
			}
			
			#############
			### for RMD_eCrosstalk.html
			$tmp_file=$output_filename;
			$tmp_file=~s/eCrosstalk.*//g;
			$tmp_file=$tmp_file."RMD_eCrosstalk.html";
			if(-e $tmp_file){
				$ajax_rmd_html_file=$tmp_file;
				$ajax_rmd_html_file=~s/^public//g;
				print STDERR "RMD_eCrosstalk locates at $ajax_rmd_html_file\n";
			}
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_txt_file;
	$c->stash(ajax_txt_file => $ajax_txt_file);
	
	# stash $ajax_rmd_html_file;
	$c->stash(ajax_rmd_html_file => $ajax_rmd_html_file);

	
  	$c->render();

}

1;
