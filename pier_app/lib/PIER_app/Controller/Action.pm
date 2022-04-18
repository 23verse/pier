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

# Render template "A2Genes.html.ep"
sub A2_genes {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $genelist = $c->param('genelist');
  	my $uploadfile = $c->req->upload('uploadfile');
	if(defined($uploadfile)) {
		$uploadfile = $uploadfile->slurp;
	}else{
		$uploadfile='';
	}
	
	my $size_range_min = $c->param('size_range_min') || 1;
	my $size_range_max = $c->param('size_range_max') || 1000000;
	my $min_overlap = $c->param('min_overlap') || 1;
	
	my $test = $c->param('test');
	my $vis = $c->param('vis');
  	
  	my $ontology="AA";

  	# The output txt file (default: '')
	my $ajax_txt_file='';
  	
  	# The output json file (default: '')
	my $ajax_json_file='';
  	
	# The output pdf file (default: '')
	my $ajax_pdf_file='';
	# The output html file (default: '')
	my $ajax_html_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $ontology.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'enrichment.Genes.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'enrichment.Genes.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
		}
		PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $background_filename="";
		if($uploadfile ne ''){
			my $my_background;
			foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $uploadfile)) {
				next if($line=~/^\s*$/);
				$my_background.=$line."\n";
			}
			$background_filename=$tmpFolder.'/'.'background.'.$rand_file.'.txt';
			PIER_app::Controller::Utils::export_to_file($background_filename, $my_background);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}else{
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
R_pipeline <- function (input.file="", background.file="", output.file="", size_range_min="", size_range_max="", min_overlap="", test="", vis="", RData.location="", ...){
	
	#############################
	## an example for coding test
	if(0){
		input.file <- "/Users/hfang/Sites/XGR/A2-site/a2_app/public/app/examples/PDL1.2columns.MI.txt"
		data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1]
		
		data("Haploid_regulators")
		## for "PDL1"
		ind <- grepl("PDL1", Haploid_regulators$Phenotype)
		df <- Haploid_regulators[ind,c("Gene","MI","FDR")]
		data <- df$Gene
		
		RData.location <- "/Users/hfang/Sites/SVN/github/bigdata_dev"
		ls_eTerm <- xEnricherGenesAdv(list_vec=data, ontologies="AA", size.range=c(10,2000), min.overlap=5, test="fisher", RData.location=RData.location)
		
		output.file <- "public/tmp/enrichment.Genes.AA.33897825.txt"
	}
	#############################
		
	# read input file
	data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1]
	
	if(background.file!=""){
		# read background file
		background <- read.delim(file=background.file, header=F, stringsAsFactors=F)[,1]	
	}else{
		background <- NULL
	}
	
	size.range <- as.numeric(c(size_range_min, size_range_max))
	min.overlap <- as.numeric(min_overlap)
	
	# perform enrichment analysis
	ls_eTerm <- xEnricherGenesAdv(list_vec=data, background=background, ontologies="AA", size.range=size.range, min.overlap=min.overlap, test=test, RData.location=RData.location, plot=F, ...)
	
	if(class(ls_eTerm)=="ls_eTerm"){
		# save enrichment results to the output file
		#res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
		res <- ls_eTerm$df
		res$CI <- paste0("[",res$CIl,",",res$CIu,"]")
		### by default
		res_f <- data.frame(term=rownames(res), res, link=rep("",nrow(res)))
		### for AA only
		if(vis!="no"){
			passflag <- input.file
			passflag <- gsub("^public/tmp/","",passflag)
			passflag <- gsub(".txt$","",passflag)
			
			y <- read.delim(file=input.file, header=F, stringsAsFactors=F)
			if(ncol(y)>=2){
				y[,2] <- as.numeric(y[,2])
				y <- y[!is.na(y[,2]),]
				z <- y[,2]
				if(max(z) <=0){
					colormap <- "deepskyblue-white"
				}else if(min(z) >=0){
					colormap <- "white-darkorange"
				}else{
					colormap <- "deepskyblue-white-darkorange"
				}
			}else{
				colormap <- "NULL"
			}
			
			if(0){
				## de novo graphml
				passflag <- paste0(passflag,"_",colormap)
				#####
				### passflag is: input.file "_" colormap
				#####
				link <- paste("/A2/explorer/", ls_eTerm$df$id, "/", passflag, sep="")
				link[ls_eTerm$df$adjp >= 0.05] <- ""
				res_f <- data.frame(term=rownames(res), res, link=link)
			
			}else{
				## pre-created graphml
				data <- y
				if(ncol(data)==2){
					colnames(data) <- c("label","lfc")
				}else if(ncol(data)==1){
					data <- cbind(data,rep(1,nrow(data)))
					colnames(data) <- c("label","lfc")
				}else{
					colnames(data) <- c("label","lfc","fdr")
				}
				
				if(colormap=="NULL"){
					z <- data[,2]
					if(max(z) <=0){
						colormap <- "deepskyblue-white"
					}else if(min(z) >=0){
						colormap <- "white-darkorange"
					}else{
						colormap <- "deepskyblue-white-darkorange"
					}
				}
				
				out <- passflag
				out <- gsub("data.Genes.AA.","",out)
				#####
				### out is: $rand_flag
				#####
				df_enrichment <- ls_eTerm$df
				if(vis=="yes_sig"){
					df_enrichment <- df_enrichment[df_enrichment$adjp<0.05, ]
				}
				ls_tmp <- lapply(1:nrow(df_enrichment), function(i){
					query <- df_enrichment$id[i]
					out.file <- paste0("public/",out,".",query,".graphml")
					xGraphML2AA(data=data, query=query, node.label="label", node.color="lfc", colormap=colormap, node.highlight="fdr", filename=out.file, RData.location=RData.location, verbose=F)
				})
				
				link <- paste("/A2/explorer/", ls_eTerm$df$id, "/", out, sep="")
				if(vis=="yes_sig"){
					link[ls_eTerm$df$adjp >= 0.05] <- ""
				}
				res_f <- data.frame(term=rownames(res), res, link=link)
			}
			
		}
		###
		utils::write.table(res_f[,c("id","name","adjp","zscore","or","CI","nOverlap","members_Overlap")], file=output.file, sep="\t", row.names=FALSE)
	
		# save to the json file
		res_f <- jsonlite::toJSON(res_f, pretty=T, digits=10)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file, perl=T)
		base::write(res_f, file=output.file.json)
		
		###########
		## CTree plot of enrichment results
		###########
		AA.template <- xRDataLoader("AA.template", RData.location=RData.location)
		# consensus tree
		ig <- AA.template$consensus$ig

		# output file ".pdf"
		gp <- xEnrichCtree(ls_eTerm, ig, node.size.range=c(0.5,5.5), leave.label.expansion=1.15, leave.label.size=3, leave.label.color="steelblue", limit.expansion=2, edge.width=0.5, leave.label.orientation=c("inwards","outwards")[2])
		if(!is.null(gp)){
			output.file.pdf <- gsub(".txt$", ".pdf", output.file, perl=T)
			width <- 6.5
			height <- 6.5*1.1
			pdf(output.file.pdf, width=width, height=height, compress=T)
			print(gp)
			dev.off()
		}
		
		# output file ".html"
		gp_atlas <- xEnrichCtree(ls_eTerm, ig, node.size.range=c(0.5,5.5), leave.label.expansion=1.15, leave.label.size=3, leave.label.color="steelblue", limit.expansion=2, edge.width=0.5, leave.label.orientation=c("inwards","outwards")[2])
		df <- subset(gp_atlas$data,leaf==T)
		#33897825.AA:hsa04310.graphml
		df$onclick <- sprintf("window.open(\"/A2/explorer/%s/%s\")", df$id, out)	
		#### only enriched with onclick event
		ind <- match(df$id, df_enrichment$id)
		df$onclick[is.na(ind)] <- ""
		####		
		df$tooltip <- paste0("- ", df$id, "\n- ", df$name, "\n- ", df$subcategory)
		gg <- gp_atlas + ggiraph::geom_point_interactive(data=df, aes(x=x, y=y, tooltip=tooltip, data_id=name, onclick=onclick), size=5, alpha=0.01, color="white")
		gr <- ggiraph::girafe(ggobj=gg, width_svg=6, height_svg=5, width=0.6)
		tooltip_css <- "background-color:white;color:black;font-size:10pt;font-style:normal;font-weight:bold;padding:10px;border-radius:5px;"
		hover_css <- "fill:orange;"
		gr <- ggiraph::girafe_options(gr, ggiraph::opts_tooltip(css=tooltip_css, opacity=0.75, use_fill=T), ggiraph::opts_hover(css=hover_css), ggiraph::opts_toolbar(position="topright", saveaspng=F))
		if(!is.null(gp_atlas)){
			output.file.html <- gsub(".txt$", ".html", output.file, perl=T)
			output.file.dir <- gsub(".txt$", "_files", output.file, perl=T)
			htmlwidgets::saveWidget(widget=gr, file=basename(output.file.html), background="white", selfcontained=F)
			if(base::file.exists(basename(output.file.html))){
				base::file.copy(from=basename(output.file.html), to=dirname(output.file.html), overwrite=T)
				file.remove(basename(output.file.html))
				
				base::file.copy(from=basename(output.file.dir), to=dirname(output.file.dir), overwrite=T, recursive=T)
				
			}

		}

	}
}
';

# for calling R function
$my_rscript.="
library(XGR)
library(jsonlite)
R_pipeline(input.file=\"$input_filename\", background.file=\"$background_filename\", output.file=\"$output_filename\", size_range_min=\"$size_range_min\", size_range_max=\"$size_range_max\", min_overlap=\"$min_overlap\", test=\"$test\", vis=\"$vis\", RData.location=\"$placeholder\")
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
			### for json
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.json/g;
			if(-e $tmp_file){
				$ajax_json_file=$tmp_file;
				$ajax_json_file=~s/^public//g;
				print STDERR "JSON locates at $ajax_json_file\n";
			}
			
			#############
			### for pdf
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.pdf/g;
			if(-e $tmp_file){
				$ajax_pdf_file=$tmp_file;
				$ajax_pdf_file=~s/^public//g;
				print STDERR "PDF locates at $ajax_pdf_file\n";
			}
			### for html
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.html/g;
			if(-e $tmp_file){
				$ajax_html_file=$tmp_file;
				$ajax_html_file=~s/^public//g;
				print STDERR "HTML locates at $ajax_html_file\n";
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
	
	# stash $ajax_json_file;
	$c->stash(ajax_json_file => $ajax_json_file);
	
	# stash $pdffile;
	my $pdffile;
	$pdffile->{ctree}=$ajax_pdf_file;
	$c->stash(pdffile => $pdffile);
	
	# stash $htmlfile;
	my $htmlfile;
	$c->stash(htmlfile => $ajax_html_file);
	
	my $num_genes=0;
	if(defined($genelist)){
		my @lines = split(/\r\n|\n/, $genelist);
		my %rec_genes;
		foreach my $line (@lines) {
			next if($line=~/^\s*$/);
			my $rec;
			$rec->{ID}=$line;
			$rec->{Name}=$line;
			$rec_genes{$line}=$rec;
		}
		$c->stash(rec_genes => \%rec_genes);
		$num_genes=scalar(keys %rec_genes);
	}
	
	$c->stash(num_genes => $num_genes);
	$c->stash(ontology => $ontology);
	
  	$c->render();

}

# Render template "A2Genes.html.ep"
sub A2_genes_bk {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $ontology = $c->param('ontology') || 'DO'; # by default: DO
  	my $genelist = $c->param('genelist');
  	my $uploadfile = $c->req->upload('uploadfile');
	if(defined($uploadfile)) {
		$uploadfile = $uploadfile->slurp;
	}else{
		$uploadfile='';
	}
	
	my $size_range_min = $c->param('size_range_min') || 1;
	my $size_range_max = $c->param('size_range_max') || 1000000;
	my $min_overlap = $c->param('min_overlap') || 1;
	
	my $test = $c->param('test');
	my $vis = $c->param('vis');
  	
  	# The output json file (default: '')
	my $ajax_json_file='';
  	
	# The output pdf file (default: '')
	my $ajax_pdf_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $ontology.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'enrichment.Genes.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'enrichment.Genes.'.$rand_file.'.r';
	
		my $my_input;
		if(0){
		foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $genelist)) {
			next if($line=~/^\s*$/);
			$my_input.=$line."\n";
		}
		}else{
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
		}
		}
		PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $background_filename="";
		if($uploadfile ne ''){
			my $my_background;
			foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $uploadfile)) {
				next if($line=~/^\s*$/);
				$my_background.=$line."\n";
			}
			$background_filename=$tmpFolder.'/'.'background.'.$rand_file.'.txt';
			PIER_app::Controller::Utils::export_to_file($background_filename, $my_background);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}else{
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
R_pipeline <- function (input.file="", background.file="", output.file="", ontology="", size_range_min="", size_range_max="", min_overlap="", test="", vis="", RData.location="", ...){
	
	#############################
	## an example for coding test
	if(0){
		input.file <- "/Users/hfang/Sites/XGR/A2-site/a2_app/public/app/examples/PDL1.2columns.MI.txt"
		data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1]
		
		data("Haploid_regulators")
		## for "PDL1"
		ind <- grepl("PDL1", Haploid_regulators$Phenotype)
		df <- Haploid_regulators[ind,c("Gene","MI","FDR")]
		data <- df$Gene
		
		RData.location <- "/Users/hfang/Sites/SVN/github/bigdata_dev"
		ls_eTerm <- xA2EnricherGenes(list_vec=data, ontologies="AA", size.range=c(10,2000), min.overlap=5, test="fisher", RData.location=RData.location)
	}
	#############################
		
	# read input file
	data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1]
	
	if(background.file!=""){
		# read background file
		background <- read.delim(file=background.file, header=F, stringsAsFactors=F)[,1]	
	}else{
		background <- NULL
	}
	
	size.range <- as.numeric(c(size_range_min, size_range_max))
	min.overlap <- as.numeric(min_overlap)
	
	# perform enrichment analysis
	ls_eTerm <- xA2EnricherGenes(list_vec=data, background=background, ontologies=ontology, size.range=size.range, min.overlap=min.overlap, test=test, RData.location=RData.location, plot=F, ...)
	
	if(class(ls_eTerm)=="ls_eTerm"){
		# save enrichment results to the output file
		#res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
		res <- ls_eTerm$df
		res$CI <- paste0("[",res$CIl,",",res$CIu,"]")
		### by default
		res_f <- data.frame(term=rownames(res), res, link=rep("",nrow(res)))
		### for AA only
		if(ontology == "AA" & vis!="no"){
			passflag <- input.file
			passflag <- gsub("^public/tmp/","",passflag)
			passflag <- gsub(".txt$","",passflag)
			
			y <- read.delim(file=input.file, header=F, stringsAsFactors=F)
			if(ncol(y)>=2){
				y[,2] <- as.numeric(y[,2])
				y <- y[!is.na(y[,2]),]
				z <- y[,2]
				if(max(z) <=0){
					colormap <- "deepskyblue-white"
				}else if(min(z) >=0){
					colormap <- "white-darkorange"
				}else{
					colormap <- "deepskyblue-white-darkorange"
				}
			}else{
				colormap <- "NULL"
			}
			
			if(0){
				## de novo graphml
				passflag <- paste0(passflag,"_",colormap)
				#####
				### passflag is: input.file "_" colormap
				#####
				link <- paste("/A2/explorer/", ls_eTerm$df$id, "/", passflag, sep="")
				link[ls_eTerm$df$adjp >= 0.05] <- ""
				res_f <- data.frame(term=rownames(res), res, link=link)
			
			}else{
				## pre-created graphml
				data <- y
				if(ncol(data)==2){
					colnames(data) <- c("label","lfc")
				}else if(ncol(data)==1){
					data <- cbind(data,rep(1,nrow(data)))
					colnames(data) <- c("label","lfc")
				}else{
					colnames(data) <- c("label","lfc","fdr")
				}
				
				if(colormap=="NULL"){
					z <- data[,2]
					if(max(z) <=0){
						colormap <- "deepskyblue-white"
					}else if(min(z) >=0){
						colormap <- "white-darkorange"
					}else{
						colormap <- "deepskyblue-white-darkorange"
					}
				}
				
				out <- passflag
				out <- gsub("data.Genes.AA.","",out)
				#####
				### out is: $rand_flag
				#####
				df <- ls_eTerm$df
				if(vis=="yes_sig"){
					df <- df[df$adjp<0.05, ]
				}
				ls_tmp <- lapply(1:nrow(df), function(i){
					query <- df$id[i]
					out.file <- paste0("public/",out,".",query,".graphml")
					xGraphML2AA(data=data, query=query, node.label="label", node.color="lfc", colormap=colormap, node.highlight="fdr", filename=out.file, RData.location=RData.location, verbose=F)
				})
				
				link <- paste("/A2/explorer/", ls_eTerm$df$id, "/", out, sep="")
				if(vis=="yes_sig"){
					link[ls_eTerm$df$adjp >= 0.05] <- ""
				}
				res_f <- data.frame(term=rownames(res), res, link=link)
			}
			
		}
		###
		utils::write.table(res_f, file=output.file, sep="\t", row.names=FALSE)
	
		# save to the json file
		res_f <- jsonlite::toJSON(res_f, pretty=T, digits=10)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file, perl=T)
		base::write(res_f, file=output.file.json)
		
		# output file ".pdf"
		## forest plot of enrichment results
		gp <- xEnrichForest(ls_eTerm, top_num=10, FDR.cutoff=0.05, CI.one=T, signature=F)
		if(!is.null(gp)){
			output.file.pdf <- gsub(".txt$", ".pdf", output.file, perl=T)
			width <- 7
			height <- 2 + nrow(gp$data)*0.1
			pdf(output.file.pdf, width=width, height=height, compress=T)
			print(gp)
			dev.off()
		}

	}
}
';

# for calling R function
$my_rscript.="
library(XGR)
library(jsonlite)
R_pipeline(input.file=\"$input_filename\", background.file=\"$background_filename\", output.file=\"$output_filename\", ontology=\"$ontology\", size_range_min=\"$size_range_min\", size_range_max=\"$size_range_max\", min_overlap=\"$min_overlap\", test=\"$test\", vis=\"$vis\", RData.location=\"$placeholder\")
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
			$tmp_file=~s/.txt$/.json/g;
			if(-e $tmp_file){
				$ajax_json_file=$tmp_file;
				$ajax_json_file=~s/^public//g;
				print STDERR "JSON locates at $ajax_json_file\n";
			}
			
			#############
			### for pdf
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.pdf/g;
			if(-e $tmp_file){
				$ajax_pdf_file=$tmp_file;
				$ajax_pdf_file=~s/^public//g;
				print STDERR "PDF locates at $ajax_pdf_file\n";
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
	
	# stash $ajax_json_file;
	$c->stash(ajax_json_file => $ajax_json_file);
	
	# stash $pdffile;
	my $pdffile;
	$pdffile->{forest}=$ajax_pdf_file;
	$c->stash(pdffile => $pdffile);
	
	my $num_genes=0;
	if(defined($genelist)){
		my @lines = split(/\r\n|\n/, $genelist);
		my %rec_genes;
		foreach my $line (@lines) {
			next if($line=~/^\s*$/);
			my $rec;
			$rec->{ID}=$line;
			$rec->{Name}=$line;
			$rec_genes{$line}=$rec;
		}
		$c->stash(rec_genes => \%rec_genes);
		$num_genes=scalar(keys %rec_genes);
	}
	
	$c->stash(num_genes => $num_genes);
	$c->stash(ontology => $ontology);
	
  	$c->render();

}



sub A2_mapper {
  	my $c = shift;
	
	my $dbh = PIER_app::Controller::Utils::DBConnect('pi');
	my $sth;
	
	#########
	# curated
	#########
	$sth = $dbh->prepare("select distinct path_hsa from KEGG_nodes where path_hsa like 'pi:%';");
	$sth->execute();
	my @data_hsa;
	my $i=0;
	my $switch_running=0;
	if($sth->rows!=0){
		while(my @row = $sth->fetchrow_array){
			$i++;
			my $query=$row[0];
			$query=~s/^pi/path/g;
			my $sth1 = $dbh->prepare("SELECT distinct path_hsa, name, subcategory, category FROM KEGG where path_hsa=?;");
			$sth1->execute($query);
			while(my @row1 = $sth1->fetchrow_array){
				my $rec;
				$rec->{path_hsa}=$row[0];
				$rec->{path_hsa}=~s/pi/AA/g;
				$rec->{name}=$row1[1];
				$rec->{category}=$row1[2];
				$rec->{subcategory}=$row1[3];
				if($i==1){
					$switch_running++;
					$rec->{switch}=$switch_running;
				}elsif($i==$sth1->rows){
					$switch_running++;
					$rec->{switch}=$switch_running;
				}else{
					$rec->{switch}=0;
				}
				push @data_hsa,$rec;
			}
		}
	}
	#########
	
	if(0){
	#########
	# automatic
	#########
	$sth = $dbh->prepare("SELECT distinct path_hsa, name, subcategory, category FROM KEGG where (subcategory='Organismal Systems' and category='Immune system') OR (subcategory='Environmental Information Processing') OR (subcategory='Cellular Processes') order by subcategory desc,category desc;");
	$sth->execute();
	my $pre=0;
	my $i=0;
	while (my @row = $sth->fetchrow_array) {
	
		my $switch=0;
		$i++;
		if($i==1){
			$switch_running++;
			$switch=$switch_running;
		}elsif($i==$sth->rows){
			#$switch_running++;
			#$switch=$switch_running;
		}else{
			if($pre ne $row[2]){
				$switch_running++;
				$switch=$switch_running;
			}
		}
		$pre=$row[2];
	
		my $rec;
		$rec->{path_hsa}=$row[0];
		$rec->{name}=$row[1];
		$rec->{category}=$row[2];
		$rec->{subcategory}=$row[3];
		$rec->{switch}=$switch;
		push @data_hsa,$rec;
	}
	#########
	}
	
	$sth->finish();
	$c->stash(rec_pathway => \@data_hsa);
	
	PIER_app::Controller::Utils::DBDisconnect($dbh);
	
	print STDERR "#path:".scalar(@data_hsa)."\n";
	
	my $hsa = $c->param('hsa') || '';	
	my $colormap = $c->param('colormap') || 'NULL';
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
  	my $genelist = $c->param('genelist') || '';
  	if(defined($genelist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = 'input'.(int rand 99999999);
		my $rand_file = $rand_flag.'.'.$hsa;
		my $input_filename=$tmpFolder.'/'.$rand_flag.'.txt';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
		}
		PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		if($hsa ne '' and $colormap ne ''){
			# redirect_to: used to redirect response
			my $pass_flag = $rand_flag.'_'.$colormap;
			$c->redirect_to("/A2/explorer/$hsa/$pass_flag");
		}

	}else{
		
		if($hsa ne ''){
			$c->redirect_to("/A2/explorer/$hsa");
		}
		
	}
  	
	my $num_genes=0;
	if(defined($genelist)){
		my @lines = split(/\r\n|\n/, $genelist);
		my %rec_genes;
		foreach my $line (@lines) {
			next if($line=~/^\s*$/);
			my $rec;
			$rec->{ID}=$line;
			$rec->{Name}=$line;
			$rec_genes{$line}=$rec;
		}  	
		$c->stash(rec_genes => \%rec_genes);
		$num_genes=scalar(keys %rec_genes);
	}
	
	$c->stash(num_genes => $num_genes);
	$c->stash(hsa => $hsa);
  	
  	$c->render();
}

sub A2_explorer_hsa {
	my $c = shift;
	
	my $hsa= $c->param("hsa");
	my $passflag= $c->param("passflag") || '';
	
	## to check whether graphml exists, if so, just loading
	my $file_exist="public/".$passflag.".".$hsa.".graphml";
	if(-e $file_exist){
		print STDERR "Just load the existing graphml file!\n";
		
		my $graphml_filename_colors = $file_exist;
		print STDERR "$graphml_filename_colors\n";
		$graphml_filename_colors=~s/^public//g;
		print STDERR "$graphml_filename_colors\n";
		## for graphml file downloading
		my $graphml_filename_download=$graphml_filename_colors;
		$graphml_filename_download=~s/^\///g;
		print STDERR "$graphml_filename_download\n";
		#################################################

		$c->stash(graphml_filename_colors => $graphml_filename_colors);
		$c->stash(graphml_filename_download => $graphml_filename_download);
		
		$c->render();
		
	}else{
	
		####################################
		my $hsa_now='';
		my $dbh = PIER_app::Controller::Utils::DBConnect('pi');
		my $sth;
	
		$sth = $dbh->prepare( "SELECT path_hsa, name, category, subcategory, source, target, source_hsa, target_hsa FROM KEGG WHERE path_hsa=?;" );
		$hsa_now=$hsa;
		$hsa_now=~s/pi/path/g;
		$hsa_now=~s/AA/path/g;
		$sth->execute($hsa_now);
		my $hsa_data=$sth->fetchrow_hashref;
		if(!$hsa_data->{path_hsa}){
			return $c->reply->not_found
		}else{
			$hsa_data->{target}=~s/NA//g;
			$hsa_data->{path_hsa}=~s/path://g;
			$hsa_data->{path_hsa}=~s/pi://g;
			$hsa_data->{path_hsa}=~s/AA://g;
			$hsa_data->{path_hsa}=~s/AT://g;
			$c->stash(hsa_data => $hsa_data);
		}
		$sth->finish;
		PIER_app::Controller::Utils::DBDisconnect($dbh);
		####################################
	
		my $data="NULL";
		my $colormap="wyr";
		if($passflag ne ''){
			my $file='';
			($file, $colormap)=split(/\_/,$passflag);
			
			########################
			# because of automatical action: GET "/A2/explorer/AA:hsa04659/crossdomain.xml"
			# so disable it
			########################
			if($file=~/^crossdomain/){
				return $c->reply->not_found;
			}
			########################
		
			if(!defined($colormap)){
				$colormap="wyr";
			}
		
			my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
			my $input_filename=$tmpFolder.'/'.$file.'.txt';;
			print STDERR "$input_filename\t$colormap\n";
			if(-e $input_filename){
				print STDERR "$input_filename\n";
				$data=$input_filename;
			}
		}
	
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}else{
			$placeholder="/var/www/bigdata_dev";
		}
		
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $rand_flag.'.'.$hsa;
		my $rscript_filename=$tmpFolder.'/'.$rand_file.'.r';
		my $graphml_filename_colors='public/'.$rand_file.'.graphml';	
		
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
		R_pipeline <- function (data="", query="", colormap="", output.file="", RData.location="", ...){

			if(data == "NULL"){
				message(paste(c("\nthe input file: ",as.character(data)), collapse=""), appendLF=TRUE)
				data <- NULL
			}else{
				message(paste(c("\nread the input file: ",as.character(data)), collapse=""), appendLF=TRUE)
		
				data <- utils::read.delim(file=data, header=T, row.names=NULL, stringsAsFactors=F)
		
				if(ncol(data)==2){
					colnames(data) <- c("label","lfc")
				}else if(ncol(data)==1){
					data <- cbind(data,rep(1,nrow(data)))
					colnames(data) <- c("label","lfc")
				}else{
					colnames(data) <- c("label","lfc","fdr")
				}
		
		
				message(paste(c(colnames(data)), collapse=""), appendLF=TRUE)
		
				if(colormap=="NULL"){
					z <- as.numeric(data[,2])
					z <- z[!is.na(z)]
					if(max(z) <=0){
						colormap <- "deepskyblue-lightyellow"
					}else if(min(z) >=0){
						colormap <- "lightyellow-darkorange"
					}else{
						colormap <- "deepskyblue-lightyellow-darkorange"
					}
				}
		
				message(paste(colormap, collapse=""), appendLF=TRUE)
		
			}

			xGraphML2AA(data=data, query=query, node.label="label", node.color="lfc", colormap=colormap, node.highlight="fdr", filename=output.file, RData.location=RData.location, ...)
		}
		';
	
		# for calling R function
		$my_rscript.="
		library(XGR)
		library(jsonlite)
		R_pipeline(data=\"$data\", query=\"$hsa_now\", colormap=\"$colormap\", output.file=\"$graphml_filename_colors\", RData.location=\"$placeholder\")
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
		
				if(! -e $graphml_filename_colors){
					print STDERR "Cannot find $graphml_filename_colors\n";
				}else{
					
					print STDERR "Create a new graphml file!\n";
				
					print STDERR "$graphml_filename_colors\n";
					$graphml_filename_colors=~s/^public//g;
					print STDERR "$graphml_filename_colors\n";
					## for graphml file downloading
					my $graphml_filename_download=$graphml_filename_colors;
					$graphml_filename_download=~s/^\///g;
					print STDERR "$graphml_filename_download\n";
					#################################################
					
					$c->stash(graphml_filename_colors => $graphml_filename_colors);
					$c->stash(graphml_filename_download => $graphml_filename_download);
		
					$c->render();
				}
			}
		}else{
			print STDERR "Cannot find $rscript_filename\n";
		}
		##########################################
		# END: R
		##########################################
	
	}
}

sub A2_viewer {
	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}

	my $example_graphml_file='public/app/examples/AKT_mTOR.graphml';

  	my $uploadfile = $c->req->upload('uploadfile') || '';
  	my $example_graphml = $c->param('example_graphml') || '';
	
	if($uploadfile=~/Mojo::Upload/ && $example_graphml ne '') {
		if(defined($example_graphml)){
			$uploadfile = PIER_app::Controller::Utils::read_from_file($example_graphml_file);
		}
	}elsif($uploadfile=~/Mojo::Upload/){
		$uploadfile = $uploadfile->slurp;
	}

  	# The output graph file (default: '')
	my $graphml_filename_colors='';
	my $graphml_filename_download='';
	
	if($uploadfile ne ''){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		$graphml_filename_colors=$tmpFolder.'/'.$rand_flag.'.graphml';
		PIER_app::Controller::Utils::export_to_file($graphml_filename_colors, $uploadfile);
	}
	
	## to check whether graphml exists, if so, just loading
	if(-e $graphml_filename_colors){
		print STDERR "Just load the existing graphml file!\n";
		
		print STDERR "$graphml_filename_colors\n";
		$graphml_filename_colors=~s/^public//g;
		print STDERR "$graphml_filename_colors\n";
		## for graphml file downloading
		$graphml_filename_download=$graphml_filename_colors;
		$graphml_filename_download=~s/^\///g;
		print STDERR "$graphml_filename_download\n";
	}
	
	$c->stash(graphml_filename_colors => $graphml_filename_colors);
	$c->stash(graphml_filename_download => $graphml_filename_download);
	
	$c->render();
}


# Render template "A2Regions.html.ep"
sub A2_regions {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $genelist = $c->param('genelist');
	my $build_conversion = $c->param('build');
	my $crosslink = $c->param('crosslink');
	my $crosslink_top = $c->param('crosslink_top');
	
	my $distance_max = $c->param('distance') || 50000;
	my $decay_kernel = $c->param('kernel') || 'rapid';
  	
  	my $uploadfile = $c->req->upload('uploadfile');
	if(defined($uploadfile)) {
		$uploadfile = $uploadfile->slurp;
	}else{
		$uploadfile='';
	}
	
	my $size_range_min = $c->param('size_range_min') || 1;
	my $size_range_max = $c->param('size_range_max') || 1000000;
	my $min_overlap = $c->param('min_overlap') || 1;
	
	my $test = $c->param('test');
	my $vis = $c->param('vis');
  	
  	my $ontology="AA";
  	
	my $output_xGenes = $c->param('output_xGenes');
  	
  	# The output txt file (default: '')
	my $ajax_txt_file='';
  	
  	# The output json file (default: '')
	my $ajax_json_file='';
  	
	# The output pdf file (default: '')
	my $ajax_pdf_file='';
	my $ajax_pdf_file_xGenes_GR='';
		
	# The output txt file (default: '')
	my $ajax_txt_file_xGenes='';
	my $ajax_txt_file_xGenes_GR='';
  	
	# The output html file (default: '')
	my $ajax_html_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $ontology.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Regions.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'enrichment.Regions.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'enrichment.Regions.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
		}
		PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $background_filename="";
		if($uploadfile ne ''){
			my $my_background;
			foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $uploadfile)) {
				next if($line=~/^\s*$/);
				$my_background.=$line."\n";
			}
			$background_filename=$tmpFolder.'/'.'background.'.$rand_file.'.txt';
			PIER_app::Controller::Utils::export_to_file($background_filename, $my_background);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}else{
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
R_pipeline <- function (input.file="", background.file="", output.file="", build.conversion="", crosslink="", crosslink_top="", distance.max="", decay.kernel="", ontology="", size_range_min="", size_range_max="", min_overlap="", test="", vis="", output_xGenes="", RData.location="", ...){
	
	#############################
	## an example for coding test
	if(0){
		RData.location <- "/Users/hfang/Sites/SVN/github/bigdata_dev"
		
		data.file <- "http://galahad.well.ox.ac.uk/bigdata/ImmunoBase_GWAS.bed"
		input <- read.delim(file=data.file, header=T, stringsAsFactors=F)
		data <- paste0(input$chrom, ":", (input$chromStart+1), "-", input$chromEnd)
		
		data(ImmunoBase)
		ls_df <- lapply(ImmunoBase, function(x) as.data.frame(x$variant))
		df <- do.call(rbind, ls_df)
		data <- unique(cbind(GR=paste0(df$seqnames,":",df$start,"-",df$end), Sig=df$Pvalue))
		data <- unique(data[as.numeric(data[,2])<5e-8,1])

		ls_eTerm <- xGR2xGeneAnnoAdv(list_vec=data, crosslink=c("genehancer","PCHiC_combined","PCHiC_combined_PE","GTEx_V6p_combined","eQTL_ImmuneCells_combined","eQTL_scRNAseq_combined","FANTOM5_Cell","FANTOM5_Tissue", "REG_lncRNA","REG_enhancer", "nearby")[5], crosslink.top=2000, ontologies="AA", size.range=c(10,2000), min.overlap=5, test="fisher", RData.location=RData.location)
		gp <- xEnrichForest(ls_eTerm, top_num=10, CI.one=F)
		gp
		
		df_xGenes_GR <- xGR2xGenes(data=data, format="chr:start-end", crosslink=c("genehancer","PCHiC_combined","PCHiC_combined_PE","GTEx_V6p_combined","eQTL_ImmuneCells_combined","eQTL_scRNAseq_combined","FANTOM5_Cell","FANTOM5_Tissue", "REG_lncRNA","REG_enhancer", "nearby")[5], cdf.function="original", scoring=F, RData.location=RData.location)
		
		eTerm <- xGR2xGeneAnno(data, format="chr:start-end", crosslink=c("genehancer","PCHiC_combined","PCHiC_combined_PE","GTEx_V6p_combined","eQTL_ImmuneCells_combined","eQTL_scRNAseq_combined","FANTOM5_Cell","FANTOM5_Tissue", "REG_lncRNA","REG_enhancer", "nearby")[5], crosslink.top=2000, ontology="AA", size.range=c(10,2000), min.overlap=5, test="fisher", RData.location=RData.location)
		
		## genehancer
		df_xGenes <- xGR2xGenes(data, format="chr:start-end", crosslink="REG_enhancer", crosslink.customised=NULL, cdf.function="original", scoring=T, scoring.scheme="max", scoring.rescale=F, nearby.distance.max=50000, nearby.decay.kernel="rapid", nearby.decay.exponent=2, verbose=T, RData.location=RData.location)
		df_xGenes_GR <- xGR2xGenes(data, format="chr:start-end", crosslink="genehancer", crosslink.customised=NULL, cdf.function="original", scoring=F, nearby.distance.max=50000, nearby.decay.kernel="rapid", nearby.decay.exponent=2, verbose=T, RData.location=RData.location)

		## nearby
		df_xGenes_GR <- xGR2xGenes(data, format="chr:start-end", crosslink="nearby", crosslink.customised=NULL, cdf.function="original", scoring=T, nearby.distance.max=50000, nearby.decay.kernel="rapid", nearby.decay.exponent=2, verbose=T, RData.location=RData.location)
		
	}
	#############################
		
	# read input file
	data_input <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1]
	
	if(background.file!=""){
		# read background file
		background <- read.delim(file=background.file, header=F, stringsAsFactors=F)[,1]	
	}else{
		background <- NULL
	}
	
	distance.max <- as.numeric(distance.max)
	
	size.range <- as.numeric(c(size_range_min, size_range_max))
	min.overlap <- as.numeric(min_overlap)
	
	crosslink.top <- as.numeric(crosslink_top)
	
	# perform enrichment analysis
	#ls_eTerm <- xGR2xGeneAnnoAdv(list_vec=data_input, background=background, build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=NULL, crosslink.top=crosslink.top, nearby.distance.max=distance.max, nearby.decay.kernel=decay.kernel, nearby.decay.exponent=2, ontologies="AA", size.range=size.range, min.overlap=min.overlap, test=test, RData.location=RData.location, plot=F, ...)
	eTerm <- xGR2xGeneAnno(data_input, background=background, format="chr:start-end", build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=NULL, crosslink.top=crosslink.top, nearby.distance.max=distance.max, nearby.decay.kernel=decay.kernel, nearby.decay.exponent=2, ontology="AA", size.range=size.range, min.overlap=min.overlap, test=test, RData.location=RData.location, ...)

	#if(class(ls_eTerm)=="ls_eTerm"){
	if(class(eTerm)=="eTerm"){
		# save enrichment results to the output file
		res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
		res$id <- rownames(res)
		#res <- ls_eTerm$df
		res$CI <- paste0("[",res$CIl,",",res$CIu,"]")
		### by default
		res_f <- data.frame(term=rownames(res), res, link=rep("",nrow(res)))
		### for AA only
		if(vis!="no"){
			passflag <- input.file
			passflag <- gsub("^public/tmp/","",passflag)
			passflag <- gsub(".txt$","",passflag)
			
			##########################
			if(0){
				ls_res <- lapply(res$members_Overlap, function(x){
					unlist(strsplit(x, ", "))
				})
				vec_res <- unique(unlist(ls_res))
				y <- data.frame(gene=vec_res, color=rep(1,length(vec_res)), stringsAsFactors=F)
			}else{
				y <- eTerm$crosslink[,c("Gene","Score")]
			}
			##########################
			
			if(ncol(y)>=2){
				y[,2] <- as.numeric(y[,2])
				y <- y[!is.na(y[,2]),]
				z <- y[,2]
				if(max(z) <=0){
					colormap <- "deepskyblue-white"
				}else if(min(z) >=0){
					colormap <- "white-darkorange"
				}else{
					colormap <- "deepskyblue-white-darkorange"
				}
			}else{
				colormap <- "NULL"
			}
			
			if(0){
				## de novo graphml
				passflag <- paste0(passflag,"_",colormap)
				#####
				### passflag is: input.file "_" colormap
				#####
				link <- paste("/A2/explorer/", ls_eTerm$df$id, "/", passflag, sep="")
				link[ls_eTerm$df$adjp >= 0.05] <- ""
				res_f <- data.frame(term=rownames(res), res, link=link)
			
			}else{
				## pre-created graphml
				data <- y
				if(ncol(data)==2){
					colnames(data) <- c("label","lfc")
				}else if(ncol(data)==1){
					data <- cbind(data,rep(1,nrow(data)))
					colnames(data) <- c("label","lfc")
				}else{
					colnames(data) <- c("label","lfc","fdr")
				}
				
				if(colormap=="NULL"){
					z <- data[,2]
					if(max(z) <=0){
						colormap <- "deepskyblue-white"
					}else if(min(z) >=0){
						colormap <- "white-darkorange"
					}else{
						colormap <- "deepskyblue-white-darkorange"
					}
				}
				
				out <- passflag
				out <- gsub("data.Regions.AA.","",out)
				#####
				### out is: $rand_flag
				#####
				df_enrichment <- res
				if(vis=="yes_sig"){
					df_enrichment <- df_enrichment[df_enrichment$adjp<0.05, ]
				}
				ls_tmp <- lapply(1:nrow(df_enrichment), function(i){
					query <- df_enrichment$id[i]
					out.file <- paste0("public/",out,".",query,".graphml")
					xGraphML2AA(data=data, query=query, node.label="label", node.color="lfc", colormap=colormap, node.highlight="fdr", filename=out.file, RData.location=RData.location, verbose=F)
				})
				
				link <- paste("/A2/explorer/", df_enrichment$id, "/", out, sep="")
				if(vis=="yes_sig"){
					link[df_enrichment$adjp >= 0.05] <- ""
				}
				res_f <- data.frame(term=rownames(res), res, link=link)
			}
			
		}
		###
		utils::write.table(res_f[,c("id","name","adjp","zscore","or","CI","nOverlap","members_Overlap")], file=output.file, sep="\t", row.names=FALSE)
	
		# save to the json file
		res_f <- toJSON(res_f, pretty=T, digits=10)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file, perl=T)
		base::write(res_f, file=output.file.json)
		
		###########
		## CTree plot of enrichment results
		###########
		AA.template <- xRDataLoader("AA.template", RData.location=RData.location)
		# consensus tree
		ig <- AA.template$consensus$ig

		# output file ".pdf"
		#gp <- xEnrichCtree(ls_eTerm, ig, node.size.range=c(0.5,5.5), leave.label.expansion=1.15, leave.label.size=3, leave.label.color="steelblue", limit.expansion=2, edge.width=0.5, leave.label.orientation=c("inwards","outwards")[2])
		gp <- xEnrichCtree(eTerm, ig, node.size.range=c(0.5,5.5), leave.label.expansion=1.15, leave.label.size=3, leave.label.color="steelblue", limit.expansion=2, edge.width=0.5, leave.label.orientation=c("inwards","outwards")[2])
		if(!is.null(gp)){
			output.file.pdf <- gsub(".txt$", ".pdf", output.file, perl=T)
			width <- 6.5
			height <- 6.5*1.1
			pdf(output.file.pdf, width=width, height=height, compress=T)
			print(gp)
			dev.off()
		}
		
		# output file ".html"
		#gp_atlas <- xEnrichCtree(ls_eTerm, ig, node.size.range=c(0.5,5.5), leave.label.expansion=1.15, leave.label.size=3, leave.label.color="steelblue", limit.expansion=2, edge.width=0.5, leave.label.orientation=c("inwards","outwards")[2])
		gp_atlas <- xEnrichCtree(eTerm, ig, node.size.range=c(0.5,5.5), leave.label.expansion=1.15, leave.label.size=3, leave.label.color="steelblue", limit.expansion=2, edge.width=0.5, leave.label.orientation=c("inwards","outwards")[2])
		df <- subset(gp_atlas$data,leaf==T)
		#33897825.AA:hsa04310.graphml
		df$onclick <- sprintf("window.open(\"/A2/explorer/%s/%s\")", df$id, out)	
		#### only enriched with onclick event
		ind <- match(df$id, df_enrichment$id)
		df$onclick[is.na(ind)] <- ""
		####		
		df$tooltip <- paste0("- ", df$id, "\n- ", df$name, "\n- ", df$subcategory)
		gg <- gp_atlas + ggiraph::geom_point_interactive(data=df, aes(x=x, y=y, tooltip=tooltip, data_id=name, onclick=onclick), size=5, alpha=0.01, color="white")
		gr <- ggiraph::girafe(ggobj=gg, width_svg=6, height_svg=5, width=0.6)
		tooltip_css <- "background-color:white;color:black;font-size:10pt;font-style:normal;font-weight:bold;padding:10px;border-radius:5px;"
		hover_css <- "fill:orange;"
		gr <- ggiraph::girafe_options(gr, ggiraph::opts_tooltip(css=tooltip_css, opacity=0.75, use_fill=T), ggiraph::opts_hover(css=hover_css), ggiraph::opts_toolbar(position="topright", saveaspng=F))
		if(!is.null(gp_atlas)){
			output.file.html <- gsub(".txt$", ".html", output.file, perl=T)
			output.file.dir <- gsub(".txt$", "_files", output.file, perl=T)
			htmlwidgets::saveWidget(widget=gr, file=basename(output.file.html), background="white", selfcontained=F)
			if(base::file.exists(basename(output.file.html))){
				base::file.copy(from=basename(output.file.html), to=dirname(output.file.html), overwrite=T)
				file.remove(basename(output.file.html))
				
				base::file.copy(from=basename(output.file.dir), to=dirname(output.file.dir), overwrite=T, recursive=T)
				
			}

		}
		
		##########################################
		# output file ".xGenes.txt"
		## forest plot of enrichment results
		if(output_xGenes!="no"){
			
			# pathway genes
			ls_res <- lapply(res$members_Overlap, function(x){
				unlist(strsplit(x, ", "))
			})
			vec_genes <- unique(unlist(ls_res))
			
			if(output_xGenes=="yes_xGenes_GR" | output_xGenes=="yes_xGenes"){
				output.file.txt <- gsub(".txt$", ".xGenes.txt", output.file, perl=T)
				if(0){
					df_xGenes <- xGR2xGenes(data=data_input, format="chr:start-end", build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=NULL, cdf.function="original", scoring=T, scoring.scheme="max", scoring.rescale=F, nearby.distance.max=distance.max, nearby.decay.kernel=decay.kernel, nearby.decay.exponent=2, RData.location=RData.location)
				}else{
					# restricted by crosslink.top
					df_xGenes <- eTerm$crosslink
				}
				df_xGenes <- df_xGenes[, c("Gene","Score")]
				if(0){
					ind <- match(df_xGenes$Gene, vec_genes)
					df_xGenes <- df_xGenes[!is.na(ind),]
				}
				utils::write.table(df_xGenes, file=output.file.txt, sep="\t", row.names=FALSE, quote=F)
			}
			
			if(output_xGenes=="yes_xGenes_GR"){
				output.file.txt <- gsub(".txt$", ".xGenes_GR.txt", output.file, perl=T)
				if(0){
					df_xGenes_GR <- xGR2xGenes(data=data_input, format="chr:start-end", build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=NULL, cdf.function="original", scoring=F, nearby.distance.max=distance.max, nearby.decay.kernel=decay.kernel, nearby.decay.exponent=2, RData.location=RData.location)
				}else{
					# restricted by crosslink.top
					df_xGenes_GR <- eTerm$evidence
				}
				df_evidence <- df_xGenes_GR[, c("GR","Gene","Score")]
				if(0){
					ind <- match(df_evidence$Gene, vec_genes)
					df_evidence <- df_evidence[!is.na(ind),]
				}
				utils::write.table(df_evidence, file=output.file.txt, sep="\t", row.names=FALSE, quote=F)

				if(1){
					# also restricted by crosslinked pathway genes
					ind <- match(df_evidence$Gene, vec_genes)
					df_evidence <- df_evidence[!is.na(ind),]
				}
				output.file.pdf <- gsub(".txt$", ".xGenes_GR.pdf", output.file, perl=T)
				Gene <- Score <- NULL
				mat_evidence <- tidyr::spread(df_evidence, key=Gene, value=Score)
				mat <- mat_evidence[,-1]
				rownames(mat) <- mat_evidence[,1]
				#### sort by chromosome, start and end
				ind <- xGRsort(rownames(mat))
				mat <- mat[ind,]
				
				x <- gsub(":.*|chr","",rownames(mat))
				x[x=="X"] <- 23
				x[x=="Y"] <- 24
				x <- as.numeric(x)
				rowsep <- cumsum(table(x))
				rowsep <- nrow(mat) - rowsep
				rowsep <- rowsep[-length(rowsep)]

				####
				if(ncol(mat)>=0){
					reorder <- "none"
				}else{
					reorder <- "col"
				}
				gp_evidence <- xHeatmap(mat, reorder=reorder, colormap="spectral", ncolors=64, barwidth=0.4, x.rotate=90, shape=19, size=2, x.text.size=6,y.text.size=6, na.color="transparent")
				gp_evidence <- gp_evidence + theme(legend.title=element_text(size=8), legend.position="left") + scale_y_discrete(position="right")
				gp_evidence <- gp_evidence + geom_hline(yintercept=rowsep+0.5,color="grey90",size=0.5)
				width <- 2 + ncol(mat)*0.075
				height <- 2 + nrow(mat)*0.075
				pdf(output.file.pdf, width=width, height=height, compress=T)
				print(gp_evidence)
				dev.off()
				
			}
			
		}
		##########################################

	}
}
';

# for calling R function
$my_rscript.="
library(XGR)
library(jsonlite)
R_pipeline(input.file=\"$input_filename\", background.file=\"$background_filename\", output.file=\"$output_filename\", build.conversion=\"$build_conversion\", crosslink=\"$crosslink\", crosslink_top=\"$crosslink_top\", distance.max=\"$distance_max\", decay.kernel=\"$decay_kernel\", ontology=\"$ontology\", size_range_min=\"$size_range_min\", size_range_max=\"$size_range_max\", min_overlap=\"$min_overlap\", test=\"$test\", vis=\"$vis\", output_xGenes=\"$output_xGenes\", RData.location=\"$placeholder\")
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
			### for json
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.json/g;
			if(-e $tmp_file){
				$ajax_json_file=$tmp_file;
				$ajax_json_file=~s/^public//g;
				print STDERR "JSON locates at $ajax_json_file\n";
			}
			
			#############
			### for pdf
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.pdf/g;
			if(-e $tmp_file){
				$ajax_pdf_file=$tmp_file;
				$ajax_pdf_file=~s/^public//g;
				print STDERR "PDF locates at $ajax_pdf_file\n";
			}
			#### .xGenes_GR.pdf
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.xGenes_GR.pdf/g;
			if(-e $tmp_file){
				$ajax_pdf_file_xGenes_GR=$tmp_file;
				$ajax_pdf_file_xGenes_GR=~s/^public//g;
				print STDERR "TXT locates at $ajax_pdf_file_xGenes_GR\n";
			}
			
			#############
			### for txt
			#### .xGenes.txt
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.xGenes.txt/g;
			if(-e $tmp_file){
				$ajax_txt_file_xGenes=$tmp_file;
				$ajax_txt_file_xGenes=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file_xGenes\n";
			}
			#### .xGenes_GR.txt
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.xGenes_GR.txt/g;
			if(-e $tmp_file){
				$ajax_txt_file_xGenes_GR=$tmp_file;
				$ajax_txt_file_xGenes_GR=~s/^public//g;
				print STDERR "TXT locates at $ajax_txt_file_xGenes_GR\n";
			}
			
			#############
			### for html
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.html/g;
			if(-e $tmp_file){
				$ajax_html_file=$tmp_file;
				$ajax_html_file=~s/^public//g;
				print STDERR "HTML locates at $ajax_html_file\n";
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
	
	# stash $ajax_json_file;
	$c->stash(ajax_json_file => $ajax_json_file);
	
	# stash $txtfile;
	my $txtfile;
	$txtfile->{xGenes}=$ajax_txt_file_xGenes;
	$txtfile->{xGenes_GR}=$ajax_txt_file_xGenes_GR;
	$c->stash(txtfile => $txtfile);

	# stash $pdffile;
	my $pdffile;
	$pdffile->{ctree}=$ajax_pdf_file;
	$pdffile->{xGenes_GR}=$ajax_pdf_file_xGenes_GR;
	$c->stash(pdffile => $pdffile);
	
	# stash $htmlfile;
	my $htmlfile;
	$c->stash(htmlfile => $ajax_html_file);
	
	my $num_genes=0;
	if(defined($genelist)){
		my @lines = split(/\r\n|\n/, $genelist);
		my %rec_genes;
		foreach my $line (@lines) {
			next if($line=~/^\s*$/);
			my $rec;
			$rec->{ID}=$line;
			$rec->{Name}=$line;
			$rec_genes{$line}=$rec;
		}  	
		$c->stash(rec_genes => \%rec_genes);
		$num_genes=scalar(keys %rec_genes);
	}
	
	$c->stash(num_genes => $num_genes);
	$c->stash(ontology => $ontology);
	
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



# Render template "A2Crosstalk.html.ep"
sub A2_crosstalk {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $ontology = $c->param('ontology') || 'DO'; # by default: DO
  	my $genelist = $c->param('genelist');
  	my $uploadfile = $c->req->upload('uploadfile');
	if(defined($uploadfile)) {
		$uploadfile = $uploadfile->slurp;
	}else{
		$uploadfile='';
	}
	
	my $size_range_min = $c->param('size_range_min') || 1;
	my $size_range_max = $c->param('size_range_max') || 1000000;
	my $min_overlap = $c->param('min_overlap') || 1;
	
	my $test = $c->param('test');
	my $vis = $c->param('vis');
  	
  	# The output json file (default: '')
	my $ajax_json_file='';
  	
	# The output pdf file (default: '')
	my $ajax_pdf_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $PIER_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $ontology.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'enrichment.Genes.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'enrichment.Genes.'.$rand_file.'.r';
	
		my $my_input;
		if(0){
		foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $genelist)) {
			next if($line=~/^\s*$/);
			$my_input.=$line."\n";
		}
		}else{
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
		}
		}
		PIER_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $background_filename="";
		if($uploadfile ne ''){
			my $my_background;
			foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $uploadfile)) {
				next if($line=~/^\s*$/);
				$my_background.=$line."\n";
			}
			$background_filename=$tmpFolder.'/'.'background.'.$rand_file.'.txt';
			PIER_app::Controller::Utils::export_to_file($background_filename, $my_background);
		}
		
		my $placeholder;
		if(-e '/Users/hfang/Sites/SVN/github/bigdata_dev'){
			$placeholder="/Users/hfang/Sites/SVN/github/bigdata_dev";
		}else{
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
R_pipeline <- function (input.file="", background.file="", output.file="", ontology="", size_range_min="", size_range_max="", min_overlap="", test="", vis="", RData.location="", ...){
	
	#############################
	## an example for coding test
	if(0){
		RData.location <- "/Users/hfang/Sites/SVN/github/bigdata_dev"
		
		input.file <- "/Users/hfang/Sites/XGR/A2-site/a2_app/public/app/examples/PDL1.2columns.FDR.txt"
		data <- read.delim(file=input.file, header=F, stringsAsFactors=F)
		
		cPath <- xCrosstalk(data, entity=c("Gene","GR")[1], significance.threshold=NULL, score.cap=NULL, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("genehancer","nearby"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring.scheme=c("max","sum","sequential"), nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, networks=c("KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease")[c(2,3,4,5,6,7)], seed.genes=T, subnet.significance=0.05, subnet.size=NULL, ontologies=c("KEGGenvironmental","KEGG","KEGGmetabolism","KEGGgenetic","KEGGcellular","KEGGorganismal","KEGGdisease", "REACTOME", "REACTOME_ImmuneSystem", "REACTOME_SignalTransduction", "AA")[c(1,3,4,5)], size.range=c(10,2000), min.overlap=10, fdr.cutoff=0.05, glayout=layout_with_kk, verbose=T, RData.location=RData.location)
		
		
		g <- xRDataLoader(RData.customised="ig.PCHiC", RData.location=RData.location, verbose=verbose)
		nodes <- igraph::get.data.frame(g, what="vertices")[,"name"]
		df_cGenes <- xSNP2cGenes(data=nodes, entity="chr:start-end", include.HiC="Combined", cdf.function="empirical", verbose=TRUE, RData.location=RData.location)
		crosslink.customised.cGenes <- unique(df_cGenes[,c("SNP","Gene","Sig")])
		colnames(crosslink.customised.cGenes) <- c("GR","Gene","Score")
		crosslink.customised.cGenes$Context <- "PCHiC_combined"
		
		
		
		pdf("a.pdf", width=7, height=8)
		gp_both <- gridExtra::grid.arrange(grobs=list(cPath$gp_paths,cPath$gp_heatmap), layout_matrix=cbind(c(1,1,1,1,2)))
		dev.off()
		
		df <- cPath$g$enrichment
		df$nOverlap <- df$nPath
		df$members <- df$path_members
		gp_forest <- xEnrichForest(df, top_num=5, FDR.cutoff=0.05, CI.one=F, signature=F)
		gp_treemap <- xEnrichTreemap(df, top_num=5, FDR.cutoff=0.05, CI.one=F, details=c("name","name_FDR","name_FDR_members")[3], caption=F)
		
		pdf("a.pdf", width=7, height=8)
		#gp_both <- gridExtra::grid.arrange(grobs=list(cPath$gp_paths,gp_forest), layout_matrix=cbind(c(1,1,1,1,2)))
		gp_both <- gridExtra::grid.arrange(grobs=list(cPath$gp_paths,gp_treemap), layout_matrix=cbind(c(1,1,1,1,2)))
		dev.off()
		
		
		
		
		g_tmp <- igraph::disjoint_union(ls_subg)
		if(0){
			set.seed(825)
			layouts <- lapply(ls_subg, layout_with_kk)
			glayout <- igraph::merge_coords(ls_subg, layouts)
			node.xcoord <- glayout[,1]
			node.ycoord <- glayout[,2]
			## scale into [-1,1]
			if(max(node.xcoord) != min(node.xcoord)){
				node.xcoord <- (node.xcoord - min(node.xcoord)) / (max(node.xcoord) - min(node.xcoord)) * 2 - 1
			}
			if(max(node.ycoord) != min(node.ycoord)){
				node.ycoord <- (node.ycoord - min(node.ycoord)) / (max(node.ycoord) - min(node.ycoord)) * 2 - 1
			}
			glayout <- cbind(node.xcoord, node.ycoord)
			V(g_tmp)$xcoord <- glayout[,1]
			V(g_tmp)$ycoord <- glayout[,2]
		}
		gp <- xGGnetwork(g=g_tmp, node.label="name", node.label.size=2, node.label.color="black", node.label.alpha=0.8, node.label.padding=0.1, node.label.arrow=0, node.label.force=0.01, node.shape=19, node.xcoord="xcoord", node.ycoord="ycoord", node.color="color", node.color.title=expression(-log[10]("input")), colormap="jet.top", ncolors=64, node.size.range=5, edge.color="orange",edge.color.alpha=0.3,edge.curve=0,edge.arrow.gap=0.02, title=paste0("Pathway crosstalk involving ",vcount(g_tmp)," genes"), zlim=NULL)
		gp
		##################################
		
		
		
	}
	#############################
		
	# read input file
	data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1]
	
	if(background.file!=""){
		# read background file
		background <- read.delim(file=background.file, header=F, stringsAsFactors=F)[,1]	
	}else{
		background <- NULL
	}
	
	size.range <- as.numeric(c(size_range_min, size_range_max))
	min.overlap <- as.numeric(min_overlap)
	
	# perform enrichment analysis
	ls_eTerm <- xA2EnricherGenes(list_vec=data, background=background, ontologies=ontology, size.range=size.range, min.overlap=min.overlap, test=test, RData.location=RData.location, plot=F, ...)
	
	if(class(ls_eTerm)=="ls_eTerm"){
		# save enrichment results to the output file
		#res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
		res <- ls_eTerm$df
		res$CI <- paste0("[",res$CIl,",",res$CIu,"]")
		### by default
		res_f <- data.frame(term=rownames(res), res, link=rep("",nrow(res)))
		### for AA only
		if(ontology == "AA" & vis!="no"){
			passflag <- input.file
			passflag <- gsub("^public/tmp/","",passflag)
			passflag <- gsub(".txt$","",passflag)
			
			y <- read.delim(file=input.file, header=F, stringsAsFactors=F)
			if(ncol(y)>=2){
				y[,2] <- as.numeric(y[,2])
				y <- y[!is.na(y[,2]),]
				z <- y[,2]
				if(max(z) <=0){
					colormap <- "deepskyblue-white"
				}else if(min(z) >=0){
					colormap <- "white-darkorange"
				}else{
					colormap <- "deepskyblue-white-darkorange"
				}
			}else{
				colormap <- "NULL"
			}
			
			if(0){
				## de novo graphml
				passflag <- paste0(passflag,"_",colormap)
				#####
				### passflag is: input.file "_" colormap
				#####
				link <- paste("/A2/explorer/", ls_eTerm$df$id, "/", passflag, sep="")
				link[ls_eTerm$df$adjp >= 0.05] <- ""
				res_f <- data.frame(term=rownames(res), res, link=link)
			
			}else{
				## pre-created graphml
				data <- y
				if(ncol(data)==2){
					colnames(data) <- c("label","lfc")
				}else if(ncol(data)==1){
					data <- cbind(data,rep(1,nrow(data)))
					colnames(data) <- c("label","lfc")
				}else{
					colnames(data) <- c("label","lfc","fdr")
				}
				
				if(colormap=="NULL"){
					z <- data[,2]
					if(max(z) <=0){
						colormap <- "deepskyblue-white"
					}else if(min(z) >=0){
						colormap <- "white-darkorange"
					}else{
						colormap <- "deepskyblue-white-darkorange"
					}
				}
				
				out <- passflag
				out <- gsub("data.Genes.AA.","",out)
				#####
				### out is: $rand_flag
				#####
				df <- ls_eTerm$df
				if(vis=="yes_sig"){
					df <- df[df$adjp<0.05, ]
				}
				ls_tmp <- lapply(1:nrow(df), function(i){
					query <- df$id[i]
					out.file <- paste0("public/",out,".",query,".graphml")
					xGraphML2AA(data=data, query=query, node.label="label", node.color="lfc", colormap=colormap, node.highlight="fdr", filename=out.file, RData.location=RData.location, verbose=F)
				})
				
				link <- paste("/A2/explorer/", ls_eTerm$df$id, "/", out, sep="")
				if(vis=="yes_sig"){
					link[ls_eTerm$df$adjp >= 0.05] <- ""
				}
				res_f <- data.frame(term=rownames(res), res, link=link)
			}
			
		}
		###
		utils::write.table(res_f, file=output.file, sep="\t", row.names=FALSE)
	
		# save to the json file
		res_f <- toJSON(res_f, pretty=T, digits=10)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file, perl=T)
		base::write(res_f, file=output.file.json)
		
		# output file ".pdf"
		## forest plot of enrichment results
		gp <- xEnrichForest(ls_eTerm, top_num=10, FDR.cutoff=0.05, CI.one=T, signature=F)
		if(!is.null(gp)){
			output.file.pdf <- gsub(".txt$", ".pdf", output.file, perl=T)
			width <- 7
			height <- 2 + nrow(gp$data)*0.1
			pdf(output.file.pdf, width=width, height=height, compress=T)
			print(gp)
			dev.off()
		}

	}
}
';

# for calling R function
$my_rscript.="
library(XGR)
library(jsonlite)
R_pipeline(input.file=\"$input_filename\", background.file=\"$background_filename\", output.file=\"$output_filename\", ontology=\"$ontology\", size_range_min=\"$size_range_min\", size_range_max=\"$size_range_max\", min_overlap=\"$min_overlap\", test=\"$test\", vis=\"$vis\", RData.location=\"$placeholder\")
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
			$tmp_file=~s/.txt$/.json/g;
			if(-e $tmp_file){
				$ajax_json_file=$tmp_file;
				$ajax_json_file=~s/^public//g;
				print STDERR "JSON locates at $ajax_json_file\n";
			}
			
			#############
			### for pdf
			$tmp_file=$output_filename;
			$tmp_file=~s/.txt$/.pdf/g;
			if(-e $tmp_file){
				$ajax_pdf_file=$tmp_file;
				$ajax_pdf_file=~s/^public//g;
				print STDERR "PDF locates at $ajax_pdf_file\n";
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
	
	# stash $ajax_json_file;
	$c->stash(ajax_json_file => $ajax_json_file);
	
	# stash $pdffile;
	my $pdffile;
	$pdffile->{forest}=$ajax_pdf_file;
	$c->stash(pdffile => $pdffile);
	
	my $num_genes=0;
	if(defined($genelist)){
		my @lines = split(/\r\n|\n/, $genelist);
		my %rec_genes;
		foreach my $line (@lines) {
			next if($line=~/^\s*$/);
			my $rec;
			$rec->{ID}=$line;
			$rec->{Name}=$line;
			$rec_genes{$line}=$rec;
		}
		$c->stash(rec_genes => \%rec_genes);
		$num_genes=scalar(keys %rec_genes);
	}
	
	$c->stash(num_genes => $num_genes);
	$c->stash(ontology => $ontology);
	
  	$c->render();

}

1;
