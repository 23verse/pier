# [The source codes behind the PiER](https://github.com/23verse/pier)

## @ Overview

> The [PiER](http://www.genetictargets.com) is web-based facilities that support ab initio and real-time genetic target prioritisation through integrative use of human disease genetics, functional genomics and protein interactions. By design, the PiER features two facilities: `elementary` and `combinatory`.

> The `elementary` facility is designed to perform specific tasks, including three online tools: [eV2CG](http://www.genetictargets.com/e/V2CG), utilising functional genomics to link disease-associated variants (particularly at the non-coding genome) to core genes likely responsible for genetic associations; [eCG2PG](http://www.genetictargets.com/e/CG2PG), using knowledge of protein interactions to ‘network’ core genes and additional peripheral genes, producing a ranked list of core and peripheral genes; and [eCrosstalk](http://www.genetictargets.com/e/Crosstalk), exploiting the information of pathway-derived interactions to identify highly-ranked genes mediating the crosstalk between molecular pathways. Each of elementary tasks giving results is sequentially piped to the next one. 

> By chaining together elementary tasks, the `combinatory` facility automates genetics-led and network-based integrative prioritisation for genetic targets at the gene level ([cTGene](http://www.genetictargets.com/c/TGene)) and at the crosstalk level ([cTCrosstalk](http://www.genetictargets.com/c/TCrosstalk)). 

>  A tutorial-like booklet, made available [here](http://www.genetictargets.com/PiERbooklet/index.html), describes step-by-step instructions on how to use.

## @ Development

> The PiER was developed using a next-generation Perl web framework [Mojolicious](https://www.mojolicious.org).

> The PiER was also built using [Bootstrap](https://getbootstrap.com), supporting the mobile-first and responsive webserver in all major platform browsers.

## @ Installation

Assume you have a `ROOT (sudo)` privilege on Ubuntu

### 1. Install Mojolicious and other perl modules

>
```ruby
sudo su
# here enter your password
curl -L cpanmin.us | perl - Mojolicious
perl -e "use Mojolicious::Plugin::PODRenderer"
perl -MCPAN -e "install Mojolicious::Plugin::PODRenderer"
perl -MCPAN -e "install DBI"
perl -MCPAN -e "install Mojo::DOM"
perl -MCPAN -e "install Mojo::Base"
perl -MCPAN -e "install LWP::Simple"
perl -MCPAN -e "install JSON::Parse"
perl -MCPAN -e "install local::lib"
perl -MCPAN -Mlocal::lib -e 'install JSON::Parse'
```


### 2. Install R and bioconductor packages

>
```ruby
sudo su
# here enter your password
wget http://www.stats.bris.ac.uk/R/src/base/R-4/R-4.1.3.tar.gz
tar xvfz R-4.1.3.tar.gz
cd ~/R-4.1.3
./configure
make
make check
make install
R # start R

install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c('Pi','rmarkdow','bookdown'), dependencies=T)
library(Pi)
```

### 3. Install pandoc

>
```ruby
sudo su
# here enter your password
wget https://github.com/jgm/pandoc/releases/download/2.17.1.1/pandoc-2.17.1.1-linux-amd64.tar.gz
tar xvzf pandoc-2.17.1.1-linux-amd64.tar.gz --strip-components 1 -C /usr/local/

R
rmarkdown::pandoc_available()
Sys.setenv(RSTUDIO_PANDOC="/usr/local/bin")
rmarkdown::render(YOUR_RMD_FILE, bookdown::html_document2(number_sections=F, theme="readable", hightlight="default"))
```


## @ Deployment

> morbo -l 'http://*:80/' script/pier_app


## @ Contact

> Please drop [email](mailto:fh12355@rjh.com.cn) for bug reports or enquiries.


