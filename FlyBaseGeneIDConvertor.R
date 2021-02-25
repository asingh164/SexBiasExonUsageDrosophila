# Packages required
require(tidyverse)    # General-purpose data wrangling
require(rvest)        # Parsing of HTML/XML files
require(stringr)      # String manipulation
require(rebus)        # Verbose regular expressions
require(lubridate)    # Eases DateTime manipulation
require(ratelimitr)   # Limits the rate at which a function performs tasks (to prevent my script from being 404'd)

### In house functions ###
# This function will make a URL for every FBgn gene id that you pass to R
FBgnURL <- function(FBgnIDs){
  url =  paste(as.character("http://flybase.org/reports/"), as.character(FBgnIDs$V1), sep = "")
  return(url)
}

# This function will pull the CG symbol for a user submitted list of FBgn symbols
CG_IDs <- function(FBgnIDs){
  #Sys.sleep(rate)
  apply(FBgnIDs, 1, function(x) (read_html(x[2]) %>%
                                   html_nodes(".col-sm-3") %>%
                                   html_text())[8])
}

# Script to Run ## This script is purposfully slowed down to prevent FlyBase from time me out.
FlyBaseIDs = read.delim("/plas1/amardeep.singh/RNA.Seq.Data/HomologyData/list.of.genes.unique.txt", header = FALSE, sep = "\t") # Load in list of FlyBase gene IDs
FlyBaseIDs$urls = FBgnURL(FlyBaseIDs)   # Add a url for each FlyBase gene ID
FlyBaseIDs$CG = NA
for (i in 1:nrow(FlyBaseIDs)){
  print(i)
  FlyBaseIDs$CG[i] = CG_IDs(FlyBaseIDs[i,])
  date_time<-Sys.time()
  #print(date_time)
  #print((as.numeric(Sys.time()) - as.numeric(date_time)))
  while((as.numeric(Sys.time()) - as.numeric(date_time))<(1/2)){}
}

CG.output = CG_IDs(FlyBaseIDs, 6)




## Unused Code -- Would have been a more elegent function, but I couldn't get it work so I turned to 'apply'
#FBgnToCG <- function(FBgnIDs){
#  url = as.character(FBgnIDs$urls)
#  #url = "http://flybase.org/reports/FBgn0031081"
#  #url.html = read_html(url)

#  (read_html(url) %>%
#      html_nodes(".col-sm-3") %>%
#      html_text())[8]
#}
#FlyBaseIDs$CG = FBgnToCG(FlyBaseIDs[1,])
