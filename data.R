library(readxl)
library(readr)
library(purrr)
library(stringr)
MS110 <- read_excel('../20190725 YLS SILAC/20190725YLS proteinGroups_110.xlsx')
MS211 <- read_excel('../20190725 YLS SILAC/20190725 YLSproteinGroups_211.xlsx')
MS312 <- read_excel('../20190725 YLS SILAC/20190725 YLS proteinGroups_312.xlsx')
MS74 <- read_excel('../20190725 YLS SILAC/20190725YLS proteinGroups 74.xlsx')
MS85 <- read_excel('../20190725 YLS SILAC/20190725YLS proteinGroup 85.xlsx')
MS96 <- read_excel('../20190725 YLS SILAC/20190725YLS ProteinGroup_96.xlsx')

extract.protein.name<- function(x){
  p <- parent.frame()
  if(is.data.frame(x)){
    name <- deparse(substitute(x))
    x$`Protein IDs` <- str_extract(x$`Protein IDs`, pattern = '(?<=\\|)\\w+(?=\\_HUMAN)')
    if(exists(name, envir = p, inherits = F)){
      assign(name, x, envir = p)
    } else {
      assign(name, x, envir = where(name))
      }
  }else if(is.character(x)){
    name <- x
    obj <- get(x, envir = p)
    obj$`Protein IDs` <- str_extract(obj$`Protein IDs`, pattern = '(?<=\\|)\\w+(?=\\_HUMAN)')
    if (exists(name, envir = p, inherits = F)) {
      assign(name, obj, envir = p)
    }else{
      assign(name, obj, envir = pryr::where(name))
    }
  } else {
    errorCondition("The type of argument should be 'character' or 'dataframe'")
  }
  invisible()
}
data.obj <- ls(pattern = 'MS*')
extract.protein.name(MS110)



