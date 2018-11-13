# Diversity maps using biomod2

### objective
# Apply species distribution models to a set of species and and create a set
# of projections under a set of climate change scenarios

## install and load libraries
if(!require(rgbif)){ ## automatically install rgbif package if needed
  install.packages("rgbif")
  require(rgbif)
}

## get the species list belonging to the genus Swartzia
spp_swartzia <- name_suggest(q='Swartzia', rank='species', limit = 100)

## clean up the species list
(spp_swartzia <- spp_swartzia[ grepl("^Swartzia ", spp_swartzia$canonicalName),])

## get species occurrences
occ_swartzia <- occ_search(taxonKey = spp_swartzia$key, continent='south_america',
                            fields = c('name','key','country','decimalLatitude','decimalLongitude'),
                            hasCoordinate=T , limit=500, return = 'data')
 
## join all data in a single data.frame
data <- NULL
for(sp in names(occ_swartzia)){
   if( !is.null( dim( occ_swartzia[[sp]] ) ) ){
     data <- rbind(data, occ_swartzia[[sp]]) 
   }
}

## replace " " by "." in species names
data$name <- sub(" ", ".", data$name)

## keep only species having more than 200 occurrences
table(data$name)

(spp_to_model <- names(table(data$name))[ table(data$name)>199 ] )

##delete spurious coordinates (e.g., X = 0 ou Y = 0)
ids <- which(data$decimalLongitude==0 | data$decimalLatitude==0)

data <- data[-ids,]

##plot and check remaining points
# define the extent; e.g. South America (assuming there is a South America
# shapefile in the working directory)
sa <- shapefile("South_America")

plot(data[, 3:4])

ids <- which(data$decimalLongitude < -80 | data$decimalLongitude > -20)

data <- data[-ids,]

plot(data[, 3:4])
plot(sa, add = T)

## install and load libraries
if(!require(biomod2)){ ## automatically install biomod2 package if needed
  install.packages("biomod2")
  require(biomod2)
}

if(!require(gridExtra)){ ## automatically install gridExtra package if needed
  install.packages("gridExtra")
  require(gridExtra)
}

if(!require(rasterVis)){ ## automatically install rasterVis package if needed
  install.packages("rasterVis")
  require(rasterVis)
}

if(!require(rgdal)){ ## automatically install rgdal package if needed
  install.packages("rgdal")
  require(rgdal)
}

if(!require(ade4)){ ## automatically install ade4 package if needed
  install.packages("ade4")
  require(ade4)
}

## get some worldclim environmental variables
dir.create("WorldClim_data", showWarnings = F)
## curent bioclim
download.file(url = "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_esri.zip", 
              destfile = "WorldClim_data/current_bioclim_10min.zip", 
              method = "auto")

## GCM -> BCC-CSM1-1, year -> 2050, RCP -> 4.5
download.file(url = "http://biogeo.ucdavis.edu/data/climate/cmip5/10m/bc45bi50.zip", 
              destfile = "WorldClim_data/2050_BC_45_bioclim_10min.zip", 
              method = "auto")

## GCM -> BCC-CSM1-1, year -> 2080, RCP -> 4.5
download.file(url = "http://biogeo.ucdavis.edu/data/climate/cmip5/10m/bc45bi70.zip", 
              destfile = "WorldClim_data/2070_BC_45_bioclim_10min.zip", 
              method = "auto")

## unzip climatic files
unzip( zipfile = "WorldClim_data/current_bioclim_10min.zip", 
       exdir = "WorldClim_data/current", 
       overwrite = T)
list.files("WorldClim_data/current/bio/")

unzip( zipfile = "WorldClim_data/2050_BC_45_bioclim_10min.zip", 
       exdir = "WorldClim_data/2050/BC_45",
       overwrite = T)
list.files("WorldClim_data/2050/BC_45/")

unzip( zipfile = "WorldClim_data/2070_BC_45_bioclim_10min.zip",
       exdir = "WorldClim_data/2070/BC_45",
       overwrite = T)
list.files("WorldClim_data/2070/BC_45/")

## load environmental variables within a 'RasterStack' object
stk_current <- stack( list.files(path = "WorldClim_data/current/bio/", 
                                  pattern = "bio_", 
                                  full.names = T), RAT=F )

## Mask and crop environmental variables to South Americas's extent
bioclim_sa <- mask(stk_current, sa)

bioclim_sa <- crop(bioclim_sa, sa)

## convert our environmental Stack into data.frame
current_df <- as.data.frame(bioclim_sa)
current_df <- na.omit(current_df)
 
## calculate Pearson correlations between pairs of variables 
cor_current <- cor(current_df)
 
## reformat correlation table for graphical analyse
cor_current[ upper.tri(cor_current, diag = T) ] <- NA
cor_current_resh <- na.omit( melt( cor_current ) )
colnames(cor_current_resh) <- c("var1", "var2", "correlation")
 
## only consider absolute value of correlations
cor_current_resh$correlation <- abs(cor_current_resh$correlation)
 
## make a correlation plot and select variables
gg_cor <- ggplot(cor_current_resh, aes(x = var1, y = var2 , fill = correlation) )
gg_cor <- gg_cor + geom_tile() + xlab("") + ylab("") + theme( axis.text.x  = element_text(angle=90, vjust=0.5))

print(gg_cor)

selected_vars <- c("bio_9", "bio_15", "bio_18")

## check correlations between selected variables
(cor_sel <- cor(current_df[,selected_vars]))

## keep only the selected variables (It's advisable to discuss the ecological
## relevance of the climatic variables with specialists on the taxonomic
## group under study. Also, try to keep an average of one variable per five
## occurrences to avoid overparameterization;
stk_current <- stack(subset(bioclim_sa, selected_vars))

## do the same for future scenarios
## NOTE : keep layer names and order across stacks
stk_2050_BC_45 <- stack( c( bio_9 = "WorldClim_data/2050/BC_45/bc45bi509.tif",
                             bio_15 = "WorldClim_data/2050/BC_45/bc45bi5015.tif",
                             bio_18 = "WorldClim_data/2050/BC_45/bc45bi5018.tif"),
                          RAT=F)
stk_2050_BC_45 <- stack(mask(stk_2050_BC_45, sa))
stk_2050_BC_45 <- stack(crop(stk_2050_BC_45, sa))
 
stk_2070_BC_45 <- stack( c( bio_9 = "WorldClim_data/2070/BC_45/bc45bi709.tif",
                             bio_15 = "WorldClim_data/2070/BC_45/bc45bi7015.tif",
                             bio_18 = "WorldClim_data/2070/BC_45/bc45bi7018.tif"),
                          RAT=F)
stk_2070_BC_45 <- stack(mask(stk_2070_BC_45, sa))
stk_2070_BC_45 <- stack(crop(stk_2070_BC_45, sa))

##wrapper 
biomod2_wrapper <- function(sp){
   cat("\n> species : ", sp)
   
   ## get occurrences points
   sp_dat <- data[ data$name == sp, ]
   
   ## formating the data
   sp_format <- BIOMOD_FormatingData(resp.var = rep( 1, nrow(sp_dat) ), 
                                     expl.var = stk_current,
                                     resp.xy = sp_dat[,c("decimalLongitude","decimalLatitude")],
                                     resp.name = sp,
                                     PA.strategy = "random", 
                                     PA.nb.rep = 3, 
                                     PA.nb.absences = 1000)
   ## print formatting summary
   sp_format
   
   ## save image of input data summary
   if(!exists(sp)) dir.create(sp)
   pdf(paste(sp, "/", sp ,"_data_formated.pdf", sep="" ))
   try(plot(sp_format))
   dev.off()
   
   ## define models options
   sp_opt <- BIOMOD_ModelingOptions()
   
   ## model species (see alternative models and models.eval.meth in
   ## 'help(BIOMOD_Modeling)'. If applying MAXENT, maxent.jar must be in the
   ## working directory)
   sp_model <- BIOMOD_Modeling( sp_format, 
                                models = c('GLM','FDA','RF'), 
                                models.options = sp_opt, 
                                NbRunEval = 3, 
                                DataSplit = 70, 
                                Yweights = NULL, 
                                VarImport = 3, 
                                models.eval.meth = c('TSS','ROC'),
                                SaveObj = TRUE,
                                rescal.all.models = FALSE,
                                do.full.models = FALSE,
                                modeling.id = "ex3")
   
   ## save some graphical outputs
   #### models scores
   pdf(paste(sp, "/", sp ,"_models_scores.pdf", sep="" ))
   try( gg1 <- models_scores_graph(sp_model, metrics = c("TSS","ROC"), by = 'models', plot=F) )
   try( gg2 <- models_scores_graph(sp_model, metrics = c("TSS","ROC"), by = 'data_set', plot=F) )
   try( gg3 <- models_scores_graph(sp_model, metrics = c("TSS","ROC"), by = 'cv_run', plot=F) )
   try(grid.arrange(gg1,gg2,gg3))
   dev.off()
   
   ## build ensemble models
   sp_ens_model <- BIOMOD_EnsembleModeling( modeling.output = sp_model,
                                            chosen.models = 'all',
                                            em.by = 'all',
                                            eval.metric = c('TSS'),
                                            eval.metric.quality.threshold = c(0.7),
                                            models.eval.meth = c('TSS','ROC'),
                                            prob.mean = TRUE,
                                            prob.cv = TRUE,
                                            prob.ci = FALSE,
                                            prob.ci.alpha = 0.05,
                                            prob.median = FALSE,
                                            committee.averaging = TRUE,
                                            prob.mean.weight = TRUE,
                                            prob.mean.weight.decay = 'proportional' )
   
   ## do projections
   proj_scen <- c("current", "2050_BC_45", "2070_BC_45")
   
   for(scen in proj_scen){
     cat("\n> projections of ", scen)
     
     ## single model projections
     sp_proj <- BIOMOD_Projection(  modeling.output = sp_model,
                                    new.env = get(paste("stk_", scen, sep = "")),
                                    proj.name = scen,
                                    selected.models = 'all',
                                    binary.meth = "TSS",
                                    filtered.meth = NULL,
                                    compress = TRUE,
                                    build.clamping.mask = TRUE,
                                    do.stack = FALSE,
                                    output.format = ".img" )
     
     ## ensemble model projections
     sp_ens_proj <- BIOMOD_EnsembleForecasting(EM.output = sp_ens_model,
                                               projection.output = sp_proj,
                                               binary.meth = "TSS",
                                               compress = TRUE,
                                               do.stack = FALSE,
                                               output.format = ".img")
     
   }
   
   return(paste(sp," modelling completed !", sep=""))
   
 }


##mod function
all_species_bm <- lapply(spp_to_model, biomod2_wrapper)

#alpha diversity maps
f_em_wmean_bin_current <- paste(spp_to_model,
                                 "/proj_current/individual_projections/", 
                                 spp_to_model,
                                 "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
                                 sep = "")
 
 ### sum all projections
 if( length(f_em_wmean_bin_current) >= 2 ){
   ## initialisation
   taxo_alpha_div_current <- raster( f_em_wmean_bin_current[1] ) 
   for(f in f_em_wmean_bin_current[-1]){
     taxo_alpha_div_current <- taxo_alpha_div_current + raster(f)
   }
 }
 
 ### mask by environemtnal mask
 taxo_alpha_div_current <- mask(taxo_alpha_div_current, subset(stk_current,1))
 
 ## 2050 conditons
 ### load binaries projections
 f_em_wmean_bin_2050 <- paste(spp_to_model,
                              "/proj_2050_BC_45/individual_projections/", 
                              spp_to_model,
                              "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
                              sep = "")
 
 ### sum all projections
 if( length(f_em_wmean_bin_2050) >= 2 ){
   ## initialisation
   taxo_alpha_div_2050 <- raster( f_em_wmean_bin_2050[1] ) 
   for(f in f_em_wmean_bin_2050[-1]){
     taxo_alpha_div_2050 <- taxo_alpha_div_2050 + raster(f)
   }
 }
 
 ### mask by environemtnal mask
 taxo_alpha_div_2050 <- mask(taxo_alpha_div_2050, subset(stk_2050_BC_45,1))
 
 ## 2070 conditons
 ### load binaries projections
 f_em_wmean_bin_2070 <- paste(spp_to_model,
                              "/proj_2070_BC_45//individual_projections/", 
                              spp_to_model,
                              "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
                              sep = "")
 
 ### sum all projections
 if( length(f_em_wmean_bin_2070) >= 2 ){
   ## initialisation
   taxo_alpha_div_2070 <- raster( f_em_wmean_bin_2070[1] ) 
   for(f in f_em_wmean_bin_2070[-1]){
     taxo_alpha_div_2070 <- taxo_alpha_div_2070 + raster(f)
   }
 }
 
 ### mask by environmental mask
 taxo_alpha_div_2070 <- mask(taxo_alpha_div_2070, subset(stk_2070_BC_45,1))

 ## plot the results
 levelplot( stack( c(current = taxo_alpha_div_current, 
                     in_2050 = taxo_alpha_div_2050, 
                     in_2070 = taxo_alpha_div_2070 ) ),
            main = expression(paste('Swartzia ',alpha, "-diversity")),
            par.settings = BuRdTheme)

 levelplot(taxo_alpha_div_current,
            main = expression(paste('Swartzia ',alpha, "-diversity")),
            par.settings = BuRdTheme)
 
