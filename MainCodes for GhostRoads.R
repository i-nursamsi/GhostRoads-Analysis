#### Roads Vs Wildlife Abundance
### Full models 
## Started September 27th, 2021
# Author: Ilyas Nursamsi

#install.packages("tidyverse")
install.packages("unmarked")
install.packages("chron")
install.packages("plyr")
install.packages("vegan")
install.packages("AICcmodavg")
install.packages("dplyr")

library(tidyverse)
library(unmarked)
library(chron)
library(plyr)
library(vegan)
library(AICcmodavg)
library(dplyr)

rm(list = ls())
rm(caps, meta)


## working directory for UQ computer
setwd("C:/Users/s4626542/Dropbox/Ghost roads vs Wildlife Abundance/Roads vs Abundance- Analysis/GhostroadVSabundance_R_220902")

##Working directory for laptop
#setwd("C:/Users/ilyas/Dropbox/Ghost roads vs Wildlife Abundance/Roads vs Abundance- Analysis/GhostroadVSabundance_R_220902") 
###### Load Custom Shrink Function #####

shrink.pcount <- function(matrice, nday) {
  dy <- nday
  while (dy < ncol(matrice)) {
    dy <- dy + nday
  }
  addcol <- dy - ncol(matrice)
  if (addcol != 0) {
    matNA <- matrix(NA, nrow = nrow(matrice), ncol = addcol)
    matrice <- data.frame(matrice, matNA)
  }
  
  period <- ncol(matrice)/nday
  newday <- rep(1:period, each = nday)
  
  shr <- function(vec) {
    nav <- is.na(vec)
    dom <- all(nav == T)
    if (dom == T) {
      y <- NA
    } else {
      abb <- sum(vec, na.rm = T)
      y <- ifelse(abb == 0, 0, abb)
    }
    return(y)
  }
  
  matday <- data.frame(newday, t(matrice))
  shrmat <- t(aggregate(matday[, -1], list(matday$newday), shr))
  
  return(shrmat[-1, ])
} 
#modified- Rovero and Zimmerman function to combine x days into 1 observation occasion, but modified to sum
#the number of defections rather than replace the value with a 1. 


######### Import and Clean Data #######

##### Captures
caps= read.csv("Ilyas roads resampled captures 1km_20220914.csv")
head(caps)
tail(caps)
str(caps) 


## Remove sampling units that were only active for 1 day
#caps = caps[!caps$cell_id_1km %in% c("263041_Pasoh_TEAM_2014", "853040_Danum_Valley_2019a"),]



## This new resamled captures does not contain Landscape, sampling_begin and sampling_end so i grab it from the metadata WITH THE SAME RESAMPLING SIZE
#call the resampled metadata that contains the sampling_begin and sampling_end
meta <- read.csv('Ilyas roads resampled metadata 1km_20220914.csv', header = TRUE, sep = ",")

#just merge and thin later, cuz i dont know better way lol
meta <- meta[ , c("cell_id_1km", "Sampling_begin", "Sampling_end", "Landscape")]
head(meta)
caps = merge(caps, meta, by = 'cell_id_1km')
head(caps)



#save the resampled with begin date just in case?
write.csv(caps, 'Ilyas roads resampled metadata 1km_20221116_WITHbeginEnd.csv', row.names=FALSE)


head(caps)
#Need to make dates as dates
caps<- caps %>%
  mutate(Date= as.Date(Date), Sampling_begin = as.Date(Sampling_begin), Sampling_end = as.Date(Sampling_end))


## check the caps again before thinning
head(caps)
## Thin caps to only the relevant cols
caps = dplyr::select(caps, cell_id_1km, Landscape, survey_id, Date, Species, total_indiv_records, Sampling_begin, Sampling_end)

rm(meta)

head(caps)

####### Metadata
meta<-read.csv('Ilyas roads resampled metadata 1km_20220914.csv', header=TRUE, sep=",")
head(meta)
str(meta)
names(meta)
meta<- mutate_if(meta, is.factor, as.character)

## Remove sampling units that were only active for 1 day
#meta = meta[!meta$cell_id_1km %in% c("263041_Pasoh_TEAM_2014", "853040_Danum_Valley_2019a"),]


# meta$survey_id[startsWith(meta$survey_id, "Sing")]= "Singapore"
sort(unique(meta$survey_id)) #2 surveys for now! 

## Create a shrink_id due to shrink mess up
meta$shrink_id = paste0("X", meta$cell_id_1km)

## Bring in Landscape factor covariate from caps
meta = merge(meta, distinct(caps[which(colnames(caps) %in% c("survey_id"))]), by = "survey_id")

#checks
head(meta)
## im just using the 1 and 3km for now
# Thin down to relevant columns add the circuit theory columns here
names(meta)
head(caps)
meta = dplyr::select(meta, cell_id_1km, shrink_id, survey_id, Landscape, Cell_Effort,
                     Avg_RdCount_GHOST_1000, Avg_RdCount_OSM_1000, Avg_RdCount_GRIP_1000,
                     Avg_RdLengthkm_GHOST_1000, Avg_RdLengthkm_OSM_1000, Avg_RdLengthkm_GRIP_1000, Avg_hunting_accessibility_1km)

## check again!!
head(meta)
names(meta)
# potentially change this in the future.
# We want percentages calculates at scale most close to Sampling Unit scale
## 3km data == 7.8 km2 area, 5km data == 21.65 km2 area.
## 1000 m radius covers an area of 3.14 km2, and 2000 covers an area of 12.5 km2
## 1500 m radius would be ideal.... OR 2500 m radius for 5km data
# Just going with 1000 for now b/c small and local scale analysis.


### Inspect how correlated Road vars are
cov.corr = cor(meta[,sapply(meta, is.numeric)])
cov.corr = as.data.frame(cov.corr)
cov.corr[abs(cov.corr)< abs(.6)] = "good"
cov.corr
# Not as correlated as I would have guessed! 

summary(meta) # no NA, good! 


## Thin captures to only those sites where we have road data... For now! 
caps = caps[caps$cell_id_1km %in% meta$cell_id_1km,]
unique(caps$survey_id) #all good! 


rm(cov.corr)

###### Generate Sampling Occasion Index #########

### Create a data frame with each sampling unit as a row, with start and stop dates
s = distinct(select(caps, cell_id_1km, Sampling_begin, Sampling_end))  

## split the dataframe by cell_id_1km, and add the full sequence of dates between start/stop dates. 
s2 = ddply(s, .(cell_id_1km), summarize,
           Date = as.character(seq.Date(from = min(Sampling_begin), to = (max(Sampling_end)+1), by = 1)))
## Add +1 to end because date was missing the last sampling day

#That worked, now bring back the start and stops
dim(s2) #38506 rows
t = merge(s2, s, by = "cell_id_1km")
dim(t) #38506, nice! 


## Add sequence from 1-n for each sampling unit
res = t[0,]
res$seq = numeric()

for(i in unique(t$cell_id_1km)){
  
  d = t[t$cell_id_1km == i,]
  
  d$seq = as.numeric(seq(from= 0, 
                         to = unique(as.numeric(difftime(d$Sampling_end+1, #+1 to accommodate the actual last day 
                                                         d$Sampling_begin, units = "days"))), by = 1))
  
  res = rbind(res, d)
}
rm(d,i)
head(res) 
str(res)
dim(res) #38506 rows --> same as before! Good! 



## Make sure all cams got accounted for
setdiff(res$cell_id_1km, unique(caps$cell_id_1km))
setdiff(unique(caps$cell_id_1km), res$cell_id_1km) # no difference

dim(caps)

## merge the sequence with the captures
t = merge(caps, res, by = c("cell_id_1km", "Date", "Sampling_begin", "Sampling_end"))
dim(caps) #22684
dim(t) #22684, good

## Captures are all good with updated observation sequence! 
caps = t

dim(caps)
head(caps)
tail(caps)
rm(t, s, s2, res)


sort(unique(caps$Species))


#look at the number of detected individue for each of the species within the landscapes


######## Loop to Create UMFs ##########

# Select which Species you want
table(caps$Species)

df <- table(caps$Species)



Species = c("Muntiacus_muntjak","Sus_scrofa", "Rusa_unicolor", "Macaca_nemestrina",
            "Tragulus_sp.", "Hystrix_brachyura", "Argusianus_argus", "Maxomys_sp.",
            "Muntiacus_reevesi", "Hemigalus_derbyanus", "Leopoldamys_sabanus", "Phasianidae_spp",
            "Lophura_ignita", "Atherurus_macrourus", "Tapirus_indicus", "Prionailurus_bengalensis", 
            "Helarctos_malayanus", "Echinosorex_gymnura", "Paradoxurus_musangus_or_hermaphroditus",
            "Sciuridae_spp", "Elephas_maximus", "Martes_flavigula", "Muntiacus_sp.", "Macaca_fascicularis",
            "Macaca_arctoides", "Sus_barbatus", "Trichys_fasciculata", "Herpestes_urva",
            "Neofelis_nebulosa", "Manis_javanica", "Paguma_larvata")

## Come here and change in the future! 


###### Loop to generate Count History Matrices and create UMFs

umf.list = list()

# t=2; j=1

for(t in 1:length(Species)){ #run this loop for each Species
  
  sp = (Species)[t]
  
  ### Specify which Landscape the Species has been detected (at least once!)  
  lands = unique(caps$Landscape[caps$Species == sp])
  # lands = unique(caps$survey_id[caps$Species == sp])
  
  ## Select relevant Landscapes
  c = caps[caps$Landscape %in% lands,] #subset captures
  m = meta[meta$Landscape %in% lands,] #subset metadata
  
  # c = caps[caps$survey_id %in% lands,] #subset captures
  # m = meta[meta$survey_id %in% lands,] #subset metadata
  
  
  ## standardize metadata
  m.num<- m[,sapply(m, is.numeric)] 
  m.std<- decostand(m.num, method = "standardize", na.rm = TRUE)
  m.std = m.std[,colSums(is.na(m.std)) < nrow(m.std)] #remove any columns with no data
  m.char<- m[,sapply(m, is.character)]
  m<- data.frame(m.char, m.std)
  
  ## Outline the structure of the matrix and add col and row names
  mat = matrix(NA, 
               nrow = length(unique(c$cell_id_1km)), 
               ncol = length(seq(from =0, to= max(c$seq))),
               dimnames = list(as.character(unique(c$cell_id_1km)), seq(from =0, to= max(c$seq))))
  
  ## Determine when each sampling unit was active-
  for(j in 1:length(unique(c$cell_id_1km))){
    
    a= c[c$cell_id_1km == unique(c$cell_id_1km)[j],] #subset for a specific sampling unit out of all possible options
    
    indx = seq(from = min(a$seq), to = max(a$seq)) #determine the sequence it was operational for 
    
    mat[j,indx]=0 # at row J, across all sampling occasions, put a zero there
  }
  
  ## Fill in the matrix based on sampling unit and sampling occasion
  for(j in 1:length(unique(c$cell_id_1km))){ #repeat for each sampling unit
    
    su = unique(c$cell_id_1km)[j] #specify the sampling unit
    
    a = c[c$cell_id_1km == su & c$Species == sp,] #subset captures data for specific sampling unit and for specific Species
    
    if(nrow(a)>0){ #Bypass cameras without a detection and leave them as zero
      
      for(s in 1:length(a$Date)){ #repeat for each Date detected at the sampling unit
        
        d = a$Date[s] #specify the sampling date
        
        indx = a$seq[a$cell_id_1km == su & a$Date == d]
        
        mat[su,indx]= unique(a$total_indiv_records[a$Date == d]) #in the matrix where the row = sampling unit, and column = occasion, 
        #use the total counts of the specific sampling occasion
        ## Added unique to bypass double detection bug of scrofa at "327429_BBS_frag4" on 2014-07-06
      }
    }
  }
  
  # shrink the matrix
  dh.mat = shrink.pcount(mat, 5) # using 3 day window b/c thats what Ive used in the past 
  ## COME HERE AND CHANGE! Using a larger shrink now b/c fewer sites and low detections
  
  ## limit the sampling duration in the dh matrix 
  # dh.mat =  dh.mat[, 1:25] #limit maximum sampling duration to 125 days.
  ## COME HERE! Leuser has no sites > 125 days, but others will. 
  
  # make sure meta matches the order as the Detection history matix 
  m = m[order(match(m$shrink_id, rownames(dh.mat))),]
  
  # make the umf
  umf = unmarkedFramePCount(y = dh.mat, siteCovs = m)
  
  umf.list[[t]]= umf
  names(umf.list)[t]= sp
  
  
} #10 warnings about char to factor --> all good! 
## Also surprised with how fast this loop is :)
rm(dh.mat, c, sp, lands, t, a, j, indx, su, d,s, 
   umf, mat, m, m.char, m.num, m.std)

summary(umf.list$Muntiacus_muntjak)
plot(umf.list$Muntiacus_muntjak)



summary(umf.list$Sus_scrofa)
plot(umf.list$Sus_scrofa)

summary(umf.list$Rusa_unicolor)
plot(umf.list$Rusa_unicolor)

###### ORIGINAL SPECIES LIST ########
# "Panthera_tigris","Panthera_pardus", "Muntiacus_muntjak",
# "Sus_scrofa", "Rusa_unicolor", "Capricornis_sp.", "Macaca_nemestrina",
# 'Helarctos_malayanus', "Tragulus_sp."## all worked tho!

# ########## SPECIES AFTER TRIMMING ########
# Species = c("Panthera_tigris","Panthera_pardus", "Muntiacus_muntjak",
#             "Sus_scrofa", "Rusa_unicolor")
# 
# 
########### Run N-Mixture Models ###########

### This loop will test each variable independently for each Species
### AND run a combo model w/ all count and length datas together. 




skip = c("cell_id_1km", "survey_id", "shrink_id","Landscape","Cell_Effort") #dont want pcount mods of these
mod.list = list() #store all models here

# store null predictions here --> not so useful for this project, but leaving legacy code anyway
null.prediction = data.frame(matrix(NA, nrow = 0, ncol = 6))
colnames(null.prediction) = c("Predicted", "SE", "lower", "upper", "Landscape", "Species")

# i=6; l=4
a = Sys.time()
for(i in 1:length(umf.list)){  # Repeat for every predator Species
  
  m = umf.list[[i]] #select a Species
  sp = names(umf.list)[i] #save the Species name
  
  null = pcount(~ Cell_Effort  ~ Landscape, m) #run the null model for the Species
  ## COME HERE
  ## Change survey_id to Landscape when we have more Landscapes w/ data! 
  
  
  #save null predictions
  p = distinct(predict(null, "state",
                       newdata = data.frame("Landscape" = m@siteCovs$Landscape), appendData= T))
  p$Species = sp
  null.prediction = rbind(null.prediction, p)
  
  
  res = list() #store results per Species 
  
  for(l in 1:length(names(m@siteCovs)[!names(m@siteCovs) %in% skip])){ # Repeat for each variable in SiteCovs, 
    ## but skip the covs in skip
    covs = names(m@siteCovs)[!names(m@siteCovs) %in% skip]
    
    form = as.formula(paste0(" ~ Cell_Effort ~ ", covs[l], "+ Landscape", sep = "")) 
    
    mod = pcount(form, m)
    
    res[[l]]= mod 
    names(res)[l]= covs[l]
    
  }
  
  res[[l+1]]= null
  names(res)[l+1]= "null"
  
  c = pcount(~Cell_Effort ~ Avg_RdLengthkm_GRIP_1000 + Avg_RdLengthkm_OSM_1000 + Avg_RdLengthkm_GHOST_1000 + Avg_hunting_accessibility_1km + Landscape, m)
  
  d = pcount(~Cell_Effort ~ Avg_RdCount_GRIP_1000 + Avg_RdCount_OSM_1000 + Avg_RdCount_GHOST_1000 + Avg_hunting_accessibility_1km + Landscape, m)
  
  res[[l+2]] = c
  names(res)[l+2] = "Length_combo"
  
  res[[l+3]] = d
  names(res)[l+3] = "Count_combo"
  
  mod.list[[i]]= res
  names(mod.list)[i]= sp
  
} ## Probably will take a while
b = Sys.time()
b-a # check how long this takes. 
## it took 2.5 hours!!! for these six landscapes


### After you have finished constructing and running all of your models, 
### make sure to save all of your models so you do not have to wait for them again!
save.image(file = "GhostRoads_Sixlandscapes_20230207.RData")
# If you want to load in all of your models and R environment, simply use the load function:
load("GhostRoads_Sixlandscapes_20230207.RData")
# It will take a few minutes to load.

########################################
######################################
#Extract the mean values for each univariate model
data <- read.csv("Preliminary_N-mix_model_coefficents_1km_230120[shrink6_ZIP].csv")
result <- data %>% 
  group_by(model) %>% 
  summarise_at(vars(Estimate, SE, P, AIC), mean)

write.csv(result, "coeff_filtered_univariate_230207[shrink3_ZIP].csv")

mod.list$Muntiacus_muntjak
mod.list$Sus_scrofa
mod.list$Rusa_unicolor


mod.list$Tragulus_sp.
mod.list$Macaca_nemestrina # all worked, but I can tell we lack stat power

### NEW SPECIES LIST AFTER TRIMMING 
# Species = c("Panthera_tigris","Panthera_pardus", "Muntiacus_muntjak",
#             "Sus_scrofa", "Rusa_unicolor")

rm(i,l,sp,mod,null,res,m,form, covs, skip, a, b, p, c,d)




### Extract relevant model info from models. 

coeff = list() #store model values as a DF here
mod.select = list() #store model selection info here

# i=1;l=2
for(i in 1:length(mod.list)){ ## Repeat for each Species
  
  m = mod.list[[i]] #subset for models within 1 Species 
  sp = names(mod.list)[i] #save the Species name
  
  ## Run the model selection functions and save per Species
  a = fitList(fits = m)
  mod.select[[i]]= modSel(a, nullmod = "null")
  names(mod.select)[i]= sp
  
  res = list() #store results for each model here
  
  for(l in 1:length(m)){ #repeat for each model output
    
    ## Extract coefficients from model summary
    r = summary(m[[l]])
    r = r[[1]] #only want state var
    
    ## Specify variables
    r$vars= rownames(r)
    r$model = names(m)[l] # save model name
    
    ## Standardize var names and use model name to fill in later on graph.
    # But do this differently for multi vs univariate models 
    
    if(grepl("combo", names(m)[l])){ # only do this for multi-variate mods
      
      t = separate(r, vars, into = c("a",'b',"c"), sep = "_")
      t$source = t$b
      t = select(t, -c(a,c)) %>%
        dplyr::rename(vars = b)
      # t$model = paste(t$model, t$vars, sep = "_")
      t$vars[2:4] = "var"
      t$vars[c(1,5)] = rownames(t)[c(1,5)]
      
      r = t
      
    }
    
    if(unique(r$model) != "null" & !grepl("combo", unique(r$model))){
      
      r$source = str_split(r$model[2], "_")[[1]][2][3] # save the data source! 
      r$vars[2]= "var"
      
    }
    
    if(unique(r$model) == "null"){
      
      r$source = "NA"
    }
    
    ## Save the Species name
    r$Species = sp
    
    ## Mark which vars are significant
    names(r)[4]= "P"
    r$sig = "Non-Significant"
    r$sig[r$P < 0.05]= "Significant"
    
    ## Make sure all model diagnostics are good
    r$convergence = m[[l]]@opt$convergence #should be 0 if it converged
    r$AIC = m[[l]]@AIC
    
    ## Run a Goodness of Fit test from AICmodavg package
    g = Nmix.chisq(m[[l]])
    r$chi_sq = g$chi.square
    
    ## Save results
    res[[l]]= r
    names(res)[l]= unique(r$model)
    
  }
  
  coeff[[i]]= res
  names(coeff)[i]= sp
  
  
} 
rm(i, l, r, a, m, sp, res, g, t)

(coeff$Muntiacus_muntjak$Length_combo)
(coeff$Muntiacus_muntjak$Avg_RdCount_GHOST_1000)


## Turn coeff into a df
a = do.call(rbind, coeff$Muntiacus_muntjak) 
coeff.res = a[0,]

for(i in 1:length(coeff)){
  a = do.call(rbind, coeff[[i]])
  coeff.res = rbind(coeff.res, a)
} 
rm(a, i)
rownames(coeff.res)= NULL
head(coeff.res)
coeff.res$combo = paste0(coeff.res$Species, "~", coeff.res$model, sep= "")

## Check which models didnt coverge
coeff.res$combo[coeff.res$convergence>0] #Tiger length combo model FAILED! 

# Save model coefficents
write.csv(coeff.res, "Preliminary_N-mix_model_coefficents_20221116.csv", row.names = FALSE)



###### Visualize Model Coefficients ############


######## Univariate Bar charts

### loop thru all spp 
coeff.uni.plots = list()

for(i in 1:length(unique(coeff.res$Species))){
  
  sp = unique(coeff.res$Species)[i]
  
  ## Reduce dataset to relevant info
  a = coeff.res[coeff.res$vars == "var" & 
                  !grepl("combo", coeff.res$model) &
                  coeff.res$Species == sp,]
  
  p =
    ggplot(a, aes(y = Estimate, x = reorder(model, -Estimate)))+
    geom_col(aes(fill= source), color = "black")+
    geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0)+
    geom_text(aes(y = Estimate+.1*sign(Estimate), label = ifelse(P < 0.05, "*", "")),
              size = 10)+
    theme_classic()+
    geom_hline(aes(yintercept = 0))+
    labs(y = "Estimate", x = NULL, fill = "Data Source", title = sp)+
    scale_x_discrete(labels = c("Avg_RdCount_OSM_1000" = "OSM Road Density 1km",
                                "Avg_RdCount_GRIP_1000" = "GRIP Road Density 1km",
                                "Avg_RdLengthkm_GRIP_1000" = "GRIP Road Length 1km",
                                "Avg_RdCount_GHOST_1000" = "GHOST Road Density 1km",
                                "Avg_RdLengthkm_GHOST_1000" = "GHOST Road Length 1km",
                                "Avg_RdLengthkm_OSM_1000" = "GHOST Road Length 1km",
                                "Avg_hunting_accessibility_1km" = "Hunting accessibility 1km"))+
    theme(axis.text.x = element_text(angle = 20, vjust = 0.9, hjust= .95))
  
  coeff.uni.plots[[i]]= p
  names(coeff.uni.plots)[i]= sp
  
}
rm(i, p, sp,a)

coeff.uni.plots$Muntiacus_muntjak
coeff.uni.plots$Rusa_unicolor

coeff.uni.plots$Panthera_tigris




###### Multi-variate Bar Charts

### loop thru all spp 
coeff.combo.plots = list()

for(i in 1:length(unique(coeff.res$Species))){
  
  sp = unique(coeff.res$Species)[i]
  
  ## Reduce dataset to relevant info
  a = coeff.res[coeff.res$vars == "var" & 
                  grepl("combo", coeff.res$model) &
                  coeff.res$Species == sp,]
  
  temp = list()
  for(l in 1:length(unique(a$model))){ #repeat for length and density
    
    b = a[a$model == unique(a$model)[l],]
    title = paste(sp, unique(a$model)[l], sep = "~")
    
    p =
      ggplot(b, aes(y = Estimate, x = reorder(source, -Estimate)))+
      geom_col(aes(fill= source), color = "black")+
      geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0)+
      geom_text(aes(y = Estimate+.1*sign(Estimate), label = ifelse(P < 0.05, "*", "")),
                size = 10)+
      theme_classic()+
      geom_hline(aes(yintercept = 0))+
      labs(y = "Estimate", x = NULL, fill = "Data Source", title = title)
    
    temp[[l]] = p
    names(temp)[l] = unique(a$model)[l]
    
  }
  
  coeff.combo.plots[[i]]= temp
  names(coeff.combo.plots)[i]= sp
  
}
rm(i, p, sp,a, title, l, temp)

coeff.combo.plots$Muntiacus_muntjak$Count_combo
coeff.combo.plots$Muntiacus_muntjak$Length_combo



#### SAVE BOTH!! 
setwd("C:/Users/s4626542/Dropbox/Ghost roads vs Wildlife Abundance/Roads vs Abundance- Figures")

for(i in 1:length(coeff.uni.plots)){
  
  a = coeff.uni.plots[[i]]
  sp = names(coeff.uni.plots)[i]
  
  path = paste0("Coefficent_barplots_univariate/", sp, "_sixlandscapes_20221116.png")
  ggsave(path, a, width = 7, height = 6, units = "in")
}
rm(a, sp, i, path)

for(i in 1:length(coeff.combo.plots)){
  
  a = coeff.combo.plots[[i]]
  sp = names(coeff.combo.plots)[i]
  
  for(l in 1:2){
    
    b = a[[l]]
    v = names(a)[l]
    
    path = paste0("Coefficent_barplots_multivariate/", sp, "~", v, "__sixlandscapes_20221116.png")
    ggsave(path, b, width = 7, height = 6, units = "in")
  }
}
rm(a, sp, i, path,b,v,l)




########## Visualize Model Predictions ########

## using the single variable models, I want to plot
## a scatterplot w/ each different data source 

## Produce prediction dataframes 
predict.res.count = list()
predict.res.length = list()

for(i in 1:length(mod.list)){ # repeat for each Species
  
  m = mod.list[[i]]
  sp = names(mod.list)[i]
  
  temp.count = list()
  temp.length = list()
  for(l in 1:length(m)){
    
    n = m[[l]]
    v = names(m)[l]
    
    if(!grepl("combo", v) & v != "null"){
      
      ## extract the survey w/ smallest SE for prediction new data
      s = summary(n)
      s = s$state
      s$var = rownames(s)
      s = s[grepl("Landscape", s$var),] #will miss intercept here! Need to re-run mods w/ -1 on Landscape!!! 
      z = s$var[s$SE == min(s$SE)]
      z = strsplit(z, "Landscape")[[1]][2] #remove Landscape from the name
      
      ## Create new data to predict upon
      newd = data.frame(v = seq(from = min(n@data@siteCovs[,v]), 
                                to = max(n@data@siteCovs[,v]),
                                length.out = 100),
                        "Landscape" = z)
      names(newd)[1]= v #ensure var name sticks 
      
      ## Make your predictions
      p = predict(n, "state", newd, appendData = T)
      
      ## Normalize the simulated data
      o = meta[,v]
      p$norm.cov = p[,5]*sd(o)+mean(o)
      p$combo = paste0(sp, "~", v, sep = "")
      names(p)[5] = "sim.cov"
      
      ## add the Species
      p$Species = sp
      
      ## Determine if relationship is significant or not. 
      s = coeff.res[coeff.res$combo == unique(p$combo),]
      p$sig = s$sig[s$vars == "var"]
      p$sig = factor(p$sig, levels = c("Non-Significant", 
                                       "Significant"))
      
      ## add data source
      p$source[grepl("OSM", p$combo)]= "OSM"
      p$source[grepl("GHOST", p$combo)]= "GHOST"
      p$source[grepl("GRIP", p$combo)]= "GRIP"
      
      
      ## Save in the appropriate list 
      if(grepl("Length", v)){
        
        temp.length[[l-3]] = p
        names(temp.length)[l-3] = v
        
      }
      
      if(grepl("Count", v)){
        
        temp.count[[l]] = p
        names(temp.count)[l] = v
        
        
      }
      
    }
    
  }
  
  ## Save predictions for count mods per Species 
  predict.res.count[[i]] = temp.count
  names(predict.res.count)[i] = sp
  
  ## Save predictions for length mods per Species 
  predict.res.length[[i]] = temp.length
  names(predict.res.length)[i] = sp
  
  #come back here again later
  
}
rm(i,l,o,sp,v,z,s,p,n,m,newd, temp.length, temp.count)

head(predict.res.count$Muntiacus_muntjak$Avg_RdCount_GHOST_1000)
head(predict.res.length$Muntiacus_muntjak$Avg_RdLengthkm_OSM_1000)
## All good! 


### Combine relevant dataframes 
count.dat = predict.res.count$Muntiacus_muntjak$Avg_RdCount_GHOST_1000
count.dat = count.dat[0,]

for(i in 1:length(predict.res.count)){
  
  t = predict.res.count[[i]]
  y = do.call(rbind, t)
  count.dat = rbind(y, count.dat)
}
rm(t,y,i)
unique(count.dat$combo) # all good! 

length.dat = predict.res.length$Muntiacus_muntjak$Avg_RdLengthkm_GHOST_1000
length.dat = length.dat[0,]

for(i in 1:length(predict.res.length)){
  
  t = predict.res.length[[i]]
  y = do.call(rbind, t)
  length.dat = rbind(y, length.dat)
}
rm(t,y,i)
unique(length.dat$combo) #all good! 

dat = list("length" = length.dat, "count" = count.dat)





##### Make the plots! 
p.plots = list()

for(i in 1:2){
  
  d = dat[[i]]
  n = names(dat)[i]
  
  temp = list()
  for(l in 1:length(unique(d$Species))){
    
    sp = unique(d$Species)[l]
    a = d[d$Species == sp,]
    
    if(n == "length"){
      
      if(unique(a$sig) == "Non-Significant"){
        
        p =
          ggplot(a, aes(y = Predicted, x = norm.cov))+
          geom_line(aes(color = source), linetype = "dashed", size = 2)+
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = source), alpha = .4)+
          theme_classic()+
          labs(y = "Predicted Abundance", x = "Road Length (m)", title = sp, color = "Data Source")+
          guides(aes(fill = FALSE, linetype = FALSE))
        
      }else{
        
        p =
          ggplot(a, aes(y = Predicted, x = norm.cov))+
          geom_line(aes(color = source), size = 2)+
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = source), alpha = .4)+
          theme_classic()+
          labs(y = "Predicted Abundance", x = "Road Length (m)", title = sp, color = "Data Source")+
          guides(aes(fill = FALSE, linetype = FALSE))
      }
    } # end length 
    
    if(n == "count"){
      
      if(unique(a$sig) == "Non-Significant"){
        
        p =
          ggplot(a, aes(y = Predicted, x = norm.cov))+
          geom_line(aes(color = source), linetype = "dashed", size = 2)+
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = source), alpha = .4)+
          theme_classic()+
          labs(y = "Predicted Abundance", x = "Road Density in 1km Radius", title = sp, color = "Data Source")+
          guides(aes(fill = FALSE, linetype = FALSE))
        
      }else{
        
        p =
          ggplot(a, aes(y = Predicted, x = norm.cov))+
          geom_line(aes(color = source), size = 2)+
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = source), alpha = .4)+
          theme_classic()+
          labs(y = "Predicted Abundance", x = "Road Density in 1km Radius", title = sp, color = "Data Source")+
          guides(aes(fill = FALSE, linetype = FALSE))
      }
    } # end count 
    
    if(n == "hunting"){
      
      if(unique(a$sig) == "Non-Significant"){
        
        p =
          ggplot(a, aes(y = Predicted, x = norm.cov))+
          geom_line(aes(color = source), linetype = "dashed", size = 2)+
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = source), alpha = .4)+
          theme_classic()+
          labs(y = "Predicted Abundance", x = "Huting Accessibility in 1km Radius", title = sp, color = "Data Source")+
          guides(aes(fill = FALSE, linetype = FALSE))
        
      }else{
        
        p =
          ggplot(a, aes(y = Predicted, x = norm.cov))+
          geom_line(aes(color = source), size = 2)+
          geom_ribbon(aes(ymin = lower, ymax = upper, fill = source), alpha = .4)+
          theme_classic()+
          labs(y = "Predicted Abundance", x = "Huting Accessibility in 1km Radius", title = sp, color = "Data Source")+
          guides(aes(fill = FALSE, linetype = FALSE))
      }
    } 
    
    # save the plots per Species
    temp[[l]] = p
    names(temp)[l] = sp
    
    
  }
  
  # Save the plots per data type 
  p.plots[[i]] = temp
  names(p.plots)[i] = n
  
}
rm(i,l,d,p,sp,temp,n,a)


p.plots$length$Rusa_unicolor
p.plots$count$Macaca_nemestrina


## Save them all! 
for(i in 1:length(p.plots)){
  
  a = p.plots[[i]]
  v = names(p.plots)[i]
  
  for(l in 1:length(a)){
    
    b = a[[l]]
    sp = names(a)[l]
    
    path = paste0("Predicted_abundance_plots/", sp, "~ Road_", v, "_20221116.png")
    ggsave(path, b, width = 7, height = 6, units = "in")
  }
}
rm(a, sp, i, path,b,v,l)



# maybe trying to resample the species methodology here,

resamp = data$Commodities
resamp2 = data$Classification



## try to make it more easier to read


