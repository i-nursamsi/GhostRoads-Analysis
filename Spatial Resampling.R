### Spatially re-sample all captures and covariates

rm(list = ls())

library(tidyverse)
library(sf)
library(chron)

### Set Working directories
## create a relative path
p = paste("~/Dropbox/") #this will change for everybody 

setwd(paste(p,"/CT capture histories database/", sep = "")) # this will not change for everybody
rm(p)


## Import metadata
# meta = read.csv("Asian clean CT data/ECL and Collaborator Camera Trap Metadata_20220802.csv")
meta = read.csv("AllCovs_WITHROADS_IN_220902.csv")
meta$X= NULL
head(meta)

## RIGHT AWAY, dont want to include ANY baited cams 
# meta = meta[meta$baited != "Yes",]

caps = read.csv("Asian clean CT data/ECL and Collaborator standardized and independent all captures_20220809.csv")
caps = caps[caps$camera_id %in% meta$camera_id,]

add = read.csv("Asian clean CT data/ECL and Collaborator standardized camera trap deployment data_20220727.csv")
add = add[add$camera_id %in% meta$camera_id,]
add = dplyr::select(add, camera_id, survey_id, Landscape)

test = merge(meta, add, by = "camera_id")
meta =test
rm(test)

## Create a landscape abbreviation col so cell ID isn't so long
sort(unique(meta$Landscape))
meta$land = NA

# i = unique(meta$Landscape)[1]
for(i in unique(meta$Landscape)){
  
  # Split apart landscape name
  s = str_split(i, "_")
  s = s[[1]]
  
  #Only want first letter
  l = paste(substr(s,1,1), collapse = "")
  
  #save it!
  meta$land[meta$Landscape == i]= l
  
}
rm(i,l,s)
unique(meta$land)
# #inspect some odd balls
# unique(meta$Landscape[meta$land == "SM"])
# meta$land[meta$Landscape == "Sarawak_Mulu"] = "Mulu"
# meta$land[meta$Landscape == "Sarawak_Maliau"] = "Maliau"
# meta$land[meta$Landscape == "Singapore"] = "SG"
# 
# unique(meta$Landscape[meta$land == "MS"])
# meta$land[meta$Landscape == "Menglun_Subreserve"] = "Menglun"
# meta$land[meta$Landscape == "Mengao_Subreserve"] = "Mengao"
# meta$land[meta$Landscape == "Mengla_Subreserve"] = "Mengla"
# should be good for now




########### Ultra-Loop to Generate Sampling Units #######

### This loop is fast! :)

## for loop to generate spatially delimited hexagonal grid cells per landscape
## In order to to group nearby cameras

# Create an empty dataframe to store results and add extra cols for new results. 
res = as.data.frame(matrix(NA, nrow = 0, ncol = length(colnames(meta))+4))
names(res) = c(names(meta), "cell_id_1km", "cell_id_3km", "cell_id_5km", "cell_id_10km")


p = "BBS_01Z_CAM_15" #problem cam that intesects two cells at 3km resolution
# p = "07_sorted" #problem cam at 5km res
i = 4

for(i in 1:length(unique(meta$Landscape))){ #for each unique landscape 
  
  #select a single landscape and subset the data
  l = unique(meta$Landscape)[i]
  dat = meta[meta$Landscape == l,]
  
  ## Convert metadata to a spatial object using the coordinates
  shape = st_as_sf(dat, coords = c("Long", "Lat"), 
                   crs = "+proj=longlat +datum=WGS84")
  
  # Ensure EPSG is 4087 via transformation
  shape = st_transform(shape, 4087)
  
  ### Make Hexagonal grids w/ with areas that cover 1km^2, 3km^2, 5km^2, and 10 km^2
  ## cellsize = hexagon short diagonal: https://www.omnicalculator.com/math/hexagon
  ## For 1km^2 hexagon, use cellsize = 1074.6
  ## For 3km^2 hexagon, use cellsize = 1861.2
  ## For 5km^2 hexagon, use cellsize = 2403
  ## For 10km^2 hexagon, use cellsize = 3398
  
  ## 1km hex
  hex1 = st_make_grid(shape, cellsize = 1074.6, square = FALSE) %>% #make the hexagon grid
    st_sf() %>% #convert to simple feature object for easier manipulation
    rowid_to_column("cell_id_1km") #create ID column called cell_id
  
  # Add landscape name to cell_id for all hex's 
  hex1$cell_id_1km = paste0(unique(shape$land), "_", hex1$cell_id_1km)
  
  ## 3km hex
  hex3 = st_make_grid(shape, cellsize = 1861.2, square = FALSE) %>%
    st_sf() %>% 
    rowid_to_column("cell_id_3km")
  # Add landscape name to cell_id
  hex3$cell_id_3km = paste0(unique(shape$land), "_", hex3$cell_id_3km)
  
  ## 5km hex
  hex5 = st_make_grid(shape, cellsize = 2403, square = FALSE) %>%
    st_sf() %>% 
    rowid_to_column("cell_id_5km")
  # Add landscape name to cell_id
  hex5$cell_id_5km = paste0(unique(shape$land), "_", hex5$cell_id_5km)
  
  ## 10km hex
  hex10 = st_make_grid(shape, cellsize = 3398, square = FALSE) %>%
    st_sf() %>% 
    rowid_to_column("cell_id_10km")
  # Add landscape name to cell_id
  hex10$cell_id_10km = paste0(unique(shape$land), "_", hex10$cell_id_10km)
  
  # # # Simple example of how to visualize hexagons and points together.  
  # ggplot()+
  #   geom_sf(data = hex5)+
  #   geom_sf(data = shape)#, aes(color = source))#[shape$camera_id == p,])
  
  # Join hexagons and points by matching points that intersect polygons and remove the geometry
  t1 = st_join(shape, hex1, join = st_intersects) %>%
    st_set_geometry(NULL)
  t1 = dplyr::select(t1, camera_id, cell_id_1km) #only interested in cam_id and cell_id
  
  t3 = st_join(shape, hex3, join = st_intersects) %>%
    st_set_geometry(NULL)
  t3 = dplyr::select(t3, camera_id, cell_id_3km)
  
  t5 = st_join(shape, hex5, join = st_intersects) %>%
    st_set_geometry(NULL)
  t5 = dplyr::select(t5, camera_id, cell_id_5km)
  
  t10 = st_join(shape, hex10, join = st_intersects) %>%
    st_set_geometry(NULL)
  t10 = dplyr::select(t10, camera_id, cell_id_10km)
  
  ## Merge all together on the basis of camera_id
  t = list(t1,t3,t5,t10) %>%
    reduce(full_join, by = "camera_id")
  
  if(length(unique(duplicated(t$camera_id))) > 1){ #add conditional statement for cams that intersect two cells
    
    ## Isolate the cams w/ duplicated cells
    du = t$camera_id[duplicated(t$camera_id)==T]
    
    ## Create empty df to fill in new (single) cell_id's
    r = t[0,]
    
    for(u in 1:length(du)){ # Repeat for each duplicated cell id
      
      # select a single unique cam
      dup = t[t$camera_id == du[u],]
      
      #just take the first cell_id, for simplicity sake. 
      dupl = dup[1,] 
      
      #save it
      r = rbind(r, dupl)
      
    }
    
    ## remove duplicated cams from t
    t = t[!t$camera_id %in% r$camera_id,]
    
    ## and bind new single cell id
    t = rbind(r, t)
    
  }
  
  ## run distinct to verify there are no repeats
  t = distinct(t)
  
  ## Merge w/ OG metadata
  dat = merge(t, dat, by = "camera_id")
  
  ## Save results!
  res = rbind(res, dat)
  
  # which landscape was finished?
  print(l)
  
}
rm(i,t,t1,t3,t5,t10,l,du,dup,dupl,r,
   hex1,hex3,hex5,hex10,dat, u, shape)


## Inspect
head(res) #good! 
anyNA(res$cell_id_1km) 
anyNA(res$cell_id_3km)
anyNA(res$cell_id_5km)
anyNA(res$cell_id_10km) #False for all! 

length(unique(res$cell_id_1km))
length(unique(res$cell_id_3km))
length(unique(res$cell_id_5km))
length(unique(res$cell_id_10km))


dim(res)  
dim(meta)# We match!


## save it! 
meta = res

rm(res)




################ Load and Prepare Captures and Metadata to be Re-Sampled ###### 


### Load data

# rec = read.csv("Asian clean CT data/ECL and Collaborator standardized and independent all captures_20220727.csv") # calebe's code
# anyNA(rec) # There should never by ANY NA values in the records EVER! 

rec = read.csv("Asian clean CT data/ECL and Collaborator standardized and independent all captures_20220809.csv")
rec = rec[rec$camera_id %in% meta$camera_id,] #thin to match relevant covs
anyNA(rec) # There should never by ANY NA values in the records EVER! 

## add missing start/stop dates to meta
add = distinct(select(rec, camera_id, camera_start.date, camera_end.date))
test = merge(meta, add, by = "camera_id")
meta = test
rm(test, add)

# Create a back-up incase meta gets messed up (Calebe's object name)
ECL_meta = meta 


## First save Hex ID's so we can look at spatially close cams later
ECL_meta$Polygon1km = ECL_meta$cell_id_1km;ECL_meta$Polygon3km = ECL_meta$cell_id_3km;
ECL_meta$Polygon5km = ECL_meta$cell_id_5km;ECL_meta$Polygon10km = ECL_meta$cell_id_10km


## Then add survey_id to cell_id to keep cells temporally separated. 
ECL_meta$cell_id_1km = paste(ECL_meta$cell_id_1km, ECL_meta$survey_id, sep = "_")
ECL_meta$cell_id_3km = paste(ECL_meta$cell_id_3km, ECL_meta$survey_id, sep = "_")
ECL_meta$cell_id_5km = paste(ECL_meta$cell_id_5km, ECL_meta$survey_id, sep = "_")
ECL_meta$cell_id_10km = paste(ECL_meta$cell_id_10km, ECL_meta$survey_id, sep = "_")


length(unique(ECL_meta$cell_id_1km))
length(unique(ECL_meta$cell_id_3km))
length(unique(ECL_meta$cell_id_5km))
length(unique(ECL_meta$cell_id_10km))
# decreasing, good. 



# check which records are missing from meta
setdiff(rec$camera_id, ECL_meta$camera_id) #All the danum_2020 SR cams (baited) --> remove them from caps!
setdiff(ECL_meta$camera_id, rec$camera_id) #0 diff! 

setdiff(rec$survey_id, ECL_meta$survey_id) #0 diff! 
setdiff(ECL_meta$survey_id, rec$survey_id) #0 diff! 

# Delete data without record for now
delete = unique(rec$camera_id)[!unique(rec$camera_id) %in% unique(ECL_meta$camera_id)]
rec = rec[!rec$camera_id %in% delete,]
rm(delete)

## Deletes nothing if 0 diff in setdiff()'s above



### Add sequential Julian date to the records and metadata
rec$seq_date = julian(as.Date(rec$Photo.Date)) ## This is the only needed date transformation
head(rec)
anyNA(rec)

ECL_meta$deployment_seq_date = julian(as.Date(ECL_meta$camera_start.date))
ECL_meta$retrival_seq_date = julian(as.Date(ECL_meta$camera_end.date))









### Add cell ids across all scales to the records

for(i in ECL_meta$camera_id){ #select a single camera_id
  
  # and add the corresponding cell_id at all scales
  rec$cell_id_1km[rec$camera_id == i] = ECL_meta$cell_id_1km[ECL_meta$camera_id == i]
  rec$cell_id_3km[rec$camera_id == i] = ECL_meta$cell_id_3km[ECL_meta$camera_id == i]
  rec$cell_id_5km[rec$camera_id == i] = ECL_meta$cell_id_5km[ECL_meta$camera_id == i]
  rec$cell_id_10km[rec$camera_id == i] = ECL_meta$cell_id_10km[ECL_meta$camera_id == i]
  
  
}
rm(i)

## Ensure there are no NAs
anyNA(rec)





########### Ultra-Loop to Re-Sample the Metadata ####### 

## Select which cols we want to average (i.e. by selecting all except which we dont want!)
avg_col = colnames(ECL_meta)[!colnames(ECL_meta) %in% 
                               c("cell_id_1km","Polygon1km","Cell_Effort","survey_id","camera_id",
                                 "cell_id_3km","Polygon3km", "cell_id_5km","Polygon5km",
                                 "cell_id_10km","Polygon10km", "land", 
                                 "Landscape","camera_start.date","camera_end.date",
                                 "deployment_seq_date","retrival_seq_date")]
## This^^ makes sense when we have many cols, but I only have a few at the moment

# avg_col = c("elevation", "forest_integrity", "human_foot", "Latitude", "Longitude",
#             "Canopy_closure", "Canopy_Height")

## Make a list to save all metadata across all scales
meta.results = list()

# i=3;l=2
for(i in 1:length(names(ECL_meta[grepl("cell_id_", names(ECL_meta))]))){ # Repeat for each spatial scale
  
  ## Select a spatial scale
  s = names(ECL_meta[grepl("cell_id_", names(ECL_meta))])[i]
  
  ## Create an empty list to fill in results per sampling unit per scale
  r = list()
  
  for(l in 1:length(unique(ECL_meta[,s]))){ # Repeat for each sampling unit @ that scale
    
    # Select a single SU
    t = ECL_meta[ECL_meta[[s]] == unique(ECL_meta[,s])[l],] # Use double brackets to index variable col name! 
    
    ## Create a new line to fill in the re-sampled metadata per cell 
    new = as.data.frame(matrix(NA, nrow = 1, ncol = ncol(ECL_meta)))
    colnames(new) = colnames(ECL_meta)
    
    ## Calculate means for specific covarites 
    for(a in 1:length(avg_col)){
      
      ## Add conditional if-else statement to retain NAs for cols that only have NA
      if(all(is.na(t[,colnames(t) == avg_col[a]]))){
        
        new[1,colnames(new) == avg_col[a]] = NA
        
        ## But still keep the column present for rbinding later
        colnames(new)[colnames(new) == avg_col[a]] = paste("Avg", avg_col[a], sep = "_")
        
      }else{
        
        ## if not all NA, calculate the mean
        new[1,colnames(new) == avg_col[a]] = mean(t[,colnames(t) == avg_col[a]], na.rm = T)
        
        ## Change col name to reflect new avg value
        colnames(new)[colnames(new) == avg_col[a]] = paste("Avg", avg_col[a], sep = "_")
      }
    }
    
    #Calculate effort in trap-nights for all cams in cell
    effort = 0
    
    for(b in 1:length(unique(t$camera_id))){
      
      ## Use difftime() to calculate duration cam was active from start-end date. 
      effort = effort + as.numeric(difftime(t$camera_end.date[t$camera_id == unique(t$camera_id)[b]], 
                                            t$camera_start.date[t$camera_id == unique(t$camera_id)[b]], units = "days"))
    }
    
    # Save total effort as cell_effort
    new$Cell_Effort = effort
    
    ### Calculate new start/stop dates to accommodate extra cams in single cells
    ## Sampling begins at the minimum start date for all cams included
    new$Sampling_begin = unique(t$camera_start.date[t$deployment_seq_date == min(t$deployment_seq_date)])
    
    ## And sampling ends at the mamimum end date for all cams included 
    new$Sampling_end = unique(t$camera_end.date[t$retrival_seq_date == max(t$retrival_seq_date)])
    
    ## Make a direct copy of character variables
    new$cameras_included = paste(sort(unique(t$camera_id)),collapse = "-")
    new$survey_id = unique(t$survey_id)
    new$Landscape = unique(t$Landscape)
    new$source = unique(t$source)
    new$camera_type = paste(sort(unique(t$camera_type)), collapse = "-")  # COME HERE! make this an observation level covariate
    new$trail_status = paste(sort(unique(t$trail_status)), collapse = "-") # COME HERE! make this an observation level covariate
    new$notes = paste(sort(unique(t$notes)), collapse = "-")
    new$Forest_type = paste(sort(unique(t$Forest_type)),collapse = "-") # COME HERE!
    new$habitat = paste(sort(unique(t$habitat)), collapse = "-")
    new$baited = paste(sort(unique(t$baited)), collapse = "-")  # COME HERE! make this an observation level covariate
    
    ## Save cell_id name, but variable style to accommodate different scales
    new[,s] = unique(t[[s]])
    
    ## Specify variable polygon scale
    id = str_split(s, "_")[[1]][3]
    p = colnames(t[grepl("Poly", colnames(t))])
    p = p[endsWith(p, id)]
    #and save it
    new[[p]] = unique(t[[p]])
    
    ## create vector of averaged cols
    avg = colnames(new[grepl("Avg", colnames(new))])
    ## potentially come here and remove some values from avg
    
    ## Select only the relevant info
    new = dplyr::select(new, 
                        ## Calculated values
                        all_of(s), all_of(p), all_of(avg), 
                        Cell_Effort, Sampling_begin, Sampling_end, 
                        survey_id, Landscape, cameras_included)
    # 
    # 
    # ## Character cols 
    # cameras_included, survey_id, Landscape, source, 
    # camera_type, trail_status, notes, Forest_type,  
    # habitat, baited)
    
    ## save it! 
    r[[l]] = new
    
    
  }
  
  ## Combine all the elements of the temp list into a df
  res = do.call(rbind, r)
  
  ## Save the entire DF
  meta.results[[i]] = res
  names(meta.results)[i] = s
  
}
## Takes ~5 min to re-sample the metadata (5820 cams) across all spatial scales! 
rm(i,l,a,b,p,s,avg,id,t,res,r,new,avg_col,effort)

head(meta.results$cell_id_1km)
head(meta.results$cell_id_5km)
tail(meta.results$cell_id_10km)
tail(meta.results$cell_id_3km)
## it all looks great!! Should be ready to save as is!






########### Ultra-loop to Re-Sample the Records ######


## Make an empty DF to save all records across all scales
rec.results = list()

## Interesting test: Two tiger detections on July 14th, 2014, same camera 1.25 hours apart
u=1;i = "BBSNP_4181_TEAM_BBS_2013"; d = 16265 ; sp = "Panthera_tigris"
i = "KYNP_34_KhaoYai2019"
## Ideally should find example where same species was detected at 2 diff cams
## in same cell on same day to visualize how re-sampling works. 


start = Sys.time()
for(u in 1:length(names(rec[grepl("cell_id_", names(rec))]))){ # Repeat for each spatial scale
  
  ## Select a spatial scale
  s = names(rec[grepl("cell_id_", names(rec))])[u]
  
  ## Create an empty df to fill in results 
  r = rec[0,]
  
  for(i in unique(rec[,s])){ # Repeat for each sampling unit @ that scale
    
    # Subset data for a single SU
    t = rec[rec[[s]] == i,]
    
    for(d in unique(t$seq_date)){ # Repeat for each day there was a detection
      
      # Subset data for a single day
      t2 = t[t$seq_date == d,]
      
      for(sp in unique(t2$Species)){ # Repeat for each different species detected
        
        # Subset data for the specific species
        t3 = t2[t2$Species == sp,]
        
        
        #create the new line to fill with re-sampled data
        new = as.data.frame(matrix(NA, nrow = 1, ncol = ncol(rec)))
        colnames(new) = colnames(rec)
        
        #Sampling Unit-level information (t)
        new[,s] = unique(t[[s]]) ## Cell_id that can change scales
        new$survey_id = unique(t$survey_id) # survey_id
        
        #Date-level information (t2)
        new$Date = unique(t2$Photo.Date)
        new$active_cams_at_date = paste(sort(unique(t2$camera_id)), collapse = " - ") # useful for referencing w/ OG captures alter
        new$num_cams_active_at_date = length(unique(t2$camera_id)) 
        ## This^^IS a really good observation covaraite in detection mods! 
        ## Its like effort, but varies per date
        
        ### COME HERE and continue to add observation-level covarites
        ## e.g. trail, cam type, baited, etc 
        
        #species-level information (t3)
        new$Species = sp
        new$independent_events = nrow(t3) # Mayb this could also be an observation covariate in detection mods?
        new$total_indiv_records = sum(t3$Individuals)
        ## should sum() or max(), idk 
        
        ## Select only the relevant info
        new = dplyr::select(new, all_of(s), survey_id, active_cams_at_date, 
                            num_cams_active_at_date, Date, Species,
                            total_indiv_records, independent_events)
        
        #save via rbind
        r = rbind(r, new)
        
      } # end per species
    } # end per date
  } # end per sampling unit
  
  ## Save per spatial scale
  rec.results[[u]] = r
  names(rec.results)[u] = s
  
} # end per spatial scale
end = Sys.time()
rm(r,s,c,d,i,s,sp,u,r,t,t2,t3, new)

### Inspect
end-start 
# takes ~2.4 hours to work on entire dataset (As of August 2nd 2022, >198,000 independent records)
rm(end,start)

head(rec.results$cell_id_1km)
anyNA(rec.results$cell_id_1km) # good! 

tail(rec.results$cell_id_3km)
anyNA(rec.results$cell_id_3km) # good! 

head(rec.results$cell_id_5km)
anyNA(rec.results$cell_id_5km) # good! 

tail(rec.results$cell_id_10km)
anyNA(rec.results$cell_id_10km) # good! 




###################### Final Inspection & Save ######

### Ensure cell_ids are the same across all scales

# 1km
setdiff(meta.results$cell_id_1km$cell_id_1km, 
        unique(rec.results$cell_id_1km$cell_id_1km)) # no diff! 
setdiff(unique(rec.results$cell_id_1km$cell_id_1km),
        meta.results$cell_id_1km$cell_id_1km) # no diff! 

# 3km
setdiff(meta.results$cell_id_3km$cell_id_3km, 
        unique(rec.results$cell_id_3km$cell_id_3km)) # no diff! 
setdiff(unique(rec.results$cell_id_3km$cell_id_3km),
        meta.results$cell_id_3km$cell_id_3km) # no diff! 

# 5 km
setdiff(meta.results$cell_id_5km$cell_id_5km, 
        unique(rec.results$cell_id_5km$cell_id_5km)) # no diff! 
setdiff(unique(rec.results$cell_id_5km$cell_id_5km),
        meta.results$cell_id_5km$cell_id_5km) # no diff! 

# 10 km
setdiff(meta.results$cell_id_10km$cell_id_10km, 
        unique(rec.results$cell_id_10km$cell_id_10km)) # no diff! 
setdiff(unique(rec.results$cell_id_10km$cell_id_10km),
        meta.results$cell_id_10km$cell_id_10km) # no diff! 


## Save em!
#1km
write.csv(rec.results$cell_id_1km,
          "Asian clean CT data/For Ilyas/Ilyas roads resampled captures 1km_20220914.csv", row.names = F)
write.csv(meta.results$cell_id_1km,
          "Asian clean CT data/For Ilyas/Ilyas roads resampled metadata 1km_20220914.csv", row.names = F)

#3km
write.csv(rec.results$cell_id_3km,
          "Asian clean CT data/For Ilyas/Ilyas roads resampled captures 3km_20220914.csv", row.names = F)
write.csv(meta.results$cell_id_3km,
          "Asian clean CT data/For Ilyas/Ilyas roads resampled metadata 3km_20220914.csv", row.names = F)

#5km
write.csv(rec.results$cell_id_5km,
          "Asian clean CT data/For Ilyas/Ilyas roads resampled captures 5km_20220914.csv", row.names = F)
write.csv(meta.results$cell_id_5km,
          "Asian clean CT data/For Ilyas/Ilyas roads resampled metadata 5km_20220914.csv", row.names = F)

#10km
write.csv(rec.results$cell_id_10km,
          "Asian clean CT data/For Ilyas/Ilyas roads resampled captures 10km_20220914.csv", row.names = F)
write.csv(meta.results$cell_id_10km,
          "Asian clean CT data/For Ilyas/Ilyas roads resampled metadata 10km_20220914.csv", row.names = F)





