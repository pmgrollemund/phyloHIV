########################################################################### ----
##################### Load Rscript options #####################################
########################################################################### ----
require(rjson)
require(phyloHIV)
require(big.phylo)
require(memisc)
require(stringr)

# args <- commandArgs(trailingOnly=TRUE)
# args <- c("./extjson/options.json")
args <- system.file(package="phyloHIV","extpreprocess","options_preprocess.json")
if(length(args)>0){
 json_file <- args[1]
 json_data <- fromJSON(file=json_file)

 for(i in 1:length(json_data)){
  assign(names(json_data)[i], json_data[[i]])
 }
}

########################################################################### ----
##################### Initialization ###########################################
########################################################################### ----
### Redefine some names and paths ----
path_RDS <- paste(normalizePath(path_RDS),"/",sep="")
path_sequences_meta <- paste(path_RDS,RDS_names[1],sep="")
path_person <- paste(path_RDS,RDS_names[2],sep="")
path_Country <- paste(path_RDS,RDS_names[3],sep="")

########################################################################### ----
##################### Preprocess LANL sequences ################################
########################################################################### ----
### Save the RDS file first version ----
file.copy(path_sequences_meta,
          gsub(x=path_sequences_meta,
               pattern=".rds",
               replacement="_save.rds"))
file.copy(path_person,
          gsub(x=path_person,
               pattern=".rds",
               replacement="_save.rds"))

### Deal with the metadata ----
sequences_meta <- data.table(readRDS(path_sequences_meta))
person         <- data.table(readRDS(path_person))

indexes        <- grep("PR/RT",sequences_meta$seqID)
sequences_meta <- sequences_meta[indexes,]
IDS            <- sequences_meta$newnum
{
 person[,race2 := cases(
  race == "" | race == "(9)Unknown"                 -> "Unknown",
  race == "(1)Hispanic, All races"                  -> "Hispanic",
  race == "(2)Not Hispanic, Amer Indian/Alaska Nat" -> "Amer. Ind.",
  race == "(3)Not Hispanic, Asian"                  -> "Asian",
  race == "(4)Not Hispanic, Black"                  -> "Black",
  race == "(5)Not Hispanic, Nat Hawaiian/Pac Isl"   -> "Pac. Isl.",
  race == "(6)Not Hispanic, White"                  -> "White",
  race == "(8)Not Hispanic, Multi-race"             -> "Multi-race"
 )]

 main_races <- names(which(table(person$race2)>min_occur))
 person[  race2 %in% main_races , main_race := race2  ]
 person[!(race2 %in% main_races), main_race := "Other"]

 birthCou_tmp <- substring(person$birthCountry,2,4)
 {
  person[,main_birthCountry := cases(
   birthCou_tmp %in% c("AFG","AZE","BHR","GEO","IRN","IRQ","ISR","KAZ","LBN","LKA",
                       "PAK","QAT","SAU","SLE","SYR","TUR","UZB","YEM"
   )  -> "Middel East",
   birthCou_tmp %in% c("AGO","BDI","BEN","BFA","BWA",'CIV',"CMR","COD","COG","ERI",
                       "ETH","GAB","GHA","GIN","GMB","GNB","KEN","LBR","LSO","MLI",
                       "MOZ","MWI","NAM","NER","NGA","RWA","SDN","SEN","SOM","SWZ",
                       "TCD","TGO","TON","TZA","UGA","ZAF",'ZMB',"ZWE","MRT"
   )  -> "SSA",
   birthCou_tmp %in% c("ALB","AUT","BEL","BGR","BIH","CHE","CZE","DEU","DJI","DNK",
                       "ESP","FIN","FRA","GBR","GRC","HUN","IRL","ISL","ITA","LUX",
                       "LVA","MDA","NLD","NOR","POL","PRT","ROU","RUS","SVK","SWE",
                       "UKR","YUG"
   )  -> "Europe",
   birthCou_tmp %in% c("ARG","ATG","BHS","BLZ","BOL","BRA",'BRB',"CHL","COL","CRI",
                       "CUB","DMA","DOM","ECU","GTM","GUF","GUY","HND","HTI","JAM",
                       "NIC","PAN","PER","PRI","SLV",'TTO',"URY","VEN","VIR","MTQ"
   )  -> "S/M America",
   birthCou_tmp %in% c("ASM","AUS","FJI","FSM","IDN","NZL","PHL","PNG","PRY"
   )  -> "Ocenia",
   birthCou_tmp %in% c("BMU","GUM","MHL","MNP","UMI","WSM","X98"
   )  -> "Other",
   birthCou_tmp %in% c("CAN","MEX"
   )  -> "Neighbors",
   birthCou_tmp %in% c("CHN","HKG","IND","JPN","KHM","KOR","LAO","MMR","MNG","MYS",
                       "NPL","SGP","THA","TWN","VNM"
   )  -> "Asia",
   birthCou_tmp %in% c("DZA","EGY","LBY","MAR","PRK"
   )  -> "North Africa",
   birthCou_tmp %in% c("USA"
   )  -> "USA",
   birthCou_tmp %in% c("X99",""
   )  -> "Unknown"
  )]
  rm(birthCou_tmp)
 }

 person$transm2  <- substring(person$transm,3)
 person[transm2 == "HETEROp",transm2 := "HETERO"]

 main_transms <- names(which(table(person$transm2)>min_occur))
 person[,main_transm := cases(
  transm2 %in% main_transms  -> transm2,
  !(transm2 %in% main_transms) -> "OTHER"
 )]

 person[,main_transm2 := cases(
  transm2 %in% c("IDU","MSM/IDU")  -> "IDU",
  transm2 == "MSM"                 -> "MSM",
  transm2 == "HETERO" & sex == "M" -> "HSX_M",
  transm2 == "HETERO" & sex == "F" -> "HSX_F",
  !(transm2 %in% c("IDU","MSM/IDU","MSM","HETERO")) -> "OTHER"
 )]
}

person$cur_county_name <- str_to_lower(person$cur_county_name)
person$cur_county_name <- gsub(person$cur_county_name,
                               pattern=" co.",replacement="")
person$rsh_county_name <- str_to_lower(person$rsh_county_name)
person$rsh_county_name <- gsub(person$rsh_county_name,
                               pattern=" co.",replacement="")
person$rsa_county_name <- str_to_lower(person$rsa_county_name)
person$rsa_county_name <- gsub(person$rsa_county_name,
                               pattern=" co.",replacement="")
person[,kc_case :=
        cur_county_name == "king" ||
        rsh_county_name == "king" ||
        rsa_county_name == "king"
       ]

saveRDS(person,file=path_person)
saveRDS(sequences_meta,file=path_sequences_meta)

### Country database ----
{
 Europe <- c("DE","CZ","FR","CH","RS","BE","DK","SE","NO","AT","CY",'FI',"GB",
             "SI","RU","PL","NL","IT","ES","GR","AL","ME","RO","LT","PT","LU",
             "SK","BY","IE","BG","UA")
 Asia   <- c("CN","JP","KR","HK","TH","VN","TW","SG","MM","MY","IN","NP")
 Mid_Or <- c("KW","TR","PK","AF","IL","YE","GE","AM")
 Ocean  <- c("AU","ID","PH")
 Nor.Af <- c("TN")
 SSA    <- c("CF","CM","SN","TD","GA","CD","UG","KE","TZ","ZM","RW",'BI',"BW",
             "ZA","MR","TG","NG","BJ","CG","MZ","MW","ZW","ET","DJ","ER","ML",
             "SO","GQ")
 Neighb <- c("CA","MX","GL")
 SM.Am. <- c("BR","CU","PE","VE","UY","AR","CL","DO","HN","SV","JM","PA","CO",
             "PR","TT","BS")
 Unkno. <- c("-")
 US     <- c("US")
 Country_db <- list(Europe = Europe,
                    Asia   = Asia,
                    Mid_Or = Mid_Or,
                    Ocean  = Ocean,
                    Nor.Af = Nor.Af,
                    SSA    = SSA,
                    Neighb = Neighb,
                    SM.Am. = SM.Am.,
                    Unkno. = Unkno.,
                    US     = US)
}
saveRDS(Country_db,file=path_Country)
