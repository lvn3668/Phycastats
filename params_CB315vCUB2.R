### ========= input parameters section ===========

## rarely changed params:

    # === image format parameters ===
    # the suffix <FRM> of saveImageAs<FRM> variable stands for
    # an image format. The following list is a subset
    # of all formats supported by a2ping.
    # Here is a complete list of supported a2ping outputs: BMP EPS GIF JPEG PBM PCL5 PDF PDF1 PGM PNG PPM PS TIFF XPM XWD markedEPS markedPS
    # If you think that I've ommited some important format,
    # please let me know or simply add the appropriate flag (and add it to imageFormat table in core.R)
    saveImageAsPDF  <- TRUE
    saveImageAsPS   <- TRUE  # this needs to be TRUE for the latex report file
    saveImageAsPBM  <- FALSE
    saveImageAsPGM  <- FALSE
    saveImageAsPPM  <- FALSE
    saveImageAsPNG  <- FALSE
    saveImageAsXWD  <- FALSE
    saveImageAsBMP  <- FALSE
    saveImageAsTIFF <- FALSE
    saveImageAsJPEG <- FALSE
    saveImageAsGIF  <- FALSE
    saveImageAsXPM  <- FALSE

    doAbtTables <- FALSE

    useHelv <- FALSE

    dumpAbuTables <- TRUE
    dumpAnnTbl <- FALSE # use annotate_tables.py instead of g3anntbl.csv
    dumpBinaryTables <- TRUE
    dumpDistTables <- TRUE
    #within_between_analysis <- TRUE
    alternate_pch_Thld <-1000

## environment/user specific config:
    # === output directory ===
    # try to avoid dashes in directory names, they are interpolated into dots
dataDir <- "/Users/lalithaviswanathan/Documents/Project/Cubi_LChesnel_EPAN12_0120/PhyCA.Stats.Out"  # directory for storing all PhyCA-Stats output data
sysTime <- format(Sys.time(), "%H.%M.%S")                     # directory for storing all this particular run of PhyCA-Stats
outDir <- paste(dataDir,"/",Sys.Date(),"_",sysTime,".PhyCA.Stats.Out",sep="")
outDir <- gsub("-",".",outDir)

nCPUs <- 3  # number of CPUs
# Note, that right now, snowfall parallel set up for each fragment of the pipeline is done separately. In production, I would consider only one initialization.

# === binaries ===
binExtractAnnotation <- "/Users/lalithaviswanathan/Documents/Project/PhyCA_Stats/code/extractAnnotation.pl"  # full path to extractAnnotation.pl script
binComUnderChars <- "/Users/lalithaviswanathan/Documents/Project/PhyCA_Stats/code/commentOutUnderscoreChars.sh" # shell script to comment out underscore characters in latex tabels
binExtHeader <- "/Users/lalithaviswanathan/Documents/Project/PhyCA_Stats/code/extractHeader.pl" # path to perl script extracting headers of data files (formatted for unifrac)
unifracBin <- "macqiime beta_diversity.py" # path to qiime pkg for unifrac, works with qiime svn revision 2507
# Note that ps2pdf is used for ps-to-pdf conversion
# and a2ping for conversion to other formats
# Right now it is assumed that ps2pdf and a2ping are in the PATH variable

# === script paths ===
coreLibPath <- "/Users/lalithaviswanathan/Documents/Project/PhyCA_Stats/code/coreLib.R"

# === input files ===
otuAnnFile <- "/Users/lalithaviswanathan/Documents/Project/Cubi_LChesnel_EPAN12_0120/abund_rep_set.GG98_FL.taxonomy.80.csv" # combo of previous seqTaxFile and seqAnnFile
pathToTree <- "/Users/lalithaviswanathan/Documents/Project/Cubi_LChesnel_EPAN12_0120/gg_97_otus_4feb2011.tre"

## project specific params:
# === input files ===
bFile <- "/Users/lalithaviswanathan/Documents/Project/Cubi_LChesnel_EPAN12_0120/ref97_gg98_filtered.txt" # FULL path to presence/absence table
aFile <- "/Users/lalithaviswanathan/Documents/Project/Cubi_LChesnel_EPAN12_0120/ref97_gg98_filtered.txt" # FULL path to abundance table
mFile <- "/Users/lalithaviswanathan/Documents/Project/Cubi_LChesnel_EPAN12_0120/metadata.csv"   # FULL path to metadata table; check the format of this file - get rid of OTU_ID field
# It is assumed that the second column of mFile, short_name, consists of short names of the files - no longer than 10 characters. Also, metadata column names must start with letters not numbers, else R puts an "X" in front.
unequal.scale.profiles <- T # profiles of abundance vs. sample for an OTU are scaled to that OTU's max if TRUE


## phyca-stats run params here to end

# R needs to know types of metadata columns. In particular, categorical variables have to to be identified by the user.
# There are two ways to specify categorical variables
# 1. write the names of categorical variabls into a comma delimeted file (mmFile)
# 2. specify the names of categorical variabls here (mFactors)
# mFactors <- NULL # set mFactors to NULL if the user does not specify categorical variables
mFactors <- c('CB315vCUB2');

# here define which factor(s) will be used for filtering.
mFilterFactors <- c('CB315vCUB2') # what is actually used in the anylsis only in DataRedn3 or 5 l 56 ff


# in doDataRedn6 the 'subjects' and 'pairings' columns must be specified here
# they must also be specified in the mFactors list
# currently the t.test is done but wilcox is also available as is FDR adjustment (paramters and conditional statements need to be added)
# subjectsFactor <- c('PatientID')   # s1, s2, etc
# pairingsFactor <- c('CvCtrl')  # pre, post etc (must be only two levels) # turn on for filter 6 and place in factor

# For bioenv, Mantel and Adonis tests need a selection of continuous
# variables.
# There are two options:

# 1. User specifies only categorical variables and all others are
#    treated as continuous (unless R detects a character variable and
#    excludes it from the list of continuous variables) This is a
#    good choice if the user wants to use bioenv() to identify
#    categorical variables associated with distances of environmental
#    variables.

# 2. User specifies continuous variables that she wants to use in
#    Mantel and Adonis tests.

# Note that contVar has to be defined, even if you don't want to specify any
# continuous variables. In this case, just set
# contVar <- NULL
contVar <- c('gDNAlab')

# === normalization paramters ===
normalize  <- T   # normalization flag, if set to TRUE, the normalization of all data reduction tables is performed.
normTarget <- 1e6     # target normalization value (used only for totalSum)
whichNorm  <- "totalSum"   #"totalSum", "rankNorm", "uqNorm"


# === flags for different data reduction methods ===
# Data reduction methods are numbered as in SOW document.
doDataRedn1 <- F
	Redn1SampleRemoval <- FALSE  # if a sample has zero OTUs in bt1 then remove sample
								#   may need to set to true if you want bt1 ordinations/clusterings if data set had a clean negative control 
								#  this can create problems with sample labels

doDataRedn2 <- F

doDataRedn3 <- F
	permute.filter3 <- TRUE 
	filter3.num.permutes <- 100

doDataRedn4 <- F
	rangeQtThld <- 0.9   # quantile treshold for the max/min ratios; only taxons with the max/min ratios in the 90% quantile will be selected in the data reduction method #4

doDataRedn5 <- TRUE 
	useANOVA <- TRUE
	useWelch <- F
	useKW    <- F
	pValThld <- .05   # threshold of ANOVA type p-values. pValThld <- FALSE will auto-select >= 50 taxa in filter 5
	filterMin <- 10 # if pValThld <- FALSE, phyca-stats will select >= this many OTUs, rounding the p value to the nearest greater power of 10 
	permute.filter5 <- 1    # 0=off, 1=shuffle samples into groups, 2=shuffle all abundandance data
	filter5.num.permutes <- 20
	
doDataRedn6 <- F
doDataRedn7 <- F
doDataRedn8 <- F
doDataRedn9 <- F
# 9 in one per correlated cluster per species/subfamily

# === parameters for distance methods ===
createSorensenDist    <- F
createUnifracDist     <- T
createWgtdUnifracDist <- T
createBrayCurtisDist  <- F
createEcldnDist       <- F

# === ordination parameters ===
doPCoA   <- T # Principal Coordinate Analysis
doPCA    <- FALSE # Principal Component Analysis
doNMDS   <- T # Non metric Multi-Dimensional Scaling
hexbinThld <- 330    # threshold for switching to a hexbin plot in ordination scatterplot plotOrd2dAllFtrs() routine

# === biplots ===
doEnvfit <- FALSE # used for biplots
contVar_Envfit <- 'nseqs'  # just pick one

ordTrellis1 <- FALSE # if TRUE, single categorical variable trellis plots are generated
ordTrellis2 <- FALSE # if TRUE, two categorical variables trellis plots are generated

# === hierarchical clustering parameters ===
doHCNN   <- FALSE # hierarchical clustering using single linkage
doHCAN   <- T # hierarchical clustering using average linkage
doHCFN   <- FALSE # hierarchical clustering using complete linkage
doHCWard <- FALSE # hierarchical clustering using Ward linkage
sepThickness <- 3;

# === singnificance testing parameters ===
# There are three significance type of tests, all from vegan package:
# bioenv(), mantel() and adonis(). Here are the flags that control
# which of them will be run.
doBioEnv <- FALSE
doMantel <- FALSE
doAdonis <- T
doAdonisPairs <- F


# params used for bioenv and mantel only
# The dissimilarity index used for community data in 'vegdist'
distIndex <- "bray" # possible values (specified in the appropriate
                    # pull down menu) are: manhattan, euclidean,
                    # canberra, bray, kulczynski, jaccard, gower,
                    # altGower, morisita, horn, mountford, raup ,
                    # binomial or chao.

# used for bioenv and mantel only
# The correlation method used in 'cor', 'bioenv' and 'mantel'subject of unknown paternal ethnic origin
corMethod <- "spearman" # possible values (specified in the
                        # appropriate pull down menu) : pearson,
                        # kendall, spearman


# === PAM ===
doPAM <- F # run PAM classifier for each categorical variable to classify samples with missing data within the selected variable
pamMinNoORFs <- 10
pamMaxNoORFs <- 20
pamfilter <- 1 # usually filter 9 is used (feb 2012)
## in PAM plots groups will be ordered alphabeticaly. So if you want to see "pre" before "post" use "1_pre" and "2_post"


# === iTOL ===
# itol requires filter 5, make sure doDataRedn5 <- TRUE
doITOL <- FALSE
itolBaseline <- 'healthy'
itolTaxonLevel <- 'family'


# === compositional, horizontal barcharts ===
doTaxonPropBarCharts <- F # generate bar charts of proportions at different taxonomic levels
doTaxonAbundancePropBarCharts <-TRUE
nTaxaInPropBarCharts <- 8 # maximum number of taxa to display in these bar charts.3
# there will be an extra black category for 'Other', for a total of
# nTaxaInPropBarCharts + 1 elements per bar
dumpBarChartTables <- TRUE

# === richness plots ===
doTaxonRichnessCharts <- F
dumpRichnessChartTables <- FALSE
richnessTaxonLevels <- c("phylum","class","order","family","genus","species")
# full options are c("phylum","class","order","family","genus","species")

### ========= end of input parameters section ===========## 
