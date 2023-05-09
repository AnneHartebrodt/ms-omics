setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')

renv::activate('r-single-cell-new')
renv::repair()
set.seed(11)


require(data.table)
library(GWENA)
library(magrittr)
require(stringr)
require(ComplexHeatmap)


bulk.data.dir<- '/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/bulk'

out.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/bulk/')
dir.create(out.dir)


annotation<-fread(file.path(bulk.data.dir, 'GSE138614_series_matrix.txt'), fill = TRUE)
annotation<-annotation[V1 %in% c('!Sample_source_name_ch1', '!Sample_characteristics_ch1', '')]
annotation<-t(annotation[,2:99])
annotation<-as.data.table(annotation)
colnames(annotation)<-c('V1','patient', 'diagnosis', 'tissue', 'lesion_type', 'individual')
annotation<-annotation[,c('patient', 'diagnosis', 'lesion_type')]


annotation$lesion_type<-as.character(sapply(annotation$lesion_type, 
                                            function(x) gsub('\\)', '', 
                                                             gsub('\\(', '', 
                                                                  str_extract_all(x, "\\(..\\)|\\(....\\)" )))))

data<-fread(file.path(bulk.data.dir, 'GSE138614_countMatrix.txt'))
rownames<-colnames(data[,2:ncol(data)])
columnnames<-rownames(data)
exp.mat<-t(data[, 2:ncol(data)])
colnames(exp.mat)<-columnnames
rownames(exp.mat)<-rownames


M1<-c("TGFBR3", "ARHGAP29", "F3", "S100A6", "ATP1A2", "LAMC1", "CFH", "PRELP", "KCNK2", "NID1", "CYP1B1", "FHL2", "RND3", "PLA2R1", "COBLL1", "TFPI", "IGFBP5", "KCNE4", "RBMS3", "PDZRN3", "AC092957.1", "TP63", "ANTXR2", "USP53", "AC097528.1", "ARHGAP10", "C7", "ITGA1", "NR2F1", "EBF1", "COL12A1", "LAMA4", "LAMA2", "COL1A2", "AC073114.1", "RBPMS", "GADD45G", "COL15A1", "LAMC3", "ITIH5", "VIM", "SVIL", "CXCL12", "BICC1", "PAPSS2", "PDLIM1", "IGFBP6", "DCN", "ADGRD1", "LHFPL6", "FLVCR2", "FLRT2", "ALDH1A2", "ISLR", "CEMIP", "LINC01197", "LINC00924", "NR2F2-AS1", "NR2F2", "CEP112", "ABCA9", "ABCA6", "COLEC12", "PLCB4", "SLC19A1", "TIMP3", "FBLN1", "AC002074.1", "FMO2", "LHCGR", "COL6A3", "TBX18", "AKR1C2", "SLC16A12", "COL3A1", "SLC6A20", "OLAH", "DCDC2C", "ADH1B", "AC093772.1")
M2<-c("TTC34", "C1orf87", "TCTEX1D1", "SPAG17", "SPATA17", "TOGARAM2", "DNAH6", "MAP3K19", "KIAA2012", "AC093865.1", "DAW1", "LTF", "ADAMTS9-AS1", "PLOD2", "ZBBX", "CCDC39", "CFAP299", "TTC29", "AC097518.2", "DNAH5", "ADGB", "FNDC1", "EGFR", "CD36", "AC083870.1", "AC069133.1", "C8orf37-AS1", "DNAI1", "ARMC3", "ARMC4", "RGR", "CFAP43", "CFAP46", "DCDC1", "CFAP54", "CCDC60", "LRRC9", "AL357093.2", "VWA3A", "DNAAF1", "TEKT1", "DNAH9", "CFAP47", "TMEM164", "AC020718.1", "AC105450.1", "FAM183A", "AC011632.1", "DYDC2", "BMP3")
M3<-c("KAZN", "IGSF21", "CSMD2", "GRIK3", "DAB1", "NEGR1", "SLC44A5", "PLPPR5", "COL11A1", "AL355306.2", "KCND3", "BCAN", "BX284613.2", "TNR", "ASTN1", "BRINP2", "LINC01344", "NMNAT2", "BRINP3", "KCNT2", "CRB1", "ESRRG", "SLC35F3", "RGS7", "PLD5", "KIF26B", "SMYD3", "PXDN", "AC064875.1", "AC068286.2", "ALK", "NRXN1", "LRRTM4", "LINC01965", "ST6GAL2", "CNTNAP5", "GPR17", "TMEM163", "LRP1B", "GALNT13", "SCN3A", "SCN2A", "SCN1A", "SCN9A", "GAD1", "RAPGEF4", "ZNF804A", "MAP2", "UNC80", "ERBB4", "AC068051.1", "VWC2L", "AC007563.2", "SCG2", "SPHKAP", "PID1", "CHL1", "CNTN6", "CNTN4", "LRRN1", "GRM7", "SATB1-AS1", "THRB", "ARPP21", "ERC2", "AC126121.3", "TAFA1", "CNTN3", "CADM2", "GRAMD1C", "GAP43", "LSAMP", "STXBP5L", "SEMA5B", "TMEM108", "LINC01322", "NLGN1", "KCNMB2", "KCNMB2-AS1", "AC068308.1", "IGF2BP2", "DGKG", "IL1RAP", "FGF12", "PLAAT1", "STK32B", "KCNIP4", "PCDH7", "PDGFRA", "LINC02283", "ADGRL3", "EPHA5", "PRKG2", "FAM13A", "CCSER1", "GRID2", "DKK2", "NDST4", "RNF150", "AC107223.1", "PDGFC", "GRIA2", "TLL1", "SEMA5A", "TRIO", "BASP1", "BASP1-AS1", "LINC02223", "CDH18", "AC091946.1", "CDH10", "PDZD2", "RGS7BP", "AC093523.1", "PDE8B", "VCAN", "AC113383.1", "MCTP1", "FBN2", "STK32A", "SOX4", "AL445250.1", "COL9A1", "RIMS1", "AL391840.1", "SNAP91", "GRIK2", "SLC16A10", "SLC35F1", "AL133346.1", "UST", "AKAP12", "RGS17", "MMD2", "NXPH1", "AC004852.2", "AC004160.1", "AC011287.1", "ETV1", "DGKB", "TMEM196", "BMPER", "AMPH", "POU6F2", "INHBA", "AC005537.1", "HECW1", "VWC2", "CALN1", "SEMA3E", "CDK14", "AC079760.1", "LHFPL3", "NRCAM", "AC002066.1", "KCND2", "PTPRZ1", "AC006148.1", "PLXNA4", "TMEM178B", "DPP6", "PTPRN2", "CSMD1", "MSRA", "XKR6", "AF279873.3", "FGFR1", "SNTG1", "XKR4", "FAM110B", "TOX", "CLVS1", "MIR2052HG", "RALYL", "AF121898.1", "MMP16", "AC090578.1", "NCALD", "BAALC", "RIMS2", "ZFPM2", "CSMD3", "TNFRSF11B", "KHDRBS3", "LINC02055", "GLDC", "MLLT3", "TEK", "SHC3", "PLPPR1", "BRINP1", "ADARB2", "ITGA8", "KIAA1217", "GPR158", "ZEB1-AS1", "ZEB1", "PCDH15", "KCNMA1", "SLIT1", "SORCS3", "SORCS1", "ATRNL1", "GFRA1", "PLPP4", "AC023282.1", "SOX6", "NELL1", "ANO5", "LUZP2", "MPPED2", "SLC1A2", "LRRC4C", "NEAT1", "GRM5", "GRIA4", "GUCY1A2", "PKNOX2", "OPCML", "IGSF9B", "AC005906.2", "PTPRO", "SLC2A13", "LRRK2", "PDZRN4", "PRICKLE1", "TAFA2", "BEST3", "MYRFL", "TRHDE", "AC078923.1", "NTN4", "NOS1", "KSR2", "TMEM132B", "TMEM132C", "TMEM132D", "FRY", "FREM2", "KLHL1", "KLF12", "HS6ST3", "FGF14", "FAM155A", "MYO16", "AL110292.1", "AL392023.2", "LRFN5", "MDGA2", "TRIM9", "SAMD4A", "SYT16", "KCNH5", "SMOC1", "SLC8A3", "MEG3", "MEG8", "GABRB3", "MEGF11", "CTXND1", "NTRK3", "ABHD2", "ADAMTS17", "SHISA9", "U91319.1", "XYLT1", "GSG1L", "AC106793.1", "AC106729.1", "CDH13", "CRISPLD2", "NTN1", "ASIC2", "HAP1", "CA10", "CACNG5", "CACNG4", "LINC01483", "LINC00511", "DLGAP1", "ARHGAP28", "CABLES1", "CHST9", "NOL4", "RIT2", "DCC", "RAB27B", "CDH7", "AC114689.3", "DOK6", "NETO1", "GALR1", "AC124254.2", "CACNA1A", "FERMT1", "PLCB1", "SNAP25", "MACROD2", "SLC24A3", "CST3", "SULF2", "COL20A1", "MYT1", "GRIK1", "DSCAM", "SEZ6L", "MIAT", "CACNG2", "NLGN4X", "CNKSR2", "AC112493.1", "TSPAN7", "OPHN1", "NEXMIF", "DACH2", "PCDH11X", "PAK3", "KLHL13", "GRIA3", "TENM1", "HS6ST2", "AFF2", "NLGN4Y", "TTTY14", "AC068643.1", "AC107419.1")
M4<-c("PIK3R3", "ADGRL4", "ADGRL2", "PALMD", "ACKR1", "HMCN1", "PLA2G4A", "NR5A2", "AL929091.1", "SNTG2", "CRIM1", "EPAS1", "IL1R1", "NOSTRIN", "ITGA6", "CAVIN2", "EPHA4", "SGPP2", "OSBPL10", "CMTM8", "PTPRG", "PLA1A", "TM4SF1", "WWTR1", "MECOM", "TFRC", "LDB2", "AC114757.1", "ABCG2", "EMCN", "LEF1", "VEGFC", "SORBS2", "ARL15", "ZNF366", "LINC02147", "NEDD9", "EDN1", "RNF144B", "ADGRF5", "ADGRG6", "PLEKHG1", "PDE10A", "THSD7A", "IGFBP3", "ABCB1", "RUNDC3B", "NAMPT", "CADPS2", "RUNX1T1", "SLC1A1", "EGFL7", "ST8SIA6", "NRP1", "GALNT18", "FLI1", "ANO2", "VWF", "A2M", "CPNE8", "IRAK3", "PTPRB", "FLT1", "SLC7A1", "DACH1", "TBC1D4", "PRKCH", "HIF1A-AS3", "ATP10A", "CGNL1", "SMAD6", "THSD4", "MRTFB", "IL4R", "SLC7A5", "NXN", "PECAM1", "PTPRM", "LAMA3", "DOCK6", "ID1", "TGM2", "TSHZ2", "AF130417.1", "CYYR1", "ERG", "ETS2", "BACE2", "CLDN5", "ITM2A", "PGR-AS1", "SOX18", "IL1RL1", "TPO")
M5<-c("HES4", "ALPL", "ROR1", "COL24A1", "KIRREL1", "SLC30A10", "KLHL29", "CYTOR", "MIR4435-2HG", "TEX41", "RBMS1", "AC019197.1", "SLC38A11", "SCN7A", "LINC01117", "MYO1B", "FN1", "COL4A4", "RFTN1", "PTH1R", "EPHA3", "PHLDB2", "P2RY14", "LINC00504", "AC093791.2", "IGFBP7", "UBA6-AS1", "ANXA3", "ENPEP", "LINC01060", "AC008825.1", "CDH6", "PRLR", "GHR", "STARD4-AS1", "MCC", "EGR1", "PDGFRB", "SLIT3", "BTNL9", "CASC15", "CCN2", "AL138828.1", "UTRN", "SAMD5", "AL033504.1", "SMOC2", "GPR85", "GRM8", "LZTS1", "CEBPD", "MSC-AS1", "HEY1", "SNTB1", "ANXA1", "FRMD3", "CENPP", "NR4A3", "LINC02664", "VCL", "ACTA2", "RBM20", "GRK5", "RRAS2", "ANO3", "AC021723.2", "C11orf96", "PRSS23", "NOX4", "TRPC6", "CACNA1C", "MGP", "RERG", "PDE3A", "TMTC1", "ANO6", "NDUFA4L2", "AC027288.3", "TESC", "TRPC4", "COL4A1", "COL4A2", "LINC02315", "FRMD6-AS2", "SYNE2", "ACTN1", "PGF", "AC104574.2", "AC012409.2", "CX3CL1", "RASD1", "PLXDC1", "GJC1", "FAM20A", "LINC01482", "AC105094.2", "NOTCH3", "KLF2", "JAG1", "MYL9", "NHS", "HTR2C", "CARMN", "CD34", "ADAMTS15", "SLC38A4", "GJA4", "TWIST2", "TBX2", "MXRA5", "CYSLTR2")
M6<-c("PTGER3", "LINC01725", "RALGPS2", "RHEX", "KCNS3", "IGKC", "ACOXL", "LYPD6B", "AC078883.1", "AC009315.1", "EPHB1", "FNDC3B", "LINC01208", "AC092546.1", "SEL1L3", "RASGEF1B", "BANK1", "AC022325.2", "AC078881.1", "NREP", "BMP6", "AL591518.1", "PRDM1", "AL589693.1", "CNBD1", "LINC02237", "SLC30A8", "PIP5K1B", "RTKN2", "CHST15", "AC087280.2", "ZNF215", "ABTB2", "LINC02745", "IFNG-AS1", "AL442636.1", "DCT", "IGHG4", "IGHGP", "IGHG1", "IGHG3", "IGHM", "AC068875.1", "ACAN", "MCTP2", "HS3ST3B1", "CHODL", "IGLC1", "IGLC2", "IGLC3", "XBP1", "AC002072.1", "IRF4", "LINC00958", "LINC02739", "TENT5C", "AC005042.2", "LINC00607", "JCHAIN", "MZB1", "AC008080.4", "AP000676.5", "AC006206.2", "IGHA1", "IGHV1-69D", "CD79A", "IGLV1-51", "SLAMF7", "FCRL5", "RASSF6", "IGLV2-18", "FCRL2", "AL353753.1", "FAM30A")
M7<-c("SLC6A17", "LINC00970", "C1orf115", "GREM2", "VSNL1", "AC110614.1", "AC073091.3", "KIAA1211L", "AC013265.1", "ZNF385B", "LINC01821", "PTH2R", "LINC01811", "LRTM1", "AC093607.1", "NWD2", "LINC01378", "FAM160A1", "RXFP1", "GLRA3", "LINC02112", "FSTL4", "AC008415.1", "LY86-AS1", "MLIP", "COL19A1", "CNR1", "AL591519.1", "MCHR2", "GRM1", "MPP6", "AC021613.1", "NEFL", "LINC01414", "SULF1", "AC060765.2", "AC004083.1", "ADAMTSL1", "AL391117.1", "LINC02645", "LINC01435", "AL136119.1", "LINC02755", "SLC22A10", "NRGN", "TAC3", "LINC02822", "LINC00507", "SIAH3", "AL356295.1", "LINC02306", "CHRM5", "ANKRD34C-AS1", "LINC02254", "HS3ST2", "CALB2", "PPM1E", "EFCAB3", "NPTX1", "AP000829.1", "DLGAP1-AS4", "MIR924HG", "CBLN2", "SLC17A7", "NNAT", "AL117329.1", "PIK3C2G", "LINC01776", "LINC02263", "AL096799.1", "SYT10", "AL445213.2", "AP000146.1")
M8<-c("CDC14A", "RGS16", "CHIT1", "AC018742.1", "AC007100.1", "IL1B", "LINC01266", "KCTD8", "CAMK2D", "MEF2C-AS1", "AC104123.1", "AC008591.1", "LINC01933", "LINC02542", "TAGAP", "GPNMB", "AC011586.2", "AC016074.2", "P4HA1", "DISC1FP1", "ARHGAP20", "LINC02698", "AP001977.1", "CD163L1", "FAR2", "AC090630.1", "EPSTI1", "AL049828.1", "LINC01500", "LINC02328", "LINC01146", "BCL2A1", "AC104041.1", "HS3ST3A1", "AC034268.2", "BRIP1", "AC025887.2", "NEDD4L", "MAFB", "AC079142.1", "SAT1", "MIR222HG", "AC011751.1", "AC005699.1", "CYP19A1", "IL21R", "RAB42", "CHRNA1", "ZNF385D-AS1", "AC087273.2", "HTRA4", "AC044893.1", "AL355838.1", "AQP9", "IL1A", "AC064834.1", "LINC02057", "AC006974.2", "AC090796.1")
M9<-c("GRHL3", "AL583808.1", "PAPPA2", "AL390957.1", "RGS2", "PROX1-AS1", "DNAH14", "WNT9A", "AL390860.1", "LINC00276", "AC079352.1", "TMEM178A", "AC012494.1", "AC007364.1", "AC023469.1", "LINC01876", "PDE1A", "AC007319.1", "C2CD6", "CPS1", "AC034195.1", "MOBP", "SAMMSON", "AC068633.1", "GPR149", "SST", "AC110809.1", "COL25A1", "AC093765.2", "AC115622.1", "LINC02511", "AC109927.1", "LINC01098", "LINC01099", "LINC02119", "AC008945.2", "AC112206.2", "ENC1", "LUCAT1", "LINC02163", "LINC01170", "AC008571.2", "RANBP17", "AL033523.1", "AL133405.2", "ZSCAN31", "TSBP1-AS1", "EYS", "HS3ST5", "AL590550.1", "RNF217-AS1", "AL591115.1", "STXBP5-AS1", "SOD2", "AC004949.1", "AGMO", "LINC01445", "SERPINE1", "AC009264.1", "AC246817.1", "LINC00529", "SLC39A14", "AC009623.1", "AC103770.1", "LINC01608", "LINC00536", "AC037486.1", "LINC01151", "PCAT1", "LINC01505", "SVEP1", "RNLS", "LIPN", "OPALIN", "TPH1", "GDPD4", "CARD18", "CRYAB", "AC090015.1", "AC016152.1", "LINC02343", "AL354809.1", "AL136441.1", "AL590807.1", "AL355916.2", "AC008050.1", "NRXN3", "PWRN1", "SCNN1G", "SLC5A11", "CNGB1", "CNTNAP4", "LINC02073", "GREB1L", "MIR4527HG", "MYO5B", "LINC01416", "LINC01924", "CNDP1", "GADD45B", "JUNB", "AL050403.2", "LBP", "AJ009632.2", "AP000477.3", "LINC01695", "LINC00313", "PPEF1", "AC136489.1", "LINC01204", "AC233296.1", "SPANXA2-OT1", "AC068759.1", "AC022433.1", "SLCO6A1", "GLYAT", "LINC02335", "AL354994.1", "AC004870.4", "LINC00616", "AC093912.1", "AC015912.3", "TPH2")
M10<-c("CHD5", "SFN", "AGBL4", "AC119674.1", "AL513166.1", "ST6GALNAC5", "CCN1", "LINC02607", "MIR137HG", "AL445433.2", "OLFM3", "AL136114.1", "CACNA1E", "AL136456.1", "KCNH1", "RYR2", "CHRM3", "KMO", "MYT1L", "LINC01250", "GALNT14", "AC009975.1", "LINC01122", "BCL11A", "AC007389.1", "ANKRD30BL", "THSD7B", "NXPH2", "LYPD6", "KCNJ3", "SLC4A10", "KCNH7", "SNHG31", "MARCH4", "DIRC3", "AC019211.1", "NYAP2", "AC019068.1", "SYN2", "SGO1-AS1", "ZNF385D", "CCK", "CADPS", "SYNPR", "SUCLG2-AS1", "ROBO2", "EPHA6", "CLSTN2", "AC117386.2", "AC073365.1", "AC098829.1", "SLIT2", "GABRA4", "KIT", "KIAA1211", "CXCL2", "FRAS1", "FSTL5", "GALNTL6", "AC108169.1", "TENM3", "CDH12", "CDH9", "LINC02224", "HCN1", "RAB3C", "SV2C", "KIAA0825", "LIX1-AS1", "LINC01340", "ADAMTS19", "CHSY3", "CXCL14", "KCTD16", "AC132803.1", "GRIA1", "PTTG1", "GABRB2", "GABRA1", "GABRG2", "TENM2", "FAM153CP", "PHACTR1", "VEGFA", "PTCHD4", "HCRTR2", "KHDRBS2", "KCNQ5", "HTR1E", "EPHA7", "EYA4", "AL031056.1", "AC073332.1", "AC003044.1", "ABCA13", "GALNT17", "PCLO", "ZNF804B", "DLX6-AS1", "RELN", "POT1-AS1", "CHRM2", "CNTNAP2", "DLGAP2", "SGCZ", "NRG1", "UNC5D", "ZMAT4", "KCNB2", "STMN2", "GEM", "AC104248.1", "ARC", "FREM1", "SH3GL2", "AL445623.2", "ELAVL2", "LINGO2", "GABBR2", "GRIN3A", "GRIN1", "CACNA1B", "CACNB2", "GAD2", "AC013287.1", "CH25H", "TCERG1L", "NRIP3", "KIAA1549L", "CNTN5", "PDGFD", "EMP1", "GRIN2B", "AC024901.1", "NELL2", "GRIP1", "PTPRR", "KCNC2", "SYT1", "AC138123.1", "IGF1", "BTBD11", "RPH3A", "SRRM4", "AC007368.1", "RIMBP2", "ATP8A2", "MTUS2", "AL445255.1", "SLC35F4", "AC004817.5", "SNHG14", "GABRG3", "LINC02853", "UNC13C", "AGBL1", "LINC00923", "RBFOX1", "GRIN2A", "CACNG3", "MYLK3", "CDH8", "CFAP52", "SHISA6", "RBFOX3", "CELF4", "AC009899.1", "CCBE1", "AC010624.5", "PAK5", "PCSK2", "PTPRT", "KCNJ6", "SYN3", "CD99", "FRMPD4", "AL591501.1", "SYTL5", "IL1RAPL2", "FGF13", "AC004946.1", "TENT5B", "AC007920.1", "Z97205.2", "AC073968.2", "AL731571.1", "LINC02351", "BEND4", "AL445207.1", "AL358335.2", "LINC02009", "AC020651.2")
M11<-c("PRDM16", "ID3", "MAN1C1", "AL603840.1", "JUN", "LINC01748", "AC099792.1", "DOCK7", "CACHD1", "GNG12-AS1", "WLS", "ZRANB2-AS2", "DDAH1", "AL627316.1", "NTNG1", "AHCYL1", "CHI3L2", "KCNN3", "DDR2", "F5", "PRRX1", "NPL", "C1orf21", "TGFB2", "FMN2", "LINC00299", "ID2", "LTBP1", "EML6", "EFEMP1", "CCDC85A", "ANTXR1", "CTNNA2", "TCF7L1", "NPAS2", "DPP10", "ARHGEF4", "AC016766.1", "AC073050.1", "CERS6", "OSBPL6", "TTN", "AC092640.1", "PPP1R1C", "DNAH7", "PARD3B", "IQCA1", "ACKR3", "SNED1", "CHL1-AS2", "AC077690.1", "AC092422.1", "RARB", "LRRC3B", "GASK1A", "CACNA2D3", "FAM107A", "ADAMTS9", "ADAMTS9-AS2", "LRIG1", "ABI3BP", "BOC", "AC026341.1", "AC092691.1", "KALRN", "ALDH1L1", "PLSCR4", "CP", "MED12L", "LINC02006", "ARHGEF26-AS1", "ARHGEF26", "PLCH1", "SERPINI2", "WDR49", "TNIK", "MASP1", "BCL6", "ATP13A4", "HES1", "LINC01182", "CD38", "AC097515.1", "ADGRA3", "LINC02506", "GABRG1", "GABRA2", "GABRB1", "SLC4A4", "ADAMTS3", "SHROOM3", "LINC01094", "LINC01088", "MAPK10", "SPARCL1", "BMPR1B", "NPNT", "ETNPPL", "SYNPO2", "AC073475.1", "FGF2", "SLC7A11", "AC002460.2", "IQCM", "DCLK2", "AC079298.3", "DCHS2", "LRAT", "AC093817.2", "AC084740.1", "CPE", "PALLD", "GPM6A", "PDLIM3", "F11-AS1", "PLEKHG4B", "ADCY2", "CTNND2", "ADAMTS12", "RANBP3L", "OSMR-AS1", "MAST4", "AC114971.1", "AC091826.2", "ADGRV1", "AC124854.1", "NR2F1-AS1", "EFNA5", "TMEM232", "KCNN2", "PRR16", "SNCAIP", "GRAMD2B", "MARCH3", "MEIKIN", "PPP2R2B", "DPYSL3", "AC113414.1", "WWC1", "HULC", "GMPR", "AL022068.1", "ID4", "CARMIL1", "COL21A1", "AL606923.2", "B3GAT2", "OGFRL1", "CD109", "AL355612.1", "FUT9", "UFL1-AS1", "AL589740.1", "NR2E1", "FAM184A", "AL096854.1", "GJA1", "TRDN", "TPD52L1", "MOXD1", "LINC00271", "PDE7B", "EZR", "AL078604.4", "PACRG", "AL078602.1", "DNAH11", "AQP1", "AC083864.5", "SUGCT", "GLI3", "AEBP1", "TPST1", "SEMA3A", "SEMA3D", "AC002069.2", "AC002429.2", "CDHR3", "HILPDA", "MFHAS1", "SLC7A2", "AC087854.1", "ADRA1A", "CLU", "RGS20", "CA8", "AC023095.1", "NKAIN3", "C8orf34", "EYA1", "CRISPLD1", "PKIA-AS1", "PKIA", "AC083837.1", "SPAG1", "RNF19A", "ANGPT1", "ADCY8", "FAM135B", "GLIS3", "IL33", "CNTNAP3B", "FAM189A2", "TRPM3", "TMC1", "ALDH1A1", "RORB", "PCSK5", "GNA14", "AL157886.1", "NTRK2", "DAPK1", "WNK2", "ABCA1", "DEC1", "TNC", "AL160272.1", "PFKP", "LINC02649", "AL392086.3", "SLC39A12", "MALRD1", "NEBL", "LINC00836", "PARD3", "FRMPD2", "PRKG1", "SLC16A9", "LINC01515", "CAMK2G", "NRG3", "PPP1R3C", "LGI1", "PLCE1", "SORBS1", "HPSE2", "ADD3-AS1", "ABLIM1", "AL137025.1", "CASC2", "IFITM3", "STK33", "ADM", "MRVI1", "SPON1", "BBOX1", "AL358944.1", "PAX6", "PAUPAR", "CD44", "PAMR1", "TENM4", "MIR4300HG", "FAT3", "ARHGAP42", "YAP1", "MGST1", "LMO3", "RERGL", "PLEKHA5", "SLCO1C1", "ABCC9", "ST8SIA1", "AC053513.1", "SOX5", "LMNTD1", "OVCH1-AS1", "CNTN1", "SLC38A1", "AC008014.1", "AC068305.2", "NAV3", "ACSS3", "LRRIQ1", "MGAT4C", "KITLG", "AC009522.1", "RMST", "AC079385.1", "RFX4", "SDS", "HSPB8", "DCLK1", "LINC00598", "LINC00378", "OBI1-AS1", "GPC5", "GPC6", "FARP1", "AL137139.2", "ITGBL1", "AL390755.1", "AL135878.1", "PRKD1", "EGLN3", "LINC00609", "AL132857.1", "DAAM1", "AL359232.1", "RGS6", "STON2", "SLC24A4", "SERPINA3", "RYR3", "AC013652.1", "AC012405.1", "AC073941.1", "RORA", "CA12", "ACSBG1", "ADAMTSL3", "RGMA", "AC091078.1", "AC110023.1", "ARRDC4", "LINC02251", "MT2A", "MT1E", "MT1M", "MT1G", "MT1X", "GINS3", "AC092131.1", "AC105411.1", "ATP1B2", "SLC47A2", "CCL2", "GFAP", "ANKFN1", "PITPNC1", "CDC42EP4", "LAMA1", "AQP4-AS1", "AQP4", "TTR", "DTNA", "FHOD3", "LINC00907", "SLC14A1", "ZBTB7C", "MAPK4", "AC022031.2", "AC006305.1", "ANGPTL4", "COL5A3", "CYP4F12", "CPAMD8", "NCAN", "ZNF98", "APOE", "HIF3A", "AL109809.5", "AL121757.2", "EYA2", "DOK5", "CDH4", "MIR99AHG", "SLC25A18", "ZNRF3", "GYG2", "MID1", "AL807742.1", "PTCHD1-AS", "MAOB", "XIST", "SYTL4", "GPC4", "AC114964.2", "LINC02327", "LINC02141", "RASSF9")
M12<-c("VAV3", "AC096637.1", "MIR3681HG", "AC097512.1", "AC097480.1", "LINC00499", "PURPL", "AC079465.1", "NPY", "VGF", "EGR3", "AC103796.1", "GPRC5A", "HMGA2", "OTX2-AS1", "LINC01727", "THBD", "CITED1", "LINC02046", "AC069410.1", "AC020584.1", "LINC02109", "LINC00395", "CCL3", "KIF18B", "HES5", "DTL", "NCAPG", "TXK", "HIST1H2BH", "AC002454.1", "PRF1", "MKI67", "KCNJ12", "CCL3L1", "TOP2A", "AC069277.1", "LINC00698", "AC106771.1")
M13<-c("IFI44L", "CD247", "LINC01934", "ITGA4", "STAT4", "HRH1", "ARHGEF3", "CD96", "DTHD1", "RHOH", "INPP4B", "EMB", "PARP8", "IQGAP2", "CAMK4", "TNFAIP8", "F13A1", "RIPOR2", "THEMIS", "SAMD3", "AL035446.2", "SYTL3", "MDFIC", "AC100849.1", "LINC01242", "MRC1", "MS4A4A", "CD163", "CD69", "CLEC2B", "TESPA1", "LINC01619", "DIAPH3", "AL162493.1", "GNG2", "TC2N", "LGMN", "BCL11B", "CCL5", "IKZF3", "SKAP1", "TTC39C", "SIGLEC1", "GRAP2", "DPP4", "AC011029.1", "CD2", "IL7R", "SCML4", "LINC00944", "TRAC", "SLFN12L", "CCL4", "LINC00470", "PYHIN1", "AC006369.1", "AC022217.3", "CALCR", "AL136962.1", "LILRB5", "FAM9B", "GPR174")
M14<-c("SLC2A5", "TNFRSF1B", "PDPN", "TMEM51", "LINC01141", "C1QC", "C1QB", "CSF3R", "SMAP2", "HIVEP3", "CD53", "LINC02798", "SRGAP2B", "S100A8", "RCSD1", "MPZL1", "RASAL2", "NIBAN1", "RGS1", "PTPRC", "SRGAP2", "AL392172.2", "DISP1", "DISC1", "SLC8A1", "ZFP36L2", "PLEK", "ARHGAP25", "MGAT4A", "NCK2", "SH3RF3", "MERTK", "KYNU", "ARHGAP15", "CACNB4", "SLC11A1", "SP100", "INPP5D", "PPARG", "PLCL2", "TGFBR2", "CMTM7", "CX3CR1", "STAB1", "CACNA1D", "HCLS1", "CD86", "P2RY12", "LINC00578", "GNB4", "KLHL6", "ST6GAL1", "TPRG1", "AC006230.1", "TMEM156", "RBM47", "LINC02232", "BMP2K", "ARHGAP24", "AC131944.1", "SPP1", "GPRIN3", "ALPK1", "LINC01091", "MAML3", "TLR2", "GUCY1A1", "FAM149A", "OTULINL", "C5orf17", "AC114930.1", "SLC1A3", "FYB1", "S100Z", "LHFPL2", "MEF2C", "ST8SIA4", "CD14", "ARHGAP26", "CSF1R", "CD74", "AC008691.1", "DOCK2", "LCP2", "KCNIP1", "RASGEF1C", "LY86", "RREB1", "GCNT2", "CD83", "HLA-DRA", "HLA-DRB1", "FGD2", "TRERF1", "PLA2G7", "AL357522.1", "MAN1A1", "ARHGAP18", "MAP3K5", "NHSL1", "MTHFD1L", "IPCEF1", "RNASET2", "CARD11", "SDK1", "SCIN", "HDAC9", "OSBPL3", "SKAP2", "CPVL", "WIPF3", "AOAH", "TNS3", "IKZF1", "DOCK4", "FOXP2", "TFEC", "CPED1", "TBXAS1", "PRKAG2", "CTSB", "AC068587.4", "MSR1", "FGL1", "CSGALNACT1", "ADAM28", "AC120193.1", "LYN", "CPA6", "PAG1", "ZFPM2-AS1", "OXR1", "CCDC26", "KCNQ3", "DENND3", "DOCK8", "BNC2", "SYK", "ROR2", "TGFBR1", "AL162414.1", "SFMBT2", "CELF2", "CAMK1D", "FRMD4A", "AC044781.1", "PLXDC2", "APBB1IP", "ALOX5", "WDFY4", "LNCAROD", "SRGN", "PALD1", "LRMDA", "CERNA2", "LINC01374", "MYOF", "ENTPD1", "BLNK", "PIK3AP1", "KCNQ1", "PDE3B", "AP001636.3", "MS4A6A", "MS4A7", "SLCO2B1", "CTSC", "MAML2", "IL18", "SORL1", "UBASH3B", "GRAMD1B", "LINC02712", "OLR1", "ETV6", "LRMP", "ITPR2", "PCED1B", "NCKAP1L", "RASSF3", "AC025569.1", "LINC02391", "CHST11", "LRCH1", "LPAR6", "DLEU1", "ABCC4", "GPR183", "AL163541.1", "ZFP36L1", "FOXN3", "KCNK13", "RIN3", "CYFIP1", "FMN1", "ATP8B4", "GLDN", "AKAP13", "MEF2A", "LRRK1", "ATF7IP2", "CIITA", "SYT17", "HS3ST4", "ITGAX", "LPCAT2", "ZFHX3", "IRF8", "ABR", "PIK3R5", "RHBDF2", "RAB31", "LDLRAD4", "RTTN", "SOCS6", "C3", "VAV1", "MYO1F", "HAMP", "FCGBP", "POU2F2", "APOC1", "APOC2", "LILRB1", "SYNDIG1", "HCK", "AL139351.3", "NFATC2", "MIR646HG", "SAMSN1", "LINC01684", "BACH1", "RUNX1", "MX2", "TRPM2", "HMOX1", "PDGFB", "CSF2RA", "ANOS1", "ARHGAP6", "CYBB", "VSIG4", "EDA", "DIAPH2", "XACT", "LINC00278", "AC078845.1", "LINC02642", "AC021504.1")
M15<-c("PKN2-AS1", "AC114485.1", "AL445426.1", "MAGI3", "CHI3L1", "AL356108.1", "SLC8A1-AS1", "AC007402.1", "AFF3", "AC011246.1", "SCHLAP1", "AC021851.2", "CALCRL", "BHLHE40", "NEK10", "PRICKLE2", "LINC02008", "COL8A1", "ARHGAP31", "ZFYVE28", "PROM1", "AC024230.1", "AC093895.1", "CXXC4-AS1", "AC093765.3", "ANKRD33B", "PLCXD3", "PDE4D", "ADAMTS6", "PAM", "LINC02208", "AC008568.1", "PRELID2", "AC109466.1", "AL353138.1", "TENT5A", "AL356124.1", "SYNE1", "HSPB1", "TAC1", "ANKRD7", "STC1", "PLAT", "BAALC-AS1", "SAMD12-AS1", "SUSD1", "AFAP1L2", "BAG3", "LINC02550", "SLC2A3", "WIF1", "LINC02458", "AC002070.1", "LINC00351", "NALCN-AS1", "AL139023.1", "SLC25A21", "TTC6", "AL365295.1", "RHOJ", "FOS", "DIO2", "SMAD3", "AC104035.1", "AC025810.1", "ADAMTS18", "AC090403.1", "SETBP1", "AC091576.1", "PCAT19", "SLCO4A1", "STS", "AC073529.1", "AL050309.1", "AL008633.1", "PCDH11Y", "LINC01324", "AC097467.3", "LINC02464", "MLPH", "LINC02211", "METTL7B", "LINC01708", "MMRN1", "DSP", "TNFRSF10C", "PKHD1L1", "LINC02185", "SLC52A3", "RNF17", "AC012078.2", "LAX1", "PTGS2", "AC110611.2", "FLT4", "ARHGEF15", "AP001037.1")


modnames<-c('PC_MS',
            'AS_MS',
            'OPC',
            'Endo_MS',
            'PC_MS',
            'B_MS',
            'N_MS',
            'IM_MS',
            'OLIGO',
            'N_MS',
            'AS',
            #'AS_IM_MS',
            'B_MS',
            'IM',
            'Endo_MS'
)
modules<-list(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, 
              #M12, 
              M13, M14, M15)

modnames<-paste0('M', c(1:11, 13:15), '_', modnames)

library(biomaRt)
data$V1<-gsub('\\.[0-9]', '', data$V1)
mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")
z <- getBM(c("ensembl_gene_id","hgnc_symbol"), "ensembl_gene_id", data$V1, mart)




data<- merge(data, z, by.x = 'V1', by.y = 'ensembl_gene_id')

for(m in 1:length(modules)){

sub<- data[hgnc_symbol %in% modules[[m]]]

mat<-as.matrix(sub[,2:(ncol(data)-1)])
rownames(mat)<-sub$hgnc_symbol
mat<-scale(mat)

lesion_color  = c('NAWM' = '#00798c', 'AL'='#d1495b', 'CA'='#edae49', 'WM'='#66a182', 'RL'='#2e4057', 'IL' = '#c02ad1')

column_ha = HeatmapAnnotation(diagnosis = annotation$diagnosis, lesion_type = annotation$lesion_type, 
                              col = list(diagnosis = c('diagnosis: Multiple sclerosis'='#3D4849', 'diagnosis: Control'='#D7D9D9'), lesion_type = lesion_color),
                              simple_anno_size = unit(1.5, "cm")
)

k= 4
hm<-Heatmap(mat, show_column_names  = T, 
            column_km  = k, 
            top_annotation = column_ha,
            show_row_names = T,
            heatmap_legend_param = list(legend_height = unit(4, "cm"), legend_width = unit(4, "cm")))
pdf(file = file.path(out.dir, paste0(modnames[m], '.pdf')), width = 30, height = 20)
ht = draw(hm)
dev.off()

}



M6_sub<- c('IGHG1', 'ABTB2', 'IGCK', 'TENT5C')
sub<- data[hgnc_symbol %in% M6_sub]


mat<-as.matrix(sub[,2:(ncol(data)-1)])
rownames(mat)<-sub$hgnc_symbol
mat<-scale(mat)


lesion_color  = c('NAWM' = '#00798c', 'AL'='#d1495b', 'CA'='#edae49', 'WM'='#66a182', 'RL'='#2e4057', 'IL' = '#c02ad1')

annotation$ms<-ifelse(annotation$diagnosis=='diagnosis: Control', 'C', 'MS')
column_ha = HeatmapAnnotation(Diagnosis = annotation$ms, 
                              `Lesion type` = annotation$lesion_type, 
                              col = list(Diagnosis = c('MS'='#3D4849', 'C'='#D7D9D9'), `Lesion type` = lesion_color),
                              annotation_legend_param = list(
                                Diagnosis = list(direction = "horizontal", nrow=1),
                                `Lesion type` = list( direction='horizontal', nrow=2),
                                labels_gp = gpar(fontsize=12)),
                              annotation_name_gp = gpar( fontsize = 18),
                              simple_anno_size = unit(1, "cm")
)

k= 3
hm<-Heatmap(mat, 
            show_column_names  = T, 
            column_labels = annotation$patient,
            column_km  = k, 
            column_title = c('', '', ''),
            top_annotation = column_ha,
            show_row_names = T,
            show_column_dend = F,
            show_row_dend = F,
            name = 'Expression Level',
            row_names_gp = gpar(fontsize = 18),
            raster_by_magick = TRUE,
            raster_quality  = 5,
            heatmap_legend_param = list(legend_height = unit(2, "cm"), 
                                        legend_width = unit(4, "cm"),
                                        direction='vertical',
                                        labels_gp = gpar(fontsize = 10)
                                        ))
pdf(file = file.path(out.dir, paste0('M6_ABTB2_IGHG1_TENT5.pdf')), width = 18, height = 2.5)

ht = draw(hm, merge_legend = TRUE)
dev.off()



############

M13_sub<-c('GNG2', 'CLEC2B', 'MRC1', 'DTHD1')
sub<- data[hgnc_symbol %in% M13_sub]

mat<-as.matrix(sub[,2:(ncol(data)-1)])
rownames(mat)<-sub$hgnc_symbol
mat<-scale(mat)

lesion_color  = c('NAWM' = '#00798c', 'AL'='#d1495b', 'CA'='#edae49', 'WM'='#66a182', 'RL'='#2e4057', 'IL' = '#c02ad1')

annotation$ms<-ifelse(annotation$diagnosis=='diagnosis: Control', 'C', 'MS')
column_ha = HeatmapAnnotation(Diagnosis = annotation$ms, 
                              `Lesion type` = annotation$lesion_type, 
                              col = list(Diagnosis = c('MS'='#3D4849', 'C'='#D7D9D9'), `Lesion type` = lesion_color),
                              annotation_legend_param = list(
                                Diagnosis = list(direction = "horizontal", nrow=1),
                                `Lesion type` = list( direction='horizontal', nrow=2),
                                labels_gp = gpar(fontsize=12)),
                              annotation_name_gp = gpar( fontsize = 18),
                              simple_anno_size = unit(1, "cm")
)

k= 6
hm<-Heatmap(mat, 
            show_column_names  = T, 
            column_km  = k, 
            column_title = c('', '', '','', '', ''),
            top_annotation = column_ha,
            show_row_names = T,
            show_column_dend = F,
            show_row_dend = F,
            name = 'Expression Level',
            row_names_gp = gpar(fontsize = 18),
            raster_by_magick = TRUE,
            raster_quality  = 5,
            column_labels = annotation$patient,
            heatmap_legend_param = list(legend_height = unit(2, "cm"), 
                                        legend_width = unit(4, "cm"),
                                        direction='vertical',
                                        labels_gp = gpar(fontsize = 10)
            ))
pdf(file = file.path(out.dir, paste0('M13_',paste0(M13_sub, collapse = '_'),'.pdf')), width = 18, height = 2.5)


ht = draw(hm, merge_legend = TRUE)
dev.off()




#########################################################################################################################
M8_sub<-c('MAFB', 'KCTD8', 'MEF2C-AS1', 'RGS16')
sub<- data[hgnc_symbol %in% M8_sub]


mat<-as.matrix(sub[,2:(ncol(data)-1)])
rownames(mat)<-sub$hgnc_symbol

mat<-scale(mat)

lesion_color  = c('NAWM' = '#00798c', 'AL'='#d1495b', 'CA'='#edae49', 'WM'='#66a182', 'RL'='#2e4057', 'IL' = '#c02ad1')

annotation$ms<-ifelse(annotation$diagnosis=='diagnosis: Control', 'C', 'MS')
column_ha = HeatmapAnnotation(Diagnosis = annotation$ms, 
                              `Lesion type` = annotation$lesion_type, 
                              col = list(Diagnosis = c('MS'='#3D4849', 'C'='#D7D9D9'), `Lesion type` = lesion_color),
                              annotation_legend_param = list(
                                Diagnosis = list(direction = "horizontal", nrow=1),
                                `Lesion type` = list( direction='horizontal', nrow=2),
                                labels_gp = gpar(fontsize=12)),
                              annotation_name_gp = gpar( fontsize = 18),
                              simple_anno_size = unit(1, "cm")
)

k= 4
hm<-Heatmap(mat, 
            show_column_names  = T, 
            column_km  = k, 
            #column_title = c('', '', '','', ''),
            top_annotation = column_ha,
            show_row_names = T,
            show_column_dend = F,
            show_row_dend = F,
            name = 'Expression Level',
            row_names_gp = gpar(fontsize = 18),
            raster_by_magick = TRUE,
            raster_quality  = 5,
            column_labels = annotation$patient,
            heatmap_legend_param = list(legend_height = unit(2, "cm"), 
                                        legend_width = unit(4, "cm"),
                                        direction='vertical',
                                        labels_gp = gpar(fontsize = 10)
            ))
pdf(file = file.path(out.dir, paste0('M8_',paste0(M8_sub, collapse = '_'),'.pdf')), width = 18, height = 2.5)

ht = draw(hm, merge_legend = TRUE)
dev.off()

