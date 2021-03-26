IsThisWebServer = FALSE
################################################################################
# Copyright (C) 2009+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the follwing link:             #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################
# This is the interfacing script to the CH3Shift program, that predicts        #
# 1H and 13C chemical shifts of protein methyl groups.                         #
################################################################################

# Flags:
doval   = TRUE
x11lib  = TRUE           

# # # FLAGS # # # FLAGS # # # FLAGS # # # FLAGS # # # FLAGS
if(IsThisWebServer==TRUE) {

  if(is.na(xlimmin) | is.na(xlimmax) | is.na(ylimmin) | is.na(ylimmax)){
    xlim <- ylim <- NULL
  } else {
    xlim <- c(xlimmax, xlimmin)
    ylim <- c(ylimmax, ylimmin)
  }
  noconstraints=NULL  # if NULL, default will be used!
  minRerefNum  =NULL
  outputname   ="out.txt"
  zooming      =FALSE   
  interactive  =FALSE
  viewplot     =FALSE
  saveplot     =TRUE                   
  plotname     ="hsqc.jpg"         
  Cmargin      =NULL      
  Hmargin      =NULL
  default.csCrange=NULL
  default.csHrange=NULL
  cex1=NULL;  lwd1=NULL
  cex2=NULL;  lwd2=NULL
  cex3=NULL;  lwd3=NULL
  # # #
  prefix.ws = "../../../"
  source("../../../CH3Shift_lib/Defaults.R")
  load("../../../CH3Shift_lib/CH3Shift.ini")

} else {

  source("command.cmd")
  load("CH3Shift_lib/CH3Shift.ini")
  # generates and reads the commands based on created cmd line
    eval(parse(text=paste("CH3Shift.cmd(",cmd,")")))
    load("CH3Shift.cmd")
    file.remove("CH3Shift.cmd")
  # 
  source("CH3Shift_lib/Defaults.R")
  prefix.ws <- NULL   #"../../../"  for web server version
  
}
# # # END of FLAGS # # # END of FLAGS # # # END of FLAGS # # #



# START
OUT <- "NOTE: * * * CH3Shift * * * CH3Shift * * *  "
OUT[2] <- paste("NOTE: Title -- ", title, sep="")
write(OUT, file="process_info.txt")



# # # # # # # # ####### PDB file examination: ####### # # # # # # # # # #
  isNMR <- isNMRpdb(pdbname)
  OUT <- c(OUT, paste("NOTE: Hydrogen atom presence - ",isNMR$isH,sep=""))
         print(OUT[length(OUT)]) # PRINTOUT
  OUT <- c(OUT, paste("NOTE: Number of structures in the model - ",isNMR$Nmodels,sep=""))
         print(OUT[length(OUT)]) # PRINTOUT
  models2analyse <- 1:isNMR$Nmodels

  if(isNMR$isH!=TRUE){
    OUT <- c(OUT,"ERROR: Please add hydrogens to your molecule!")
           stop(OUT[length(OUT)]) # PRINTOUT
  }
  if(isNMR$Nmodels > 1 & examine==0){
    OUT <- c(OUT, paste("NOTE: All the ",isNMR$Nmodels,
             " models will be analyzed with the shifts averaged!",sep=""))
             print(OUT[length(OUT)]) # PRINTOUT
    models2analyse <- 1:isNMR$Nmodels
  }

  if(isNMR$Nmodels > 1 & examine > 0 & examine <= isNMR$Nmodels){
    OUT <- c(OUT,paste("NOTE: Only the structure ",examine," will be analyzed.",sep=""))
           print(OUT[length(OUT)]) # PRINTOUT
    isNMR$Nmodels <- 1
    models2analyse <- examine
  }
  
  write(OUT[3:length(OUT)], file="process_info.txt", append=TRUE)
# # # # # # # # ######## End of PDB file examination: ####### # # # # # # 



# Loading the parameters:
OUT <- c(OUT, "NOTE: Loading the parameters..."); print(OUT[length(OUT)])
write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
load(paste(prefix.ws,"CH3Shift_lib/CH3Shift.par",sep=""))
OUT <- c(OUT, "NOTE: Parameters are loaded."); print(OUT[length(OUT)])
write(OUT[length(OUT)], file="process_info.txt", append=TRUE)



#########################################################################
#########################################################################
count.tried.model <- 1
for(model in models2analyse){

  # Generating the file for a particular model (model.pdb)
  file.copy(pdbname, "model.pdb", overwrite=TRUE)

  #--
  if(model==1){                                                                          ###15OCT12###
    resName.ForQuickCheck <- pdbparser.simple("model.pdb")$resName                       ###15OCT12###
    wrong.naming.ind <- c(which(resName.ForQuickCheck=="LYP"),                           ###15OCT12###
                          which(resName.ForQuickCheck=="CYN"),                           ###15OCT12###
                          which(resName.ForQuickCheck=="CYS2"),                          ###15OCT12###
                          which(resName.ForQuickCheck=="HID"),                           ###15OCT12###
                          which(resName.ForQuickCheck=="HIE"),                           ###15OCT12###
                          which(resName.ForQuickCheck=="HISE"),                          ###15OCT12###
                          which(resName.ForQuickCheck=="HISD"))                          ###15OCT12###
    if(length(wrong.naming.ind)!=0){                                                     ###15OCT12###
      OUT <- c(OUT, "ERROR: The pdb file uses force field related naming for residues,") ###15OCT12###
      OUT <- c(OUT, "ERROR: such as LYP, CYN, CYS2, HID, HIE, HISE, HIDE etc. Please")   ###15OCT12###
      OUT <- c(OUT, "ERROR: change them into normal amino acid naming convention.")      ###15OCT12###
      OUT <- c(OUT, "ERROR: Also, make sure to have delta-protonated histidine in your") ###15OCT12###
      OUT <- c(OUT, "ERROR: structure (still named as HIS) for an extra precision!")     ###15OCT12###
      write(OUT[(length(OUT)-4):length(OUT)], file="process_info.txt", append=TRUE)      ###15OCT12###
      stop(OUT[(length(OUT)-4):length(OUT)])                                             ###15OCT12###
    }                                                                                    ###15OCT12###
  }                                                                                      ###15OCT12###
  #--


  isNMRpdb("model.pdb", model.num=model)

  # ###### Parsing the PDB file: ######
  if(pdbtype=="complex"){
    OUT <- c(OUT, "NOTE: Parsing with a 'complex' flag."); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    parsed.pdb <- pdbparser.complex("model.pdb")$PROTEIN
    parsed.pdb$resSeq <- parsed.pdb$resSeq + seqshift
  }
  if(pdbtype=="almost"){
    OUT <- c(OUT, "NOTE: Parsing with an 'almost' flag."); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    parsed.pdb <- pdbparser.almost("model.pdb")
    parsed.pdb$resSeq <- parsed.pdb$resSeq + seqshift
  }
  if(pdbtype=="medium"){
    OUT <- c(OUT, "NOTE: Parsing with a 'medium' flag."); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    parsed.pdb <- pdbparser.medium("model.pdb")
    parsed.pdb$resSeq <- parsed.pdb$resSeq + seqshift
  }
  if(pdbtype=="simple"){  # At present is the same as medium
    OUT <- c(OUT, "NOTE: Parsing with a 'simple' flag."); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    parsed.pdb <- pdbparser.simple("model.pdb")
    parsed.pdb$resSeq <- parsed.pdb$resSeq + seqshift
  }



  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Converting the $name convention to the one used in ALMOST, which is the
  #    first record in topology.lib:
  OUT <- c(OUT, "NOTE: Converting the PDB naming convention."); print(OUT[length(OUT)])
  write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
  # assigning "UN" for the chainIDs that are either "" or NA:                    ###18AUG11###
  repair.chain.ind <- which(parsed.pdb$chainID=="" | is.na(parsed.pdb$chainID))  ###18AUG11###
  if(length(repair.chain.ind) > 0){                                              ###18AUG11###
    parsed.pdb$chainID[repair.chain.ind] <- "UN"                                 ###18AUG11###
  }; rm(repair.chain.ind)                                                        ###18AUG11###
  #                                                                              ###18AUG11###
  record.variants <- pdbProc(parsed.pdb, topology=proc.topol)
  parsed.pdb <- convert.topology(parsed.pdb=parsed.pdb, proc.topol=proc.topol,
                                 record.variants=record.variants)

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  if(count.tried.model==1){
    # Finding out the methyl group presence:
    # Will be done for only the first model, and of course the rest are
    # assumed to be the same:
    Methyl.res <- Methyl.seq <- Methyl.chain <- NULL

    ala.pres <- amac.find(parsed.pdb, AmAc="ALA")
    if(ala.pres$num!=0){
      Methyl.res   <- c(Methyl.res,   rep(ala.pres$resName, ala.pres$num))
      Methyl.seq   <- c(Methyl.seq,   ala.pres$resSeq)
      Methyl.chain <- c(Methyl.chain, ala.pres$chainID)
      OUT <- c(OUT, paste("NOTE: ",ala.pres$info,sep="")); print(OUT[length(OUT)])
      write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    }
    val.pres <- amac.find(parsed.pdb, AmAc="VAL")
    if(val.pres$num!=0){
      Methyl.res   <- c(Methyl.res,   rep(val.pres$resName, val.pres$num))
      Methyl.seq   <- c(Methyl.seq,   val.pres$resSeq)
      Methyl.chain <- c(Methyl.chain, val.pres$chainID)
      OUT <- c(OUT, paste("NOTE: ",val.pres$info,sep="")); print(OUT[length(OUT)])
      write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    }
    leu.pres <- amac.find(parsed.pdb, AmAc="LEU")
    if(leu.pres$num!=0){
      Methyl.res   <- c(Methyl.res,   rep(leu.pres$resName, leu.pres$num))
      Methyl.seq   <- c(Methyl.seq,   leu.pres$resSeq)
      Methyl.chain <- c(Methyl.chain, leu.pres$chainID)
      OUT <- c(OUT, paste("NOTE: ",leu.pres$info,sep="")); print(OUT[length(OUT)])
      write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    }
    ile.pres <- amac.find(parsed.pdb, AmAc="ILE")
    if(ile.pres$num!=0){
      Methyl.res   <- c(Methyl.res,   rep(ile.pres$resName, ile.pres$num))
      Methyl.seq   <- c(Methyl.seq,   ile.pres$resSeq)
      Methyl.chain <- c(Methyl.chain, ile.pres$chainID)
      OUT <- c(OUT, paste("NOTE: ",ile.pres$info,sep="")); print(OUT[length(OUT)])
      write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    }
    thr.pres <- amac.find(parsed.pdb, AmAc="THR")
    if(thr.pres$num!=0){
      Methyl.res   <- c(Methyl.res,   rep(thr.pres$resName, thr.pres$num))
      Methyl.seq   <- c(Methyl.seq,   thr.pres$resSeq)
      Methyl.chain <- c(Methyl.chain, thr.pres$chainID)
      OUT <- c(OUT, paste("NOTE: ",thr.pres$info,sep="")); print(OUT[length(OUT)])
      write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    }
    met.pres <- amac.find(parsed.pdb, AmAc="MET")
    if(met.pres$num!=0){
      Methyl.res   <- c(Methyl.res,   rep(met.pres$resName, met.pres$num))
      Methyl.seq   <- c(Methyl.seq,   met.pres$resSeq)
      Methyl.chain <- c(Methyl.chain, met.pres$chainID)
      OUT <- c(OUT, paste("NOTE: ",met.pres$info,sep="")); print(OUT[length(OUT)])
      write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
    }

    maxnumofentry <- ( (2*(ala.pres$num + thr.pres$num + met.pres$num))+
                       (4*(val.pres$num + leu.pres$num + ile.pres$num)) )
    ###############
    RESID <- SEQ  <- NUCL <- CHAINID <- INTERCHAIN.firstM <- rep(NA, maxnumofentry)
    ChSh.pred <- SE_fit <- SD_train <- list(NULL)
    ChSh.pred[[maxnumofentry+1]] <- 
    SE_fit[[maxnumofentry+1]]    <- 
    SD_train[[maxnumofentry+1]]  <- NA
    ###############

  } else { #of if(model==1)
    OUT <- c(OUT, paste("NOTE: Model ",model," of ",isNMR$Nmodels,".",sep="")); print(OUT[length(OUT)])
    write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
  }
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



  # # # #  DISTANCE ####### ###### DISTANCE ###### ####### DISTANCE # # # #  
  # Generating the distance cell names, which will hold the further distance data:
  amac.list <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS",
                 "ILE", "LEU", "LYS", "MET", "PHE", "VAL", "SER", "THR",
                 "TYR", "TRP", "GLY", "PRO")                      
  obj <- NULL
  for(i in amac.list) {  
    i1 <- i
    if(i == "CYS"){i1 <- "CYN"}
    if(i == "HIS"){i1 <- "HID"}  #-- will be initialized as HIS but with atom names of HID                    ###18AUG11###
                                 #-- so that additional step is needed to create a dummy objects for "HE2",   ###18AUG11###
                                 #-- in case the users forget to convert all HIS to HID, so that HE2 remains. ###18AUG11###
    m <- topol.retrieve(AmAc=i1,proc.topol=proc.topol, nvar=1)
    resid <- m[,"residues"]
    rec <- m[,"record"]
    rec <- c(rec, "NAN")                                                                                      ###18AUG11###
    if(i == "HIS"){ rec <- c(rec,"HE2") }                                                                     ###18AUG11###
    for(k in 1:length(rec)){
      obj <- c(obj, paste(i, ".", rec[k],sep=""),
                    paste(i, ".", rec[k],".R1",sep=""),  
                    paste(i, ".", rec[k],".R3",sep=""),
                    paste(i, ".", rec[k],".R6",sep="")  )
    }
  }; distance.cell.names <- obj; rm(resid,m,rec,i,i1,k,obj)
  # # # #

  # Generating the cell names for the other data holders
  dih.angle.names     <- c("PHI","PSI","CHI1","CHI2","CHI3")
  dihedral.cell.names <- c("ROT1","ROT2","ROT3") # only rotamer data!!!
  ring.cell.names     <- c("RING.HIS","RING.TRP5","RING.TRP6","RING.TYR","RING.PHE")
  aniso.cell.names    <- c("ANISO.ASP","ANISO.ASN","ANISO.GLU","ANISO.GLN","ANISO.ARG",
                           "ANISO.PEPT")
  ef.cell.names       <- c("EF.FARCC")
  info.cell.names     <- c("RESSEQ", "AMACID", "NMRNUCL", "CHAINID")
  # # # #

  # All the cell names (besides info and dihedral angles!!!), ready for initialization.
  all.cell.names <- c(dih.angle.names, dihedral.cell.names, ring.cell.names, aniso.cell.names,
                      ef.cell.names, distance.cell.names)

  # Initializing all the cells:
  for(i in all.cell.names){
    # 2* for storing both C and H data and 4* if 2 methyl groups exist:
    eval(parse(text=paste(i, "<-rep(0,", maxnumofentry, ")",sep="")))
  }; rm(i)
  # # # #

  # # # info.cell.names is initialized and stored only once (for the first model!)
  if(count.tried.model==1){

    for(i in info.cell.names){
      # 2* for storing both C and H data and 4* if 2 methyl groups exist:
      eval(parse(text=paste(i, "<-rep(0,", maxnumofentry, ")",sep="")))
    }; rm(i)

  }
  # # # #  DISTANCE ####### ###### DISTANCE ###### ####### DISTANCE # # # #



  OUT <- c(OUT, "NOTE: Extracting the terms for prediction..."); print(OUT[length(OUT)])
  write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
  row.counter <- 1
  for(methyl.ind in 1:length(Methyl.seq)){
    dihedrals <- all.dihedral(Methyl.seq[methyl.ind], parsed.pdb, READLIB=readlib.dih, chainID=Methyl.chain[methyl.ind])
    ch3.names <- methyl.name(Methyl.res[methyl.ind], resSeqN=Methyl.seq[methyl.ind], parsed.pdb, chainID=Methyl.chain[methyl.ind]) 

    # *****Individual record (single entry)
    for(ch3ind in 1:ch3.names$ch3NUM){
      #--------------------------------
      if(count.tried.model==1){
        # Storing identification info:
        RESSEQ[row.counter]  <- RESSEQ[row.counter+1]  <- Methyl.seq[methyl.ind]   #C&H
        AMACID[row.counter]  <- AMACID[row.counter+1]  <- Methyl.res[methyl.ind]   #C&H
        CHAINID[row.counter] <- CHAINID[row.counter+1] <- Methyl.chain[methyl.ind] #C&H
        NMRNUCL[row.counter]   <- ch3.names$C[[ch3ind]]                            #C - row.counter
        NMRNUCL[row.counter+1] <- ch3.names$H[[ch3ind]][1]                         #H - row.counter+1
        #
      }
      #--------------------------------
      query.xyz  <- ch3.names$Cxyz[[ch3ind]]
      query.Hxyz <- ch3.names$Havxyz[[ch3ind]]
      areal <- getarea(xyz=query.xyz, parsed.pdb=parsed.pdb, r=6.5, charge=TRUE, proc.topol=proc.topol)
      INTERCHAIN.firstM[row.counter] <- INTERCHAIN.firstM[row.counter+1] <- areal$interchain
      DistDepos <- CH3.Merge.Same.Type.dist(arealt=areal, cutoffR=1.8)
      # -*-* Storing distance information:
      for(i in 1:length(DistDepos$name)) {
        if(DistDepos$name[i]!="OT1" &  DistDepos$name[i]!="OT2" &
           DistDepos$name[i]!="HT1" &  DistDepos$name[i]!="HT2" &
           DistDepos$name[i]!="HT3" &  DistDepos$name[i]!="OXT" &
           DistDepos$name[i]!="H1"  &  DistDepos$name[i]!="H2"  &
           DistDepos$name[i]!="H3"  &  DistDepos$name[i]!="HN1" &
           DistDepos$name[i]!="HN2" &  DistDepos$name[i]!="HN3"   ) { 
 
          if(DistDepos$name[i]!="HN2"  &  DistDepos$resName[i]!="PRO") { 
            eval(parse(text=paste(DistDepos$resName[i],".",DistDepos$name[i],"[",row.counter,"]","<-",
                                  DistDepos$resName[i],".",DistDepos$name[i],"[",row.counter+1,"]","<-",
                                  DistDepos$Sdist[i],sep=""))) 
            eval(parse(text=paste(DistDepos$resName[i],".",DistDepos$name[i],".R1[",row.counter,"]","<-",
                                  DistDepos$resName[i],".",DistDepos$name[i],".R1[",row.counter+1,"]","<-",
                                  DistDepos$SdistR1[i],sep="")))
            eval(parse(text=paste(DistDepos$resName[i],".",DistDepos$name[i],".R3[",row.counter,"]","<-",
                                  DistDepos$resName[i],".",DistDepos$name[i],".R3[",row.counter+1,"]","<-",
                                  DistDepos$SdistR3[i],sep=""))) 
            eval(parse(text=paste(DistDepos$resName[i],".",DistDepos$name[i],".R6[",row.counter,"]","<-",
                                  DistDepos$resName[i],".",DistDepos$name[i],".R6[",row.counter+1,"]","<-",
                                  DistDepos$SdistR6[i],sep="")))
          }
        }
      }; rm(i) # the distances are stored
      # -*-* # # # # # # # # # # # # # #

      # Storing dihedral angle data in a list for further smart averaging (median)
      PHI[row.counter]  <- PHI[row.counter+1]   <- dihedrals[1]
      PSI[row.counter]  <- PSI[row.counter+1]   <- dihedrals[2]
      CHI1[row.counter] <- CHI1[row.counter+1]  <- dihedrals[4]
      CHI2[row.counter] <- CHI2[row.counter+1]  <- dihedrals[5]
      CHI3[row.counter] <- CHI3[row.counter+1]  <- dihedrals[6]

      # Storing electric field data
      EF <- CH3.getEF(ch3ind, ch3.names, areal, parsed.pdb, discard.rad=1.8) 
      EF.FARCC[row.counter] <- EF.FARCC[row.counter+1] <- EF$FARCC
      rm(EF)
    
      ### For C (row.counter)
      # Storing COO and CON anisotropy data
      asp <- Aniso.Plane(query.xyz, areal, parsed.pdb, 
                         atom1="CG", atom2="OD1", atom3="OD2", Amino="ASP")$ANISO.SUM
      if(length(asp)!=0){ ANISO.ASP[row.counter] <- asp }
      #
      asn <- Aniso.Plane(query.xyz, areal, parsed.pdb, 
                         atom1="CG", atom2="OD1", atom3="ND2", Amino="ASN")$ANISO.SUM
      if(length(asn)!=0){ ANISO.ASN[row.counter] <- asn }
      #
      glu <- Aniso.Plane(query.xyz, areal, parsed.pdb, 
                         atom1="CD", atom2="OE1", atom3="OE2", Amino="GLU")$ANISO.SUM
      if(length(glu)!=0){ ANISO.GLU[row.counter] <- glu }
      #
      gln <- Aniso.Plane(query.xyz, areal, parsed.pdb, 
                         atom1="CD", atom2="OE1", atom3="NE2", Amino="GLN")$ANISO.SUM
      if(length(gln)!=0){ ANISO.GLN[row.counter] <- gln }
      #
      arg <- Aniso.Arg(query.xyz, areal, parsed.pdb, 
                         atom1="CZ", atom2="NH1", atom3="NH2", Amino="ARG")$ANISO.SUM
      if(length(arg)!=0){ ANISO.ARG[row.counter] <- arg }
      rm(asp,asn,glu,gln,arg)
  
      # Storing peptide moiety anisotropy data
      An.Ppt <- Aniso.Pept(query.xyz, Methyl.seq[methyl.ind], areal, parsed.pdb, chainID=Methyl.chain[methyl.ind]) 
      if(length(An.Ppt$ANISO.PEPT)!=0) { 
        ANISO.PEPT[row.counter] <- An.Ppt$ANISO.PEPT 
      }
    
      # Storing ring current data
      his.ring <- CH3.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="HIS")$RC.SUM 
       if(length(his.ring)!=0){ RING.HIS[row.counter] <- his.ring }
      trp5.ring <- CH3.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="TRP")$RC1.SUM 
       if(length(trp5.ring)!=0){ RING.TRP5[row.counter] <- trp5.ring }
      trp6.ring <- CH3.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="TRP")$RC2.SUM 
       if(length(trp6.ring)!=0){ RING.TRP6[row.counter] <- trp6.ring }
      tyr.ring <- CH3.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="TYR")$RC.SUM 
       if(length(tyr.ring)!=0){ RING.TYR[row.counter] <- tyr.ring }
      phe.ring <- CH3.Ring.AmAc(query.xyz, areal, parsed.pdb, rng="PHE")$RC.SUM
       if(length(phe.ring)!=0){ RING.PHE[row.counter] <- phe.ring }
      rm(his.ring, trp5.ring, trp6.ring, tyr.ring, phe.ring)

      ### For H (row.counter+1 AND query.Hxyz)
      # Storing COO and CON anisotropy data
      asp <- Aniso.Plane(query.Hxyz, areal, parsed.pdb, 
                         atom1="CG", atom2="OD1", atom3="OD2", Amino="ASP")$ANISO.SUM
      if(length(asp)!=0){ ANISO.ASP[row.counter+1] <- asp }
      #
      asn <- Aniso.Plane(query.Hxyz, areal, parsed.pdb, 
                         atom1="CG", atom2="OD1", atom3="ND2", Amino="ASN")$ANISO.SUM
      if(length(asn)!=0){ ANISO.ASN[row.counter+1] <- asn }
      #
      glu <- Aniso.Plane(query.Hxyz, areal, parsed.pdb, 
                         atom1="CD", atom2="OE1", atom3="OE2", Amino="GLU")$ANISO.SUM
      if(length(glu)!=0){ ANISO.GLU[row.counter+1] <- glu }
      #
      gln <- Aniso.Plane(query.Hxyz, areal, parsed.pdb, 
                         atom1="CD", atom2="OE1", atom3="NE2", Amino="GLN")$ANISO.SUM
      if(length(gln)!=0){ ANISO.GLN[row.counter+1] <- gln }
      #
      arg <- Aniso.Arg(query.Hxyz, areal, parsed.pdb, 
                         atom1="CZ", atom2="NH1", atom3="NH2", Amino="ARG")$ANISO.SUM
      if(length(arg)!=0){ ANISO.ARG[row.counter+1] <- arg }
      rm(asp,asn,glu,gln,arg)
 
      # Storing peptide moiety anisotropy data
      An.Ppt <- Aniso.Pept(query.Hxyz, Methyl.seq[methyl.ind], areal, parsed.pdb, chainID=Methyl.chain[methyl.ind]) 
      if(length(An.Ppt$ANISO.PEPT)!=0) { 
        ANISO.PEPT[row.counter+1] <- An.Ppt$ANISO.PEPT
      }
      
      # Storing ring current data
      his.ring <- CH3.Ring.AmAc(query.Hxyz, areal, parsed.pdb, rng="HIS")$RC.SUM 
       if(length(his.ring)!=0){ RING.HIS[row.counter+1] <- his.ring }
      trp5.ring <- CH3.Ring.AmAc(query.Hxyz, areal, parsed.pdb, rng="TRP")$RC1.SUM 
       if(length(trp5.ring)!=0){ RING.TRP5[row.counter+1] <- trp5.ring }
      trp6.ring <- CH3.Ring.AmAc(query.Hxyz, areal, parsed.pdb, rng="TRP")$RC2.SUM 
       if(length(trp6.ring)!=0){ RING.TRP6[row.counter+1] <- trp6.ring }
      tyr.ring <- CH3.Ring.AmAc(query.Hxyz, areal, parsed.pdb, rng="TYR")$RC.SUM 
       if(length(tyr.ring)!=0){ RING.TYR[row.counter+1] <- tyr.ring }
      phe.ring <- CH3.Ring.AmAc(query.Hxyz, areal, parsed.pdb, rng="PHE")$RC.SUM
       if(length(phe.ring)!=0){ RING.PHE[row.counter+1] <- phe.ring }
      rm(his.ring, trp5.ring, trp6.ring, tyr.ring, phe.ring)
      # Shifting the counter
      row.counter <- row.counter + 2
      
      if(count.tried.model==1){
        # storing only for the first model:
        RESID <- AMACID
        SEQ   <- RESSEQ
        NUCL  <- NMRNUCL
        #CHAINID <- CHAINID
        #INTERCHAIN.firstM
      }

    }
    # *****Individual record (single entry)

  } # all the methyls from a single model are recorded



  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # PREDICTING AND STORING CHEMICAL SHIFT DATA



  # Inferring the rotamer types:
  for(i in 1:length(CHI1)){
    chi1.angle <- CHI1[i]
    # # # #
    if(!is.na(chi1.angle)) {  # case when data are not NA (non-Ala)
      if(chi1.angle>-120 & chi1.angle<=0){ROT1[i]<-1}
      if(chi1.angle>0 & chi1.angle<=120){ROT2[i]<-1}
      if(ROT1[i]==0 & ROT2[i]==0){ROT3[i]<-1}
    }
    # # # #
  }; rm(i, chi1.angle)
  # # # # # # # # # # #
  
  # necessary especially for some PHI/PSI/CHIs that return NA
  PHI[is.na(PHI)] <- 0.00
  PSI[is.na(PSI)] <- 0.00
  CHI1[is.na(CHI1)] <- 0.00
  CHI2[is.na(CHI2)] <- 0.00
  CHI3[is.na(CHI3)] <- 0.00
  RING.HIS[is.na(RING.HIS)] <- 0.00
  RING.TRP5[is.na(RING.TRP5)] <- 0.00
  RING.TRP6[is.na(RING.TRP6)] <- 0.00
  RING.TYR[is.na(RING.TYR)] <- 0.00
  RING.PHE[is.na(RING.PHE)] <- 0.00
  ANISO.ASP[is.na(ANISO.ASP)] <- 0.00
  ANISO.ASN[is.na(ANISO.ASN)] <- 0.00
  ANISO.GLU[is.na(ANISO.GLU)] <- 0.00
  ANISO.GLN[is.na(ANISO.GLN)] <- 0.00
  ANISO.ARG[is.na(ANISO.ARG)] <- 0.00
  ANISO.PEPT[is.na(ANISO.PEPT)] <- 0.00
  EF.FARCC[is.na(EF.FARCC)] <- 0.00
  
  # Generating the data frame:
  DFrame <- data.frame(     RESSEQ         =as.numeric(RESSEQ),
                            AMACID         =as.character(AMACID),
                            NMRNUCL        =as.character(NMRNUCL),
                            EF.FARCC       =as.numeric(EF.FARCC),
                            ROT1           =as.numeric(ROT1),
                            ROT2           =as.numeric(ROT2),
                            ROT3           =as.numeric(ROT3),
                            ANISO.ARG      =as.numeric(ANISO.ARG),
                            ANISO.PEPT     =as.numeric(ANISO.PEPT), 
                            ANISO.AMD      =as.numeric(ANISO.ASN)+as.numeric(ANISO.GLN),
                            ANISO.ACD      =as.numeric(ANISO.ASP)+as.numeric(ANISO.GLU),
                            RING.HIS       =as.numeric(RING.HIS),
                            RING.TRP5      =as.numeric(RING.TRP5),
                            RING.TRP6      =as.numeric(RING.TRP6),
                            RING.TYR       =as.numeric(RING.TYR),
                            RING.PHE       =as.numeric(RING.PHE),
                            PHI            =(as.numeric(PHI)*pi/180),
                            PHI2           =(as.numeric(PHI)*pi/180)^2,
                            PHI3           =(as.numeric(PHI)*pi/180)^3,
                            PHI4           =(as.numeric(PHI)*pi/180)^4,
                            PSI            =(as.numeric(PSI)*pi/180),
                            PSI2           =(as.numeric(PSI)*pi/180)^2,
                            PSI3           =(as.numeric(PSI)*pi/180)^3,
                            PSI4           =(as.numeric(PSI)*pi/180)^4,
                            CHI1           =as.numeric(CHI1)*pi/180, 
                            CHI12          =(as.numeric(CHI1)*pi/180)^2,
                            CHI13          =(as.numeric(CHI1)*pi/180)^3, 
                            CHI14          =(as.numeric(CHI1)*pi/180)^4, 
                            CHI2           =as.numeric(CHI2)*pi/180,
                            CHI22          =(as.numeric(CHI2)*pi/180)^2,
                            CHI23          =(as.numeric(CHI2)*pi/180)^3,
                            CHI24          =(as.numeric(CHI2)*pi/180)^4,  
                            CHI3           =as.numeric(CHI3)*pi/180,
                            CHI32          =(as.numeric(CHI3)*pi/180)^2,
                            CHI33          =(as.numeric(CHI3)*pi/180)^3,
                            CHI34          =(as.numeric(CHI3)*pi/180)^4   )

  # Adding all the distance information (in Angstrom) to the created data frame:
  cmd <- "DFrame <- data.frame(DFrame"
  for(i in distance.cell.names) {
    cmd <- paste(cmd,", ",i,"=as.numeric(",i,")",sep="") 
  }; rm(i)
  cmd <- paste(cmd, ")",sep="")
  eval(parse(text=cmd)); rm(cmd)
  # # # #

  # Removing all the template cells except CHAINID:
  rm(list= all.cell.names[-which(all.cell.names=="CHAINID")])
  
  # Adding the missing factors:
  DFrame <- addCos(Data=DFrame)
  DFrame <- addMergeDist(DFrame)
  # Data frame is ready to be used for predictions.
  
  
  OUT <- c(OUT, "NOTE: Predicting chemical shifts..."); print(OUT[length(OUT)])
  write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
  
  
  Ncl.type <- rep(NA, length(DFrame[,1]))
  identifier <- unique(paste(as.character(DFrame$AMACID),"_",
                             as.character(DFrame$NMRNUCL),sep=""))
  for(i in identifier){
    resid <- unlist(strsplit(i,"_"))[1]   #residue of prediction
    ncl   <- unlist(strsplit(i,"_"))[2]   #nucleus of prediction
    rowind <- which(DFrame$AMACID==resid & DFrame$NMRNUCL==ncl) # rows in the DFrame
    Ncl.type[rowind] <- unlist(strsplit(ncl,""))[1]
    options(warn=-1) 
    eval(parse(text=paste("pr.lm <- predict.lmold(object=",i,"_all.methcs, newdata=DFrame[rowind,] , type='response',se.fit=TRUE)",sep="")))
    pr.lm <- CH3.filter.violation(pr.lm, ncl, resid)
    #Discarding results for MET
    if(resid=="MET"){
      pr.lm <- NULL
      pr.lm$fit <- rep(NA, length(rowind))
      pr.lm$se.fit <- rep(NA, length(rowind))
      pr.lm$residual.scale <- NA
    }
    #
    for(ii in 1:length(rowind)){
      ChSh.pred[[ rowind[ii] ]] <- c(ChSh.pred[[ rowind[ii] ]], as.vector(pr.lm$fit[ii]))
      SE_fit[[ rowind[ii] ]]    <- c(SE_fit[[ rowind[ii] ]],    as.vector(pr.lm$se.fit[ii]))
      SD_train[[ rowind[ii] ]]  <- c(SD_train[[ rowind[ii] ]],  as.vector(pr.lm$residual.scale))
    }; rm(ii)
  }; rm(i)
  # # # END OF PREDICTING AND STORING CHEMICAL SHIFT DATA
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

  count.tried.model <- count.tried.model + 1
} # all the models are recorded, but without averaging.
#########################################################################
#########################################################################
file.remove("model.pdb")



# removing last NA
ChSh.pred <- ChSh.pred[1:(length(ChSh.pred)-1)]
SE_fit    <-    SE_fit[1:(length(SE_fit)-1)]
SD_train  <-  SD_train[1:(length(SD_train)-1)]

# Averaging the data obtained for multiple conformers
Pred.CS               <- sapply(ChSh.pred, FUN=function(i){mean(i,na.rm=TRUE)}, simplify=TRUE)
Pred.sefit            <- sapply(SE_fit,    FUN=function(i){mean(i,na.rm=TRUE)}, simplify=TRUE)
Pred.SDresidual.train <- sapply(SD_train,  FUN=function(i){mean(i,na.rm=TRUE)}, simplify=TRUE)



#########################################################################
#  Filtering out the erroneous and violating predictions + reformating data
violating.shift.ind <- which(is.na(Pred.CS))
if(length(violating.shift.ind)>0){
  Pred.CS[violating.shift.ind] <- "xxxxx"
  Pred.CS[-violating.shift.ind] <- format(round(as.numeric(Pred.CS[-violating.shift.ind]),3),nsmall=3)
} else {
  Pred.CS <- format(round(Pred.CS,3),nsmall=3)
}

violating.shift.ind <- which(Pred.sefit>999.99)
if(length(violating.shift.ind)>0){
  Pred.sefit[violating.shift.ind] <- "xxxxx"
  Pred.sefit[-violating.shift.ind] <- format(round(as.numeric(Pred.sefit[-violating.shift.ind]),3),nsmall=3)
} else {
  Pred.sefit <- format(round(Pred.sefit,3),nsmall=3)
}

violating.shift.ind <- which(Pred.SDresidual.train>999.99)
if(length(violating.shift.ind)>0){
  Pred.SDresidual.train[violating.shift.ind] <- "xxxxx"
  Pred.SDresidual.train[-violating.shift.ind] <- format(round(
                                      as.numeric(Pred.SDresidual.train[-violating.shift.ind]),3),nsmall=3)
} else {
  Pred.SDresidual.train <- format(round(as.numeric(Pred.SDresidual.train),3),nsmall=3)
}
# #######################################################################

  # corrected NUCL ids for output and er.data (ie HD21 = HD2, HB1 = HB, etc.) 
  NUCL.cr <- sapply(NUCL, FUN=function(smth){
                                 smth <- unlist(strsplit(smth,""));
                                 if(length(smth)>2){
                                   smth <- paste(smth[1:3], collapse="");
                                 } else {
                                   smth <- paste(smth, collapse="");
                                 }
                                 if(smth=="HB1"){smth <- "HB"};
                                 if(smth=="HE1"){smth <- "HE"};
                                 return(as.vector(smth));
                              }, simplify=TRUE)

#########################################################################



# # # # Reading in the experimental data
if(length(experdata)!=0){
  OUT <- c(OUT, "NOTE: Processing the experimental data."); print(OUT[length(OUT)])
  write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
  exp.data <- readLines(experdata)

  # dummy and blank line exclusion is also added.
  comment.rows <- unique(c(grep("#", exp.data, fixed=TRUE),  
                           which(exp.data==""),
                           which(exp.data==" "),
                           which(exp.data=="  "),
                           which(exp.data=="   "),
                           which(exp.data=="    "),
                           which(exp.data=="     "),
                           which(exp.data=="      ")  ))

  if(length(comment.rows!=0)){ 
    exp.data <- exp.data[-comment.rows] 
  }

  exp.amac <- exp.seq <- exp.nucl <- exp.shift <- exp.chain <- NULL

  # # # # 
  for(i in 1:length(exp.data)){
    line <- linesplit(exp.data[i])
    if(length(line)==4 | length(line)==5){
      amacic    <- as.character(line[1])
      exp.amac  <- c(exp.amac, amacic)
      exp.seq   <- c(exp.seq, as.numeric(line[2]))
      nuclic    <- as.character(line[3])
      # switching to internal convention:
      if(amacic=="ALA"& nuclic=="HB"){nuclic<-"HB1"}     
      if(amacic=="VAL"& nuclic=="HG1"){nuclic<-"HG11"}
      if(amacic=="VAL"& nuclic=="HG2"){nuclic<-"HG21"}
      if(amacic=="LEU"& nuclic=="HD1"){nuclic<-"HD11"}
      if(amacic=="LEU"& nuclic=="HD2"){nuclic<-"HD21"}
      if(amacic=="ILE"& nuclic=="CD1"){nuclic<-"CD"}
      if(amacic=="ILE"& nuclic=="HD11"){nuclic<-"HD1"}
      if(amacic=="ILE"& nuclic=="HG2"){nuclic<-"HG21"}
      if(amacic=="THR"& nuclic=="HG2"){nuclic<-"HG21"}
      if(amacic=="MET"& nuclic=="HE"){nuclic<-"HE1"}
      exp.nucl  <- c(exp.nucl, nuclic)
      exp.shift <- c(exp.shift, as.numeric(line[4]))
      chainic   <- as.character(line[5])
      if(is.na(chainic)){chainic <- "UN"}
      exp.chain <- c(exp.chain, chainic)
    }
  };rm(i)
  # # # # 

  Exp.CS <- Exp_Pred.CS <- rep(NA, length(NUCL))

  for(i in 1:length(exp.shift)) {
    position <-  which( RESID   == exp.amac[i]  &
                        NUCL    == exp.nucl[i]  & 
                        SEQ     == exp.seq[i]   &
                        CHAINID == exp.chain[i]   )
    if(length(position)!=0){
      Exp.CS[position] <- round(exp.shift[i], 3)
      prediction <- as.numeric(Pred.CS[position])
      if(!is.na(prediction)){
        Exp_Pred.CS[position] <- exp.shift[i]-prediction
      }
    } else {
      OUT <- c(OUT, "ERROR: Some of the exp. data records do not match the structure,") ###15OCT12###
      OUT <- c(OUT, "ERROR: or the experimental data file is incorrectly formatted.")   ###15OCT12###
      OUT <- c(OUT, paste("ERROR: i=",i,", RESID=",exp.amac[i],", NUCL=",exp.nucl[i],   ###15OCT12###
                          ", SEQ=",exp.seq[i],", CHAINID=",exp.chain[i],sep=""))        ###15OCT12###
      write(OUT[c(length(OUT)-2,length(OUT)-1,length(OUT))],                            ###15OCT12###
           file="process_info.txt", append=TRUE)                                        ###15OCT12###
      stop(OUT[c(length(OUT)-2,length(OUT)-1,length(OUT))])                             ###15OCT12###
    } 
  }

  Exp.CS.nonfilt <- Exp.CS

  ndx <- which(is.na(Exp.CS))
  if(length(ndx)!=0){
   Exp.CS[ndx]  <- "xxxxx"
   Exp.CS[-ndx] <- format(round(as.numeric(Exp.CS[-ndx]),3),nsmall=3)
  } else {
   Exp.CS <- format(round(as.numeric(Exp.CS),3),nsmall=3)  
  }
  ndx <- which(is.na(Exp_Pred.CS))
  if(length(ndx)!=0){
   Exp_Pred.CS[ndx]  <- "xxxxxx"
   Exp_Pred.CS[-ndx] <- format(round(as.numeric(Exp_Pred.CS[-ndx]),3),nsmall=3)
  } else {
   Exp_Pred.CS <- format(round(as.numeric(Exp_Pred.CS),3),nsmall=3)  
  }
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  # # # REREFERENCING
  if(rereference==TRUE){
    Pred.CS.num <- as.numeric(Pred.CS)
    # rereferencing proton chemical shifts:
    h.reref.ind <- which(!is.na(Pred.CS.num) & !is.na(Exp.CS.nonfilt))
    exp.reref  <- Exp.CS.nonfilt[h.reref.ind]
    pred.reref <- Pred.CS.num[h.reref.ind]
    if(length(exp.reref)>minRerefNum) {
      lm.reref <- lm(formula = exp.reref ~ 1, offset = pred.reref)
      Pred.CS.num <- Pred.CS.num + lm.reref$coefficients[1]
    }
    # recalculating errors
    Exp_Pred.CS <- Exp.CS.nonfilt-Pred.CS.num

    # reformating back
    na.pred.CS <- which(is.na(Pred.CS.num))
    if(length(na.pred.CS)>0){
      Pred.CS.num[na.pred.CS] <- "xxxxx"
      Pred.CS.num[-na.pred.CS] <- format(round(as.numeric(Pred.CS.num[-na.pred.CS]),3),nsmall=3)
    } else {
      Pred.CS.num <- format(round(as.numeric(Pred.CS.num),3),nsmall=3)
    }
    Pred.CS <- Pred.CS.num
    na.err.CS <- which(is.na(Exp_Pred.CS))
    if(length(na.err.CS)>0){
     Exp_Pred.CS[na.err.CS]  <- "xxxxxx" 
     Exp_Pred.CS[-na.err.CS] <- format(round(as.numeric(Exp_Pred.CS[-na.err.CS]),3),nsmall=3)
    } else {
     Exp_Pred.CS <- format(round(as.numeric(Exp_Pred.CS),3),nsmall=3)  
    }
    #
  }
  # # # END OF REREFERENCING

  # Generating the Probs and Stars data objects that comment on probabilities
  Exp_Pred.CS.num <- as.vector(sapply(Exp_Pred.CS, fun5, simplify=TRUE))
  Probs <- rep("xxxxx", length(Exp_Pred.CS.num))
  Stars <- rep("xxx", length(Exp_Pred.CS.num))
  UNQ.SEQ.CHAINID <- paste(SEQ, CHAINID, sep="_")
  for(SEQ.CH.i in unique(UNQ.SEQ.CHAINID)) {                          
    is <- which(SEQ==as.numeric(unlist(strsplit(SEQ.CH.i,"_"))[1]) &  
                CHAINID==unlist(strsplit(SEQ.CH.i,"_"))[2] &          
                Ncl.type=="H" &                                       
                RESID!="MET"  &                                        
                !is.na(Exp_Pred.CS.num))                              
    if(length(is)!=0) {   # is can be a num or vector of nums
      # reading only the first 3 letters of the nucleus name (relevant for getP.ch3)!
      P.data <- getP.ch3(AMC=ResidueName(unique(RESID[is]), capitalize=FALSE,
                                         threeletter=TRUE, oneletter = FALSE),
                         NCL=NUCL.cr[is], De_c=Exp_Pred.CS.num[is], er.data=er.data)
      Probs[is] <- format(round(as.numeric(P.data$Prob),3),nsmall=3)  
      Stars[is] <- paste(rep("*", P.data$Stars), collapse="")
      rm(P.data)
    }
    rm(is)   
  }; rm(SEQ.CH.i)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  if(doval==TRUE){

     # Generating chimera_cmd.txt file with UCSF Chimera command to color residues
     # according to their probabilities:
     colfun <- colorRampPalette(c("red","orange","yellow","green",5,"blue"), bias=bias)
     colors <- colfun(500)
     prob.range <- seq(from=0, to=1, by=0.002001)
     # plotting pallete:
     if(x11lib==TRUE){
       jpeg(quality=100, height=500, width=500, filename="protein_prob_palette.jpg")
         plot(x=prob.range,y=rep(0,500), col=colors, cex=5, pch=3, ylim=c(-0.1,0.1),
              main=paste("Bias = ", bias,sep=""))
       dev.off()
     }
     #
     chimera.cmd <- NULL
     for(SEQ.CH.i in unique(UNQ.SEQ.CHAINID)){
       is <- which(SEQ==as.numeric(unlist(strsplit(SEQ.CH.i,"_"))[1]) & 
                   CHAINID==unlist(strsplit(SEQ.CH.i,"_"))[2] &
                   Probs!="xxxxx")
       if(length(is)!=0) {
         is <- is[1]
         prb <- as.numeric(Probs[is]) 
         colind <- which(abs(prob.range-prb) == min(abs(prob.range-prb)))  
         chimera.cmd <- c(chimera.cmd," color ",colors[colind],
                          " :",as.numeric(unlist(strsplit(SEQ.CH.i,"_"))[1]),
                           ".",unlist(strsplit(SEQ.CH.i,"_"))[2],";")
         rm(prb, colind)
       }
       rm(is)
     }; rm(SEQ.CH.i, colors, colfun, prob.range)
     
     write(paste(chimera.cmd, collapse=""), file="chimera_cmd.txt")
     rm(chimera.cmd)
  }
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Printing out the results:
  if(outorder=="bytype"){ # Prints results odered by type.
    RESULTS.rep <- 
        c("NOTE:                       (o:     PRINTING OUT THE RESULTS:     :o)                          ",
          "NOTE:         At present the results for all the CARBONS and MET should be discarded!!!        ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn   CS(ex)   D(e-p)     P     Np",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")," ; ",
          format(CHAINID, width=3, justify="centre")," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre"),"; ",
          format(NUCL.cr, width=4, justify="centre")," ; ",
          format(Pred.CS, width=6, justify="centre")," ; ",
          format(Pred.sefit, width=5, justify="centre")," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")," ; ",
          format(Exp.CS, width=6, justify="centre")," ; ",
          format(Exp_Pred.CS, width=6, justify="centre")," ; ",
          format(Probs, width=5, justify="centre")," ; ",
          format(Stars, width=3, justify="centre")
          ,sep=""))
  }
  if(outorder=="byseq"){ # Prints results odered by sequence number.
    byseq.ord <- order(SEQ)
    if(length(unique(CHAINID)) > 1){
      byseq.ord <- byseq.ord[order(CHAINID[byseq.ord])]
    }
    RESULTS.rep <- 
        c("NOTE:                       (o:     PRINTING OUT THE RESULTS:     :o)                          ",
          "NOTE:         At present the results for all the CARBONS and MET should be discarded!!!        ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn   CS(ex)   D(e-p)     P     Np",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID[byseq.ord]," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")[byseq.ord]," ; ",
          format(CHAINID, width=3, justify="centre")[byseq.ord]," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre")[byseq.ord],"; ",
          format(NUCL.cr, width=4, justify="centre")[byseq.ord]," ; ",
          format(Pred.CS, width=6, justify="centre")[byseq.ord]," ; ",
          format(Pred.sefit, width=5, justify="centre")[byseq.ord]," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")[byseq.ord]," ; ",
          format(Exp.CS, width=6, justify="centre")[byseq.ord]," ; ",
          format(Exp_Pred.CS, width=6, justify="centre")[byseq.ord]," ; ",
          format(Probs, width=5, justify="centre")[byseq.ord]," ; ",
          format(Stars, width=3, justify="centre")[byseq.ord]
          ,sep=""))
  }
  if(outorder=="byshift"){ # Prints results odered by sequence number.
    byshift.ord <- order(as.numeric(Pred.CS))
    RESULTS.rep <- 
        c("NOTE:                       (o:     PRINTING OUT THE RESULTS:     :o)                          ",
          "NOTE:         At present the results for all the CARBONS and MET should be discarded!!!        ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn   CS(ex)   D(e-p)     P     Np",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")[byshift.ord]," ; ",
          format(CHAINID, width=3, justify="centre")[byshift.ord]," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre")[byshift.ord],"; ",
          format(NUCL.cr, width=4, justify="centre")[byshift.ord]," ; ",
          format(Pred.CS, width=6, justify="centre")[byshift.ord]," ; ",
          format(Pred.sefit, width=5, justify="centre")[byshift.ord]," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")[byshift.ord]," ; ",
          format(Exp.CS, width=6, justify="centre")[byshift.ord]," ; ",
          format(Exp_Pred.CS, width=6, justify="centre")[byshift.ord]," ; ",
          format(Probs, width=5, justify="centre")[byshift.ord]," ; ",
          format(Stars, width=3, justify="centre")[byshift.ord]
          ,sep=""))
  }

  rm(Probs, Stars, Exp_Pred.CS.num, fun5, er.data, getP)

} else {  # of if(length(experdata)!=0) => experimental data does not exist!!!

  if(outorder=="bytype"){ # Prints results odered by type.
    RESULTS.rep <- 
        c("NOTE:          (o:     PRINTING OUT THE RESULTS:     :o)       ",
          "NOTE:   The results for CARBONS and MET should be discarded!!! ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")," ; ",
          format(CHAINID, width=3, justify="centre")," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre"),"; ",
          format(NUCL.cr, width=4, justify="centre")," ; ",
          format(Pred.CS, width=6, justify="centre")," ; ",
          format(Pred.sefit, width=5, justify="centre")," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")
          ,sep=""))
  }
  if(outorder=="byseq"){ # Prints results odered by thype.
    byseq.ord <- order(SEQ)
    if(length(unique(CHAINID)) > 1){
      byseq.ord <- byseq.ord[order(CHAINID[byseq.ord])]
    }
    RESULTS.rep <- 
        c("NOTE:          (o:     PRINTING OUT THE RESULTS:     :o)       ",
          "NOTE:   The results for CARBONS and MET should be discarded!!! ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID[byseq.ord]," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")[byseq.ord]," ; ",
          format(CHAINID, width=3, justify="centre")[byseq.ord]," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre")[byseq.ord],"; ",
          format(NUCL.cr, width=4, justify="centre")[byseq.ord]," ; ",
          format(Pred.CS, width=6, justify="centre")[byseq.ord]," ; ",
          format(Pred.sefit, width=5, justify="centre")[byseq.ord]," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")[byseq.ord]
          ,sep=""))
  }
  if(outorder=="byshift"){ # Prints results odered by thype.
    byshift.ord <- order(as.numeric(Pred.CS))
    RESULTS.rep <- 
        c("NOTE:          (o:     PRINTING OUT THE RESULTS:     :o)       ",
          "NOTE:   The results for CARBONS and MET should be discarded!!! ",
    paste("NOTE:    RESID   SEQ   CHAIN IN1  NUCL   CS(pr)   SEfit   SDtrn",sep=""),
    paste(rep("RESULT:", length(NUCL))," ; ",
          RESID[byshift.ord]," ; ",
          format(SEQ, nsmall=0, width=5, justify="centre")[byshift.ord]," ; ",
          format(CHAINID, width=3, justify="centre")[byshift.ord]," ;",
          format(INTERCHAIN.firstM, width=3, justify="centre")[byshift.ord],"; ",
          format(NUCL.cr, width=4, justify="centre")[byshift.ord]," ; ",
          format(Pred.CS, width=6, justify="centre")[byshift.ord]," ; ",
          format(Pred.sefit, width=5, justify="centre")[byshift.ord]," ; ",
          format(Pred.SDresidual.train, width=5, justify="centre")[byshift.ord]
          ,sep=""))
  }
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

OUT <- c(OUT, RESULTS.rep, 
         "NOTE: CH3Shift - Good luck!",
         "NOTE: For questions and bug reports contact - aleksahak@cantab.net ",
         "NOTE: N.B. Please, attach all the related files to the e-mail.   ")
write(OUT, file=outputname)
print("NOTE: Output file is written!")
write("NOTE: Output file is written!", file="process_info.txt", append=TRUE)




# RESID; SEQ; NUCL.cr; Pred.CS; Pred.sefit; Pred.SDresidual.train; Ncl.type
#########################################################################
#####                           PLOTTING                           ######
#########################################################################


### Function borrowd from the "calibrate" R package!
textxy <- function (X, Y, labs, cx = 0.5, dcol = "black", m = c(0, 0)) {
    posXposY <- ((X >= m[1]) & ((Y >= m[2])))
    posXnegY <- ((X >= m[1]) & ((Y < m[2])))
    negXposY <- ((X < m[1]) & ((Y >= m[2])))
    negXnegY <- ((X < m[1]) & ((Y < m[2])))
    if (sum(posXposY) > 0) 
        text(X[posXposY], Y[posXposY], labs[posXposY], adj = c(-0.3, 
            -0.3), cex = cx, col = dcol)
    if (sum(posXnegY) > 0) 
        text(X[posXnegY], Y[posXnegY], labs[posXnegY], adj = c(-0.3, 
            1.3), cex = cx, col = dcol)
    if (sum(negXposY) > 0) 
        text(X[negXposY], Y[negXposY], labs[negXposY], adj = c(1.3, 
            -0.3), cex = cx, col = dcol)
    if (sum(negXnegY) > 0) 
        text(X[negXnegY], Y[negXnegY], labs[negXnegY], adj = c(1.3, 
            1.3), cex = cx, col = dcol)
}
###


print("NOTE: Plotting the results...")
 write("NOTE: Plotting the results...", file="process_info.txt", append=TRUE)
Pred.CS <- as.numeric(Pred.CS)
if(length(experdata)!=0){
  Exp.CS  <- as.numeric(Exp.CS)
}

Amacid.CS <- as.character(DFrame[,"AMACID"])
Resseq.CS <- as.numeric(DFrame[,"RESSEQ"])
Nmrnucl.CS <- as.character(DFrame[,"NMRNUCL"])
Chainid.CS <- CHAINID

# Determining plot range:
if(length(experdata)!=0){
  plotrange.csC <- c(Pred.CS[which(Ncl.type=="C")],Exp.CS[which(Ncl.type=="C")])
  plotrange.csC <- range(plotrange.csC[!is.na(plotrange.csC)])
  if(plotrange.csC[1]==Inf | plotrange.csC[2]==Inf) {
    plotrange.csC <- rev(default.csCrange)
  } else {
    plotrange.csC <- rev(c(plotrange.csC[1]-Cmargin, plotrange.csC[2]+Cmargin))
  }

  plotrange.csH <- c(Pred.CS[which(Ncl.type=="H")],Exp.CS[which(Ncl.type=="H")])
  plotrange.csH <- range(plotrange.csH[!is.na(plotrange.csH)])
  if(plotrange.csH[1]==Inf | plotrange.csH[2]==Inf) {
    plotrange.csH <- rev(default.csHrange)
  } else {
    plotrange.csH <- rev(c(plotrange.csH[1]-Hmargin, plotrange.csH[2]+Hmargin))
  }
  
} else {

  plotrange.csC <- Pred.CS[which(Ncl.type=="C")]
  plotrange.csC <- range(plotrange.csC[!is.na(plotrange.csC)])
  if(plotrange.csC[1]==Inf | plotrange.csC[2]==Inf) {
    plotrange.csC <- rev(default.csCrange)
  } else {
    plotrange.csC <- rev(c(plotrange.csC[1]-Cmargin, plotrange.csC[2]+Cmargin))
  }

  plotrange.csH <- Pred.CS[which(Ncl.type=="H")]
  plotrange.csH <- range(plotrange.csH[!is.na(plotrange.csH)])
  if(plotrange.csH[1]==Inf | plotrange.csH[2]==Inf) {
    plotrange.csH <- rev(default.csHrange)
  } else {
    plotrange.csH <- rev(c(plotrange.csH[1]-Hmargin, plotrange.csH[2]+Hmargin))
  }
}
# END of determining plot range.



# # # # # # # # # # # #
if(length(xlim)==0){
  xlim <- plotrange.csH
}
if(length(ylim)==0){
  ylim <- plotrange.csC
}
# # # # # # # # # # # #

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # function to make more user friendly label names for plotting
    newlabel <- function(label){
       newlab <- NULL
       for(lblb in 1:length(newlabel)){
         labline <- label[lblb]
         if(labline=="CE") {newlab[lblb] <- "e" }
         if(labline=="CB") {newlab[lblb] <- "b" }
         if(labline=="CG1"){newlab[lblb] <- "g1"}
         if(labline=="CG2"){newlab[lblb] <- "g2"}
         if(labline=="CD") {newlab[lblb] <- "d1"}
         if(labline=="CD1"){newlab[lblb] <- "d1"}
         if(labline=="CD2"){newlab[lblb] <- "d2"}
       }
       return(newlab)
     }
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

lineX <- lineY <- pointX <- pointY <- pointLAB <- lineLAB.X <- 
lineLAB.Y <- colorsP <- colorsLX <- colorsLY <-  NULL
for(i in 1:length(Pred.CS)) {
  if(!is.na(Pred.CS[i])) {           # Pred.CS[i] is not NA
    eval(parse(text=paste("col <- col",Amacid.CS[i],sep="")))
    ######################
    if(Ncl.type[i]=="H") {
      x <- Pred.CS[i]
      Nucltip <- Nmrnucl.CS[i]
      if(Nucltip=="HE1"){Neighb.nucltip <- "CE"}
      if(Nucltip=="HB1"){Neighb.nucltip <- "CB"}
      if(Nucltip=="HG11"){Neighb.nucltip <- "CG1"}
      if(Nucltip=="HG21"){Neighb.nucltip <- "CG2"}
      if(Nucltip=="HD1"){Neighb.nucltip <- "CD"}
      if(Nucltip=="HD11"){Neighb.nucltip <- "CD1"}
      if(Nucltip=="HD21"){Neighb.nucltip <- "CD2"}
      Cor.nucl <- "C"
      Cor.CS.ind <- which(Ncl.type==Cor.nucl & 
                          Amacid.CS==Amacid.CS[i] &
                          Resseq.CS==Resseq.CS[i] &
                          Chainid.CS==Chainid.CS[i] &
                          Nmrnucl.CS==Neighb.nucltip)
      if(length(Cor.CS.ind)==1) {
        y <- Pred.CS[Cor.CS.ind]
        if(!is.na(y)){
          lab <- paste(Amacid.CS[i],"-",newlabel(Neighb.nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
          pointX <- c(pointX,x)
          pointY <- c(pointY,y)
          pointLAB <- c(pointLAB, lab)
          colorsP <- c(colorsP, col)
          #points(x=x, y=y,col=col, cex=cex1, lwd=lwd1)
          #points(x=x, y=y,col=col, cex=cex2, lwd=lwd2)
          #points(x=x, y=y,col=col, cex=cex3, lwd=lwd3)
          #textxy(X=x, Y=y, labs=lab)
        } else {
          lab <- paste(Amacid.CS[i],"-",newlabel(Neighb.nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
          lineX <- c(lineX,x)
          lineLAB.X <- c(lineLAB.X, lab)
          colorsLX <- c(colorsLX, col)
          #abline(v=x,lwd=lwd1, lty="dashed", col=col)
        }
      } else {
        lab <- paste(Amacid.CS[i],"-",newlabel(Neighb.nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
        lineX <- c(lineX,x)
        lineLAB.X <- c(lineLAB.X, lab)
        colorsLX <- c(colorsLX, col)
        #abline(v=x,lwd=lwd1, lty="dashed", col=col)
      }      
    } else {
      y <- Pred.CS[i]
      Nucltip <- Nmrnucl.CS[i]
      if(Nucltip=="CE"){Neighb.nucltip <- "HE1"}
      if(Nucltip=="CB"){Neighb.nucltip <- "HB1"}
      if(Nucltip=="CG1"){Neighb.nucltip <- "HG11"}
      if(Nucltip=="CG2"){Neighb.nucltip <- "HG21"}
      if(Nucltip=="CD"){Neighb.nucltip <- "HD1"}
      if(Nucltip=="CD1"){Neighb.nucltip <- "HD11"}
      if(Nucltip=="CD2"){Neighb.nucltip <- "HD21"}
      Cor.nucl <- "H"
      Cor.CS.ind <- which(Ncl.type==Cor.nucl & 
                          Amacid.CS==Amacid.CS[i] &
                          Resseq.CS==Resseq.CS[i] &
                          Chainid.CS==Chainid.CS[i] &
                          Nmrnucl.CS==Neighb.nucltip)
      if(length(Cor.CS.ind)==1) {
        x <- Pred.CS[Cor.CS.ind]
        if(!is.na(x)){
          lab <- paste(Amacid.CS[i],"-",newlabel(Nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
          pointX <- c(pointX,x)
          pointY <- c(pointY,y)
          pointLAB <- c(pointLAB, lab)
          colorsP <- c(colorsP, col)
          #points(x=x, y=y,col=col, cex=cex1, lwd=lwd1)
          #points(x=x, y=y,col=col, cex=cex2, lwd=lwd2)
          #points(x=x, y=y,col=col, cex=cex3, lwd=lwd3)
          #textxy(X=x, Y=y, labs=paste(Amacid.CS[i],"-",Nucltip,"-",Resseq.CS[i],sep="")  )
        } else {
          lab <- paste(Amacid.CS[i],"-",newlabel(Nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
          lineY <- c(lineY,y)
          lineLAB.Y <- c(lineLAB.Y, lab)
          colorsLY <- c(colorsLY, col)
          #abline(h=y,lwd=lwd1, lty="dashed", col=col)
        }
      } else {
        lab <- paste(Amacid.CS[i],"-",newlabel(Nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
        lineY <- c(lineY,y)
        lineLAB.Y <- c(lineLAB.Y, lab)
        colorsLY <- c(colorsLY, col)
        #abline(h=y,lwd=lwd1, lty="dashed", col=col)
      }      


    }
    ######################
  }
};rm(i)


if(interactive!=TRUE) {
  if(viewplot==TRUE) {
   x11(width=13, height=15)
  } else {
    if(saveplot==TRUE){
     if(x11lib==TRUE){
      jpeg(quality=100, height=plotheight, width=plotwidth, filename=plotname)
     }
    }
  }
} else {
  x11(width=13, height=15)
}

xlimA <- xlim; rm(xlim)
ylimA <- ylim; rm(ylim)

plottingfun <- function(xlim=xlimA, ylim=ylimA, cx=0.5){  
  plot(NA, xlab="1H chemical shift (ppm)",
           ylab="13C chemical shift (ppm)", xlim=xlim, ylim=ylim)
  points(x=pointX, y=pointY, col=colorsP, cex=cex1, lwd=lwd1)
  points(x=pointX, y=pointY, col=colorsP, cex=cex2, lwd=lwd2)
  points(x=pointX, y=pointY, col=colorsP, cex=cex3, lwd=lwd3)
  if(pointlabel==TRUE){
    textxy(X=pointX, Y=pointY, labs=pointLAB, dcol=colorsP, cx=cx, m=c(-2,-2))
  }
  for(i in 1:length(lineX)){
     abline(v=lineX[i],lwd=lwd1, lty="dashed", col=colorsLX[i])
   if(linelabel==TRUE){
     textxy(X=lineX[i], Y=(ylimA[1]+0.4), labs=lineLAB.X[i], dcol=colorsLX[i], cx=cx, m=c(-2,-2))
   }
  };rm(i)
  for(i in 1:length(lineY)){
     abline(h=lineY[i],lwd=lwd1, lty="dashed", col=colorsLY[i])
   if(linelabel==TRUE){
     textxy(X=(xlimA[1]+0.4), Y=lineY[i], labs=lineLAB.Y[i], dcol=colorsLY[i], cx=cx, m=c(-2,-2))
   }
  };rm(i)

  pointX.p <- pointX  # back-up-ing for interactive mode
  pointY.p <- pointY
  pointLAB.p <- pointLAB

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  if(length(experdata)!=0 & plotexper==TRUE){

  lineX <- lineY <- pointX <- pointY <- pointLAB <- lineLAB.X <- 
  lineLAB.Y <- colorsP <- colorsLX <- colorsLY <-  NULL
  for(i in 1:length(Exp.CS)) {
    if(!is.na(Exp.CS[i])) {           # Exp.CS[i] is not NA
      eval(parse(text=paste("col <- col",Amacid.CS[i],"exp",sep="")))
      ######################
      if(Ncl.type[i]=="H") {
        x <- Exp.CS[i]
        Nucltip <- Nmrnucl.CS[i]
        if(Nucltip=="HE1"){Neighb.nucltip <- "CE"}
        if(Nucltip=="HB1"){Neighb.nucltip <- "CB"}
        if(Nucltip=="HG11"){Neighb.nucltip <- "CG1"}
        if(Nucltip=="HG21"){Neighb.nucltip <- "CG2"}
        if(Nucltip=="HD1"){Neighb.nucltip <- "CD"}
        if(Nucltip=="HD11"){Neighb.nucltip <- "CD1"}
        if(Nucltip=="HD21"){Neighb.nucltip <- "CD2"}
        Cor.nucl <- "C"
        Cor.CS.ind <- which(Ncl.type==Cor.nucl & 
                            Amacid.CS==Amacid.CS[i] &
                            Resseq.CS==Resseq.CS[i] &
                            Chainid.CS==Chainid.CS[i] &
                            Nmrnucl.CS==Neighb.nucltip)
        if(length(Cor.CS.ind)==1) {
          y <- Exp.CS[Cor.CS.ind]
          if(!is.na(y)){
            lab <- paste(Amacid.CS[i],"-",newlabel(Neighb.nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
            pointX <- c(pointX,x)
            pointY <- c(pointY,y)
            pointLAB <- c(pointLAB, lab)
            colorsP <- c(colorsP, col)
            #points(x=x, y=y,col=col, cex=cex1, lwd=lwd1)
            #points(x=x, y=y,col=col, cex=cex2, lwd=lwd2)
            #points(x=x, y=y,col=col, cex=cex3, lwd=lwd3)
            #textxy(X=x, Y=y, labs=lab)
          } else {
            lab <- paste(Amacid.CS[i],"-",newlabel(Neighb.nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
            lineX <- c(lineX,x)
            lineLAB.X <- c(lineLAB.X, lab)
            colorsLX <- c(colorsLX, col)
            #abline(v=x,lwd=lwd1, lty="dashed", col=col)
          }
        } else {
          lab <- paste(Amacid.CS[i],"-",newlabel(Neighb.nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
          lineX <- c(lineX,x)
          lineLAB.X <- c(lineLAB.X, lab)
          colorsLX <- c(colorsLX, col)
          #abline(v=x,lwd=lwd1, lty="dashed", col=col)
        }      
      } else {
        y <- Exp.CS[i]
        Nucltip <- Nmrnucl.CS[i]
        if(Nucltip=="CE"){Neighb.nucltip <- "HE1"}
        if(Nucltip=="CB"){Neighb.nucltip <- "HB1"}
        if(Nucltip=="CG1"){Neighb.nucltip <- "HG11"}
        if(Nucltip=="CG2"){Neighb.nucltip <- "HG21"}
        if(Nucltip=="CD"){Neighb.nucltip <- "HD1"}
        if(Nucltip=="CD1"){Neighb.nucltip <- "HD11"}
        if(Nucltip=="CD2"){Neighb.nucltip <- "HD21"}
        Cor.nucl <- "H"
        Cor.CS.ind <- which(Ncl.type==Cor.nucl & 
                            Amacid.CS==Amacid.CS[i] &
                            Resseq.CS==Resseq.CS[i] &
                            Chainid.CS==Chainid.CS[i] &
                            Nmrnucl.CS==Neighb.nucltip)
        if(length(Cor.CS.ind)==1) {
          x <- Exp.CS[Cor.CS.ind]
          if(!is.na(x)){
            lab <- paste(Amacid.CS[i],"-",newlabel(Nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
            pointX <- c(pointX,x)
            pointY <- c(pointY,y)
            pointLAB <- c(pointLAB, lab)
            colorsP <- c(colorsP, col)
            #points(x=x, y=y,col=col, cex=cex1, lwd=lwd1)
            #points(x=x, y=y,col=col, cex=cex2, lwd=lwd2)
            #points(x=x, y=y,col=col, cex=cex3, lwd=lwd3)
            #textxy(X=x, Y=y, labs=paste(Amacid.CS[i],"-",Nucltip,"-",Resseq.CS[i],sep="")  )
          } else {
            lab <- paste(Amacid.CS[i],"-",newlabel(Nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
            lineY <- c(lineY,y)
            lineLAB.Y <- c(lineLAB.Y, lab)
            colorsLY <- c(colorsLY, col)
            #abline(h=y,lwd=lwd1, lty="dashed", col=col)
          }
        } else {
          lab <- paste(Amacid.CS[i],"-",newlabel(Nucltip),"-",Resseq.CS[i],"-",Chainid.CS[i],sep="")
          lineY <- c(lineY,y)
          lineLAB.Y <- c(lineLAB.Y, lab)
          colorsLY <- c(colorsLY, col)
          #abline(h=y,lwd=lwd1, lty="dashed", col=col)
        }      


      }
      ######################
    }
  };rm(i)

  if(length(expCOL)!=0) {
    points(x=pointX, y=pointY, col=expCOL, cex=cex1, lwd=lwd1)
    points(x=pointX, y=pointY, col=expCOL, cex=cex2, lwd=lwd2)
    points(x=pointX, y=pointY, col=expCOL, cex=cex3, lwd=lwd3)

    if(pointlabel==TRUE){
      textxy(X=pointX, Y=pointY, labs=pointLAB, dcol=expCOL, cx=cx, m=c(-2,-2))
    }
    for(i in 1:length(lineX)){
       abline(v=lineX[i],lwd=lwd1, lty="dashed", col=expCOL)
     if(linelabel==TRUE){
       textxy(X=lineX[i], Y=(ylimA[1]+0.4), labs=lineLAB.X[i], dcol=expCOL, cx=cx, m=c(-2,-2))
     }
    };rm(i)
    for(i in 1:length(lineY)){
       abline(h=lineY[i],lwd=lwd1, lty="dashed", col=expCOL)
     if(linelabel==TRUE){
       textxy(X=(xlimA[1]+0.4), Y=lineY[i], labs=lineLAB.Y[i], dcol=expCOL, cx=cx, m=c(-2,-2))
     }
    };rm(i)
  } else {
    points(x=pointX, y=pointY, col=colorsP, cex=cex1, lwd=lwd1)
    points(x=pointX, y=pointY, col=colorsP, cex=cex2, lwd=lwd2)
    points(x=pointX, y=pointY, col=colorsP, cex=cex3, lwd=lwd3)
  
    if(pointlabel==TRUE){
      textxy(X=pointX, Y=pointY, labs=pointLAB, dcol=colorsP, cx=cx, m=c(-2,-2))
    }
    for(i in 1:length(lineX)){
       abline(v=lineX[i],lwd=lwd1, lty="dashed", col=colorsLX[i])
     if(linelabel==TRUE){
       textxy(X=lineX[i], Y=(ylimA[1]+0.4), labs=lineLAB.X[i], dcol=colorsLX[i], cx=cx, m=c(-2,-2))
     }
    };rm(i)
    for(i in 1:length(lineY)){
       abline(h=lineY[i],lwd=lwd1, lty="dashed", col=colorsLY[i])
     if(linelabel==TRUE){
       textxy(X=(xlimA[1]+0.4), Y=lineY[i], labs=lineLAB.Y[i], dcol=colorsLY[i], cx=cx, m=c(-2,-2))
     }
    };rm(i)
   }
  }
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # results for interactive mode
  RESL <- NULL
  RESL$pointX.p <- pointX.p
  RESL$pointY.p <- pointY.p
  RESL$pointLAB.p <- pointLAB.p
  return(RESL)
}



if(zooming == TRUE){
  zoom.rev(fun = plottingfun, zoom.col = "red", delay = 0.5, xlimINI=xlimA, ylimINI=ylimA, cx=0.6)
} else {
  RESL <- plottingfun(xlim=xlimA, ylim=ylimA, cx=0.6)
}

if(interactive!=TRUE & viewplot!=TRUE & saveplot==TRUE & zooming!=TRUE) { dev.off() }

if(interactive==TRUE & zooming!=TRUE){
  print("NOTE: Now you can annotate the prediction results on the plot")
  print("NOTE:  by left-clicking on the predicted spots. The plot window")
  print("NOTE:  should not be resized and the data cannot be saved!")
  if(length(experdata)!=0 & plotexper==TRUE) {
    identify(x=c(RESL$pointX.p, pointX), y=c(RESL$pointY.p, pointY), label=c(RESL$pointLAB.p, pointLAB))
  } else {
    identify(x=RESL$pointX.p, y=RESL$pointY.p, label=RESL$pointLAB.p)
  }
}

#write(paste(outputname,",",plotname, sep=""), file="resultfiles.txt")
print("NOTE: CH3Shift is DONE!!!")
write("NOTE: CH3Shift is DONE!!!", file="process_info.txt", append=TRUE)

############ FINISHED
