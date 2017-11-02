## List : 
## 	1. Function to get choped sequences
## 	2. Read a conting file in FASTA formate 
##	3. Function to get unique sequence FASTA file with frequency count
##	4. Function to handle pairwise identy matrix of BioEdit software output
## 	5. Function to get the matching part from main sequence based on BLAST output
##	6. Clean a specific character form FASTA file 
##	7. Extract genome assemble name from fasta file 
##	8. Find duplicate gene in genome


##############################################################
## 1. Function to get choped sequences
##############################################################
## file.fasta = FASTA formated file, form = Amino acid of neocleotide start position , to = form = Amino acid of neocleotide stop position, file.out= Output 
chompSeq<-function(file.fasta,from,to,file.out = "result.fasta")
	{	sequ.names <- names(file.fasta)
		sequ <- NULL
		sequ[sequ.names] <- list(NULL)
		for(i in 1 : length(names(file.fasta)))
		{
			sequ[[i]]<-substr(file.fasta[[i]][1], from,to)
		}	
		write.fasta(sequences = sequ, names = names(file.fasta), nbchar = 80, file.out)
	}
##############################################################

###############################################################	
## 2. Read a conting file in FASTA formate 
##############################################################
## Pipeline : 	1. Read a part of a sequence in a FASTA file 
##				2. Translate it to protein sequence of six reading frame 
##				3. Return biggest protein sequence form six reading frame by consider STOP codon.
## Input : contig= ID of interest, aesembly= assemble of interest, file = Main contig file in FASTA formate, form = start, to = stop, 
## Example : 	1. Sample input file header : 
##				>cont1.62 dna:supercontig supercontig:GCA_000174475.1:cont1.62:1:32906:1
##				AAATGCCAGGGGGATTTTCHATGCAT
##				2. contig = "cont1.62", aesembly= "aesembly", file ="sample.fasta",from=116271,to=122902,strand=1
cont2prot<-function(contig="adpmb-supercont1.5",aesembly="asm",file=fast,from=116271,to=122902,strand=1)
	{	mat<-c(contig,aesembly)
		mine <- c(1,1)
		nam<-unlist(lapply(getAnnot(file),function(i)strsplit(i,">")[[1]][2]))
		fast<-file[which(apply(apply(sapply(mat, regexpr, nam, ignore.case=TRUE),1,function(i)as.numeric(i>-1)),2,identical,mine)=="TRUE")]
		att<-attributes(fast[[1]])
		att2<-att[[2]]
		sequ<-substr((fast[[1]][1]),from,to)
			
		strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
			
		sequ1<-strReverse(sequ)
		contig_strand<-as.numeric(as.character(substr(att2,nchar(att2),nchar(att2))))
			
		if((contig_strand==strand)==TRUE)
			{
			f1<-paste(getTrans(s2c(sequ),frame=0,sens="F"),collapse = '')
			f1_dna<-sequ
			f2<-paste(getTrans(s2c(sequ),frame=1,sens="F"),collapse = '')
			f2_dna<-substr(sequ,2,nchar(sequ))
			f3<-paste(getTrans(s2c(sequ),frame=2,sens="F"),collapse = '')
			f3_dna<-substr(sequ,3,nchar(sequ))
				
			r1<-paste(getTrans(s2c(sequ),frame=0,sens="R"),collapse = '')
			r1_dna<-paste(comp(s2c(sequ1),forceToLower = FALSE),collapse = '')
			r2<-paste(getTrans(s2c(sequ),frame=1,sens="R"),collapse = '')
			r2_dna<-paste(comp(s2c(substr(sequ1,2,nchar(sequ1))),forceToLower = FALSE),collapse = '')
			r3<-paste(getTrans(s2c(sequ),frame=2,sens="R"),collapse = '')
			r3_dna<-paste(comp(s2c(substr(sequ1,3,nchar(sequ1))),,forceToLower = FALSE),collapse = '')
				
			res<-list(f1,f1_dna,f2,f2_dna,f3,f3_dna,r1,r1_dna,r2,r2_dna,r3,r3_dna)
			names(res)<-c("for_1","for_1_DNA","for_2","for_2_DNA","for_3","for_3_DNA","rev_1","rev_1_DNA","rev_2","rev_2_DNA","rev_3","rev_3_DNA")
				
			res1<-list(f1,f2,f3,r1,r2,r3)
			res2<-r1
				
			## Checking which strand has max length
			se<-res1
			ind<-which(rapply(rapply(lapply(se,strsplit,"\\*"),nchar, how = "list"),max)==max(rapply(rapply(lapply(se,strsplit,"\\*"),nchar, how = "list"),max)))
			res3<-res1[[ind]]
			return(res3)
			}
		if ((contig_strand==strand)==FALSE)
			{note<-"CHECK"
			return(note)}
	}
###############################################################	

###############################################################	
## 3. Function to get unique sequence FASTA file with frequency count
##############################################################
## Input : 	all.fa=  Input fasta file,
## Output : Example: 
##			>unique.1_2529
##			MKNNLRYGIRKHKLGAASVFLGTMIVVGMGQDKEAA
## Explanation : Unique sequence 1 with frequency 2519
useq<-function(all.fa,file.out="result.fasta")
	{	all_seq<-unlist(all.fa)
		uniq_seq<-unique(unlist(all.fa))
		nam_uniq_seq<-paste("unique",1:length(uniq_seq),sep=".")
	
		fre_uniq_seq<-rep(NA,length(uniq_seq))
		lst_uniq_seq<-list(NULL)
	
		for(i in 1: length(uniq_seq) )
		{	lst_uniq_seq[[i]]	<-	names(all_seq[which(all_seq%in%uniq_seq[i])])
			fre_uniq_seq[i]		<- 	length(which(all_seq%in%uniq_seq[i]))
		}
		nam<-paste(nam_uniq_seq,fre_uniq_seq,sep="_")
		write.fasta(sequences = as.list(uniq_seq),names = nam, nbchar=60,as.string = TRUE, file.out)	
	}
###############################################################	

###############################################################	
## 4. Function to handle pairwise identy matrix of BioEdit software output
##############################################################
##	Input : data: BioEdit output of Identity matrix (Tab delimated file )
##			alignment= 	IF the input file if alligned or not. 
##						If align , BioEdit will produce a concensus sequence column
##			readLines commend can be used to read the input file 
## 			Example : fnbpa_N1_u_idt<-readLines("fnbpa_N1_u.idt")
## 	Output: Mean,Variance, Quantile of upper diagonal matrix
##	Example : 	data<-readLines("FileName")
##				result<- iden2disp(data)
iden2disp<-function(data,alignment=TRUE)
	{	data<-data[-(1:3)]
		tmp<-matrix(NA,ncol=length(data),nrow=length(data))
		for (i in 1:length(data)) tmp[i,]<-unlist(strsplit(data[i],split="\t"))
		tmp[which(tmp=="ID")]<-1

		if(alignment==TRUE)
			{	cnam<-tmp[1,][-c(1,nrow(tmp))] 	## Name of Column
				rnam<-tmp[,1][-c(1,nrow(tmp))]	## Name of Row 
				data1<-tmp[,-c(1,nrow(tmp))]
				data2<-data1[-c(1,nrow(tmp)),]
				colnames(data2)<-cnam
				rownames(data2)<-rnam
			}
		else{	cnam<-tmp[1,][-1] 				## Name of Column
				rnam<-tmp[,1][-1]				## Name of Row 
				data1<-tmp[,-1]
				data2<-data1[-1,]
				colnames(data2)<-cnam
				rownames(data2)<-rnam
			}

		data2[lower.tri(data2,diag=TRUE)] <- NA
		val<-((data2[which(is.na(data2)==FALSE)]))
		library(stringr); val<-as.numeric(str_replace(val,",","."))
		res<-c(n=length(cnam),pair=length(val),min=min(val),max=max(val),mean=mean(val),std=sd(val),quantile(val))
		res1<-list(val,res)
		names(res1)<-c("pairwise_identity","summary")
		return(res1)
	}
##############################################################

##############################################################
## 5. Function to get the matching part from main sequence based on BLAST output
##############################################################
##	Input :  	file.fasta: Main FASTA formated file 
##				file.blast: BLAST output in tabular (-outfmt 6) formate.
##				reverse: Defalt is 	FALSE 	- The match length in BLAST output will reture
##									TRUE 	- The rest of the sequence except the matching part wiill retuen
## Output : FASTA formated file  
## Limitation : reverse option "TRUE" will work if only BLAST match occure at the starting of the sequence.
blast2seq<-function(file.fasta, file.blast,reverse=FALSE,file.out = "result.fasta")
	{	pos<-match(names(file.fasta),file.blast$subject.id)
	
		sequ.names <- names(file.fasta)
		sequ <- NULL
		sequ[sequ.names] <- list(NULL)
		
		if(reverse==FALSE)
			{	for(i in 1 : length(names(file.fasta)))
					{	sequ[[i]]<-substr(file.fasta[[i]][1],as.numeric(as.character(file.blast$s.start[pos[i]])),as.numeric(as.character(file.blast$s.end[pos[i]])))
					}
			}
		
		if(reverse==TRUE)
			{	#pos<-match(file.blast$subject.id,names(file.fasta))
				leng<-unlist(lapply(file.fasta,function(i)nchar(i[[1]][1])))[match(file.blast$subject.id,names(file.fasta))]
				file.blast<-cbind(file.blast,long=leng)
		
				for(i in 1 : length(names(file.fasta)))
					{	sequ[[i]]<-substr(file.fasta[[i]][1],(as.numeric(as.character(file.blast$s.end[pos[i]])))+1,as.numeric(as.character(file.blast$long[pos[i]])))
					}
			}
		
		write.fasta(sequences = sequ, names = names(file.fasta), nbchar = 80,as.string = TRUE, file.out)
	}

##############################################################
## 6. Clean a specific character form FASTA file 
##############################################################
##
##	file.input = File name without quote (FASTA file formate)
## 	pattern = "-"
##	file.out = "filename.fa"
clnFasta<-function(file.input,pattern="-",file.out)
	{	#gg<-(lapply(file.input,function(i)strsplit(i,ele)[[1]][[1]][1]))
		gg<-lapply(file.input,function(i) gsub(pattern,"",i))
		write.fasta(sequences = gg,names = names(gg), nbchar=80, file.out, as.string = TRUE)
	}

##############################################################
## 7. Extract genome assemble name from fasta file 
##############################################################
##	Input : FASTA formated file
##############################################################
str2genome<-function(a)
	{	nam<-toString(attributes(a)[[1]])
		anno<-toString(attributes(a)[[2]])
		library(stringr)
		p<-str_extract(anno,":[A-Z]{1,9}_[0-9|.|0-9]{1,20}:")
		asmbl<-strsplit(p,":")[[1]][2]
		asmbl_1<-lapply(asmbl,function(i)strsplit(i,".",fixed=TRUE)[[1]][1])
		re<-unlist(c(nam,asmbl_1))
		return(re)
	}

##############################################################
## 8. Find duplicate gene in genome
##############################################################
## 	input : FASTA file
##			Header format : >ADC37925 pep:novel supercontig:GCA_000025145.2:CP001844:1841991:1848551:-1 gene:SA2981_1714 transcript:ADC37925 description:"Predicted cell-wall-anchored protein SasC (LPXTG motif)"
## Output : $GCA_000536175.1
##			[1] "EUP12568" "EUP14795" 
##############################################################
dupgnm<-function(a)
			{
			str2genome<-function(a)
				{	nam<-toString(attributes(a)[[1]])
					anno<-toString(attributes(a)[[2]])
					library(stringr)
					p<-str_extract(anno,":[A-Z]{1,9}_[0-9|.|0-9]{1,20}:")
					asmbl<-strsplit(p,":")[[1]][2]
					re<-c(nam,asmbl)
					return(asmbl)
				}	
			genm<-lapply(a,str2genome)
			kk<-unlist(genm)
			ta<-table(kk)
			gg<-ta[which(ta>0)]		
			mm<-kk[which(kk%in%names(gg))]
			
			li<-list()
			for(i in 1: length(unique(mm)))
				{	li[i]<-list(names(mm)[which(mm%in%unique(mm)[i])])
					names(li)[i]<-unique(mm)[i]
				}
			return(li)
			}
