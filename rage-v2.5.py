import argparse
import string
import sys
import subprocess
import os
import time
import convertgff3
import renameSeqnamegff3
import assignscore
import Combiner
import callPeaks
#nohup python rage-v2.5.py $HOME/results/2014-04-07/rnaseqtrimmed.fastq $HOME/results/data/outputsfull/filter/ribotrimmedfiltered_rRNA_tRNA.fastq $HOME/testdata/k12genome_newheader.fasta -c false -f false -s ecolik12 -d $HOME/results/data/outputfullwglimmer/ &
#python rage-v2.5.py $HOME/testdata/k12rnaseq-test-3K.fastq $HOME/testdata/k12riboseq-test-3K.fastq $HOME/testdata/k12genome_newheader.fasta -c false -f false
#nohup python rage-v2.5.py $HOME/results/data/outputfullwglimmer/bowtie/rnaseqtrimmed.sorted.bam $HOME/results/data/outputsfull/bowtie/ribotrimmedfiltered_rRNA_tRNA.sam $HOME/testdata/k12genome.fasta -c false -f false -b false -s ecolik12 -d $HOME/results/data/outputsfinal0/ &

def adaptlist():
    with open("../data/adapterlist.txt", "r") as f:
        adapts = f.read()
    adapts = adapts.strip()
    adaptlist = adapts.split('\n')
    return adaptlist

class Pipeline(object):
    def __init__(self, rnainput, riboinput, genomeinput,cutadapt,filter_RNA,bowtie,WorkDirctory,
                 rna_adapter, ribosome_adapter,weights,species):
        self.output = WorkDirctory
        self.rRNAfilterfile = "../data/rRNA_data/SSURef_NR99_115_tax_silva.fasta"
        self.tRNAfilterfile = "../data/tRNA_data/bacterial-tRNAs.fa"
        self.rnainput = rnainput
        self.riboinput = riboinput
        self.genomeinput = genomeinput
        self.cutadapt = cutadapt
        self.filter_RNA = filter_RNA
        self.bowtie = bowtie
        self.rna_adapter = rna_adapter
        self.ribosome_adapter = ribosome_adapter
        self.rnatag = "RNA-seq"
        self.ribotag = "Ribo-footprinting"
        
        self.weightSetting = weights.split("#")
        self.weightPeriod = self.weightSetting[0]
        self.weightRNA = self.weightSetting[1]
        self.weightNaive = self.weightSetting[2]
        
        self.seqname_fasta = species + "_genome"
        self.seqname_glimmer = species + "_glimmer"
        self.seqname_cem = species + "_cem"
        self.seqname_naive = species + "_naiveORFs"
        self.seqname_periodicity = species +"Ribo-footprinting"
        
        self.naivepred_output = self.output + "naive-pred/"+ os.path.basename(self.genomeinput)[:-len(".fasta")]+".predict"
        
        self.cutadapt_setting = ["-e", "0.1", "-O", "5", "-m", "15", "-q", "20", "-z"]
        self.cutadapt_rnaoutput = (self.output + "cutadapt/"
        + os.path.basename(self.rnainput) [:-len(".fastq")] + "-RNAtrimmed.fastq")
        self.cutadapt_rnareport = (self.output + "cutadapt/"
        + os.path.basename(self.rnainput)[:-len(".fastq")] + "-RNAreport.txt")
        self.cutadapt_ribooutput = (self.output + "cutadapt/"
        + os.path.basename(self.riboinput) [:-len(".fastq")] + "-ribotrimmed.fastq")
        self.cutadapt_riboreport = (self.output + "cutadapt/"
        + os.path.basename(self.riboinput)[:-len(".fastq")] + "-riboreport.txt")
        
 
    def cutadapter (self, cutadapt_args, inputpath, cutadapt_output,cutadapt_report, adapters ):
        cutadapt_args.append(inputpath)
        #set up RNA output
        cutadapt_args.append("-o")
        cutadapt_args.append(cutadapt_output)
        #set up RNA report
        cutadapt_args.append("1>")
        cutadapt_args.append(cutadapt_report)
        #set up RNA adapter
        for adapter in adapters:
            cutadapt_args.append("-a")
            cutadapt_args.append(adapter)
        #set up setting
        for setting in self.cutadapt_setting:
            cutadapt_args.append(setting)
        #commandline generation
        for i in xrange(len(cutadapt_args)):
            cutadapt_args[i] = cutadapt_args[i] + " "
        cutadapt_args = "".join(cutadapt_args)
        #call cutadapt for RNA
        subprocess.Popen(cutadapt_args, shell = True)
        
    def cutadaptWrapper(self):
        if not(os.path.isdir(self.output)):
            mkDirArgs=["mkdir", self.output]
            subprocess.call(mkDirArgs)
        if not(os.path.isdir(self.output+"cutadapt/")):
            mkDirArgs=["mkdir", self.output+"cutadapt/"]
            subprocess.call(mkDirArgs)
        if self.cutadapt == "false":
            print "skipping cutadapt..."
            self.cutadapt_rnaoutput = self.rnainput
            self.cutadapt_ribooutput = self.riboinput
        else:
            cutadapt_rnaargs = ["cutadapt"]
            print "cutadapter for RNA-seq..."
            self.cutadapter(cutadapt_rnaargs, self.rnainput, self.cutadapt_rnaoutput ,self.cutadapt_rnareport, self.rna_adapter)
            print "cutadapter for ribosomal-seq..."
            cutadapt_riboargs = ["cutadapt"]
            self.cutadapter(cutadapt_riboargs, self.riboinput, self.cutadapt_ribooutput ,self.cutadapt_riboreport, self.ribosome_adapter)
        
    def filter_rRNA_tRNA(self):
        
        if not(os.path.isdir(self.output+"filter/")):
            subprocess.call(["mkdir",self.output+"filter/"])
        #if cutadapt is not done, wait
        if self.cutadapt == "true":
            while (not(os.path.isfile(self.cutadapt_riboreport)) or os.path.getsize(self.cutadapt_riboreport) <= 5 ):
                time.sleep(2)

        self.filter_input = self.cutadapt_ribooutput
        if self.filter_RNA == "false":
            print "skipping filter..."
            self.filter_output = self.filter_input
        else:
            print "Filtering rRNA and tRNA...."
            self.filter_mid = self.output+"filter/"+os.path.basename(self.filter_input)[0:-len(".fastq")]+"filtered_rRNA.fastq"
            self.filter_output = self.output+"filter/"+os.path.basename(self.filter_input)[0:-len(".fastq")]+"filtered_rRNA_tRNA.fastq"        
            #build the index for rRNA dataset
            index_rRNA_Output = self.output+"filter/"+os.path.basename(self.rRNAfilterfile)[:-len(".fasta")]+"index"
            if not(os.path.isfile(index_rRNA_Output+".1.ebwt")):
                buildIndexArgs_rRNA = ["bowtie-build",self.rRNAfilterfile,index_rRNA_Output]
                subprocess.call(buildIndexArgs_rRNA)
            #filter the reads and get the unaligned reads
            filter_rRNA_Args=["bowtie","-S","--un",self.filter_mid,index_rRNA_Output,self.filter_input]
            subprocess.call(filter_rRNA_Args)

            #build the index for tRNA dataset
            index_tRNA_Output=self.output+"filter/"+os.path.basename(self.tRNAfilterfile)[0:-len(".fa")]+"index"
            if not(os.path.isfile(index_tRNA_Output+".1.ebwt")):
                buildIndexArgs_tRNA=["bowtie-build",self.tRNAfilterfile,index_tRNA_Output]
                subprocess.call(buildIndexArgs_tRNA)
   
            #filter the reads and get the unaligned reads
            filter_tRNA_Args=["bowtie","-S","--un",self.filter_output,index_tRNA_Output,self.filter_mid]
            subprocess.call(filter_tRNA_Args)

    def bowtieAlign(self, genomeData, seqData):
        #building genome index
        indexOutput=self.output+"bowtie/"+os.path.basename(genomeData)[0:-len(".fasta")]+"index"
        buildIndexArgs=["bowtie-build", genomeData, indexOutput]
        subprocess.call(buildIndexArgs)
        
        #align reads to genome    
        alignOutput=self.output+"bowtie/"+os.path.basename(seqData)[0:-len(".fastq")]+".sam"
        alignArgs=["bowtie", "-S", indexOutput, seqData, alignOutput]
        subprocess.call(alignArgs)

        #convert to bam using samtools
        bamOutput=alignOutput[0:len(alignOutput)-len(".sam")]+".bam"
        bamArgs=["samtools", "view", "-bS", "-o", bamOutput, alignOutput]
        subprocess.call(bamArgs)

        #sort bam
        sortedOutput=bamOutput[0:len(bamOutput)-len(".bam")]+".sorted"
        sortArgs=["samtools", "sort", bamOutput, sortedOutput]
        subprocess.call(sortArgs)

        #index sorted bam file
        sortedOutput+=".bam"
        
        indexSortArgs=["samtools", "index", sortedOutput]
        subprocess.call(indexSortArgs)

        #remove all unnecesary files
        for filename in os.listdir(self.output+"bowtie/"):
            if ".ebwt" in filename:
                os.remove(self.output+"bowtie/"+filename)
        return sortedOutput
                
    def bowtieAlignWrapper(self):
        if self.bowtie == 'false':
            print "skipping Bowtie..."
            self.bowtie_rnaoutput = self.rnainput
            self.bowtie_ribooutput = self.riboinput
        else:
            self.bowtie_rnainput = self.cutadapt_rnaoutput
            self.bowtie_riboinput = self.filter_output
            if not(os.path.isdir(self.output+"bowtie/")):
                mkDirArgs=["mkdir", self.output+"bowtie/"]
                subprocess.call(mkDirArgs)
            print "Mapping to genome..."
            self.bowtie_rnaoutput = self.bowtieAlign(self.genomeinput, self.bowtie_rnainput)
            self.bowtie_ribooutput = self.bowtieAlign(self.genomeinput, self.bowtie_riboinput)
            self.bowtie_ribooutput = self.bowtie_ribooutput[:-len(".sorted.bam")]+".sam"
        
    
    def CEMassembly(self, filepath):
        subprocess.call(["runcem.py", filepath])
    
    def CEMwrapper(self):
        print "Assembling RNAreads..."
        self.cem_rnainput = self.bowtie_rnaoutput
        self.cem_rnaoutput = self.output + "bowtie/" + os.path.basename(self.cem_rnainput)+ ".instance.pred.gtf"
        self.CEMassembly(self.cem_rnainput)

    
    def gtf2gff3(self, filepath,seqtype):
        tools = "gffread"
        outputfile = self.output + "transcripts/" + os.path.basename(filepath)[0:-len(".gtf")]+"mid.gff3"
        Args = [tools,filepath,"-t",seqtype, "-o", outputfile]
        subprocess.call(Args)
    def reSeqname(self,filename,outputName,seqname):
        files = renameSeqnamegff3.renameSeq(filename,outputName,seqname)
        files.renameSeqnameGFF3()
    
    def gtf2gff3wrapper(self):
        self.gtf2gff3_rnainput = self.cem_rnaoutput
        self.gtf2gff3_rnaoutput = self.output + "transcripts/" + os.path.basename(self.gtf2gff3_rnainput)[0:-len(".gtf")]+".gff3"
        self.gtf2gff3_rnamid = self.output + "transcripts/" + os.path.basename(self.gtf2gff3_rnainput)[0:-len(".gtf")]+"mid.gff3"
        print "converting gtf to gff3..."
        if not(os.path.isdir(self.output+"transcripts/")):
            subprocess.call(["mkdir",self.output+"transcripts/"])
        while not(os.path.isfile(self.gtf2gff3_rnainput)):
            time.sleep(2)
        self.gtf2gff3(self.gtf2gff3_rnainput,self.rnatag)
        self.reSeqname(self.gtf2gff3_rnamid,self.gtf2gff3_rnaoutput,self.seqname_cem+self.rnatag )
        
    def naiveORFfind(self):
        print "Naive long ORFs finding..."
        if not(os.path.isdir(self.output + "naive-pred/")):
            mkDirArgs=["mkdir", self.output + "naive-pred/"]
            subprocess.call(mkDirArgs)
        predArgs=["long-orfs","-n" ,"-t", "1.15", self.genomeinput, self.naivepred_output]
        subprocess.call(predArgs)
        self.naivepred_output2gff3 = self.naivepred_output[:-len(".predict")]+".gff3"
        self.predict2gff3(self.naivepred_output, self.naivepred_output2gff3,self.seqname_naive)


#####
#
# We'll need pointer to our periodicity verified training data
# and stuff here
    def predict2fasta(self, genomeData, inFile,outFile):
        subprocess.Popen("extract "+genomeData+" "+inFile+" > "+outFile,shell = True)
        
    def callPeaks(self):
        if not(os.path.isdir(self.output+"peaks/")):
            mkDirArgs=["mkdir", self.output+"peaks/"]
            subprocess.call(mkDirArgs)
        
        print "Analyzing ribo-footprinting data..."
        self.periodicity_input = self.bowtie_ribooutput
        self.periodicity_output = self.output+"peaks/"+os.path.basename(self.bowtie_ribooutput)[0:-len(".sam")]+"-PeaksTest.predict"
        callPeaks.run( self.periodicity_input, self.genomeinput,self.periodicity_output)
        self.periodicity_output2gff3 = self.periodicity_output[:-len(".predict")]+".gff3"
        self.periodicity_output2fasta = self.periodicity_output[:-len(".predict")]+".fasta"
        self.predict2gff3(self.periodicity_output, self.periodicity_output2gff3,self.seqname_periodicity)
        self.predict2fasta(self.genomeinput, self.periodicity_output, self.periodicity_output2fasta)
    
        
#
#####

    def icmBuild(self, inFile, outFile):
        subprocess.Popen("build-icm " + outFile + " < "+ inFile,shell = True)
    
    def glimmerRun(self, genomeData, icmFile, identifier):
        subprocess.call(["glimmer3",genomeData,icmFile,identifier])


    def predict2gff3(self, filename,outputName,seqname):
        convertgff3.convertToGFF3(filename,outputName,seqname)

    def glimmerWrapper(self):
        if not(os.path.isdir(self.output + "glimmer/")):
            mkDirArgs=["mkdir", self.output + "glimmer/"]
            subprocess.call(mkDirArgs)
        print "Constructing Glimmer Training Data..."
        self.glimmer_icm_input  = self.periodicity_output2fasta #= self. This needs to be our training data
        tag = self.output+"glimmer/"+os.path.basename(self.glimmer_icm_input)[:-len(".fasta")]
        self.glimmer_icm  = tag +".icm" #= self. This should be our traning data with *.icm
        self.icmBuild(self.glimmer_icm_input, self.glimmer_icm)
        while (not(os.path.isfile(self.glimmer_icm)) or os.path.getsize(self.glimmer_icm) <= os.path.getsize(self.glimmer_icm_input)/4 ):
            time.sleep(2)
        ## We need to set a tag here
        print "Running Glimmer..."
        self.glimmerRun(self.genomeinput, self.glimmer_icm, tag )
        # The output files are going to be tag.predict and tag.detail
        self.glimmer_predict  = tag+".predict"#= our tag.predict file
        self.glimmer_predict2gff3 = tag+".gff3"
        self.predict2gff3(self.glimmer_predict, self.glimmer_predict2gff3,self.seqname_glimmer)

    def Combiner(self):
        print "Assigning scores to evidence..."
        assignscore.AssignScore(self.periodicity_output2gff3,self.weightPeriod)
        assignscore.AssignScore(self.naivepred_output2gff3,self.weightNaive)
        assignscore.AssignScore(self.gtf2gff3_rnaoutput,self.weightRNA)
        print "Combining evidence..."
        results = Combiner.combine(self.periodicity_output2gff3,self.naivepred_output2gff3,self.gtf2gff3_rnaoutput,self.glimmer_predict2gff3)
        self.outputgff3 = self.genomeinput[:-len(".fasta")]+".gff3"
        print "Writing results..."
        assignscore.writeResults2gff3(results,self.outputgff3,self.seqname_fasta,self.seqname_Combiner)
        
        
#####
#
#####
    
parser = argparse.ArgumentParser(description='Process adapter trimming.')
parser.add_argument('rnainput', type=str , help='input RNA-seq fastq file')
parser.add_argument('riboinput', type=str , help='input Ribosome fastq file')
parser.add_argument('genomeinput',type=str , help='input Genome fasta file')
parser.add_argument('-c','--cutadapt',type=str , choices=['true', 'false'],default = 'true',help='whether use cutadapt,\
                    if choose false, please input the clean RNA-seq and ribosome fastq file; default is true')
parser.add_argument('-f','--filter_rRNA_tRNA',type=str , choices=['true', 'false'],default = 'true',help='whether filter rRNA and tRNA,\
                    if choose false, please input the clean ribosome fastq file; default is true')
parser.add_argument('-b','--bowtie',type=str , choices=['true', 'false'],default = 'true',help='whether use bowtie to map,\
                    if choose false, please input sam file for Ribo-footprinting and sorted .bam for RNA-seq; default is true')
parser.add_argument('-d','--WorkDirctory', type=str , default = "../results/data/outputs/", help='working directory,\
                    default in results/data/outputs, if user specified, please use full path.')
parser.add_argument('-a', "--rna_adapter", type = str, action = 'append', default = adaptlist(),
                    help='input custom adapt sequence for RNA-seq')
parser.add_argument('-i', "--ribosome_adapter", type = str, action = 'append',
        default = ['CTGTAGGCACCATCAAT'], help='input custom adapt sequence for Ribosome profiling')
parser.add_argument('-w','--weights',type = str, default = "10#5#1",help='Assign points for ribo-footprinting, RNAtranscripts and naiveORFs; format XX#XX#XX')
parser.add_argument('-s','--species',type = str,default = 'speciesX', help='Species name of the data')

 
args = parser.parse_args()
pipe = Pipeline(args.rnainput, args. riboinput, args. genomeinput,args.cutadapt,
                args.filter_rRNA_tRNA,args.bowtie,args.WorkDirctory, args.rna_adapter, args.ribosome_adapter,
                args.weights,args.species)
pipe.cutadaptWrapper()
pipe.filter_rRNA_tRNA()
pipe.naiveORFfind()
pipe.bowtieAlignWrapper()
pipe.callPeaks()
pipe.CEMwrapper()
pipe.gtf2gff3wrapper()

pipe.glimmerWrapper()
pipe.Combiner()
