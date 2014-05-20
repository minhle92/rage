import subprocess
import sys

def pathhelp():
    path = []
    s = raw_input('Please input full path for where you put our folder:')
    if s[-1] !="/":
        s = s+"/"
    path0 = "export PATH=$PATH:/"
    path1 = ["tools/cutadapt-1.4.2/bin","tools/glimmer3.02/bin","tools/cufflinks-2.2.0.Linux_x86_64 2","tools/cem/bin","tools/bowtie-1.0.1","tools/samtools-0.1.19"]
    p=raw_input("\nPlease enter the path for your system's bash_profile file:")
    if p[-1]!="/":
       p=p+"/"
    fo=open(p+".bash_profile","a+") 
    path_bash=p+"~/.bash_profile"
    for item in path1:
         if not path0+s+item in fo.readlines():
              subprocess.call(["echo","path0+s+item",">>","path_bash")
    fo.close()
    
    
pathhelp()
