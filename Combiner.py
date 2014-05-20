from operator import itemgetter, attrgetter
import sys

def mycmp(l1, l2):
    # compare s1 and s2
    # result > 0   ==>  s1 > s2
    # result == 0  ==>  s1 == s2
    # result < 0   ==>  s1 < s2
    if float(l1[2]) - float(l2[2]) != 0:
      return int(float(l1[2]) - float(l2[2]))
    return int(l1[0]) - int(l2[0])



def ParsingGff3 (filename,genelist,origin):
    with open(filename, "r") as f:
        contents = f.read()
    lines = contents.split("\n")
    #print lines
    for i in range (0, len(lines)-1):
        #ignore all other headers
        if lines[i][0] == ">": 
                continue
        cols = lines[i].split("\t")
        if int(cols[3]) > int(cols[4]):
            genelist += [[int(cols[4]),int(cols[3]),float(cols[5])]]
            origin += [[int(cols[4]),int(cols[3]),float(cols[5])]]
        else:
            genelist += [[int(cols[3]),int(cols[4]),float(cols[5])]]
    genelist.sort(mycmp)
    return

def make2dList(rows, cols, value):
    a=[]
    for row in xrange(rows):
        a += [[value]*cols]
    return a

def maxscore(instdic):
    opt = []
    select = []
    for i in xrange(len(instdic)):
      if i == 0:
        o1 = 0
      else:
        o1 = opt[i-1]
      k = -1
      for ok in xrange(i-1,-1,-1):
        if instdic[ok][1] < instdic[i][0]:
          k = ok
          break
      if k == -1 :
        o2 =instdic[i][2]
      else:
        o2 = opt[k] + instdic[i][2]
      if o1 < o2:
        if k == -1:
          ok = 0
        else:
          ok = opt[k]
        opt.append(ok + instdic[i][2])
        select.append(True)
      else:
        opt.append(opt[i-1])
        select.append(False)
    return (opt,select)

def outputk(opt,select,instdic,origin,results):
    i = len(select) - 1
    while i >= 0:
      if select[i]:
        if instdic[i] in origin:
          results += [[instdic[i][1],instdic[i][0],instdic[i][2]]]
        else:
          results += [[instdic[i][0],instdic[i][1],instdic[i][2]]]
        k = -1
        for ok in xrange(i-1,-1,-1):
          if instdic[ok][1] < instdic[i][0]:
            k = ok
            i = ok
            break
        if k == -1:
          break
      else:
        i -= 1




def combine(periodicity_gff3,naivepred_gff3,gff3_rnaoutput,glimmer_gff3):
    genelist = []
    origin = []
    results = []
    ParsingGff3(glimmer_gff3,genelist,origin)
    ParsingGff3(gff3_rnaoutput,genelist,origin)
    ParsingGff3(periodicity_gff3,genelist,origin)
    ParsingGff3(naivepred_gff3,genelist,origin)
    (opt,select) = maxscore(genelist)
    total_score = opt[-1]
    outputk(opt,select,genelist,origin,results)
    results.sort(key = itemgetter(2),reverse = True)
    return results

#print combine("a.txt","naive.gff3","rnapred.gff3","test.gff3")
    
