import os,sys,re
path1=sys.argv[1]
path2=sys.argv[2]
path3=sys.argv[3]
outf=sys.argv[4]
fw=open(outf,'w')
pathLs = [path1,path2,path3]
fDt={}
for i in pathLs:
    for root,dirs,files in os.walk(i):
        for f in files:
            if re.search('accepted',f):
                absP=os.path.join(root,f)
                k=absP.split('/')[-2]
                fDt.setdefault(k,[]).append(absP) 
                #fLs.append(os.path.join(root,f))
print(len(fDt.keys()))
for k in sorted(fDt.keys()):
    fw.write(k+','+fDt[k][0]+'\n')
