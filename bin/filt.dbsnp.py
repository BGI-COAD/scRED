import os,sys,re
import pandas as pd

dbsnp=sys.argv[1]
inf=sys.argv[2]
outf=sys.argv[3]

#print(outf)
df_dbsnp=pd.read_csv(dbsnp,header=None)
df=pd.read_csv(inf,sep="\t",header=None)

#df_dbsnp['site']=df_dbsnp[0]+"_"+df_dbsnp[1].astype(str)
df['site']=df[0]+"\t"+df[1].astype(str)

snp_ls=df_dbsnp[0].tolist()
site=df['site'].tolist()
vld_site = list(set(site).difference(set(snp_ls)))

df1=df[df['site'].isin(vld_site)]

df2=df1.drop(columns=['site'],axis=1)
df2.to_csv(outf,sep="\t",index=False,header=None)
