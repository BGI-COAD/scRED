import os,sys,re
import pandas as pd
inf=sys.argv[1]
outf=sys.argv[2]

df=pd.read_csv(inf,sep="\t")
df['ratio']=df['Alteration_support_reads']/df['Depth']

df1=df[(df['ratio']<1) & (df['Genotype']!=".")]
df1=df1.drop(columns=['ratio'],axis=1)
df1.to_csv(outf,index=False,sep="\t")
