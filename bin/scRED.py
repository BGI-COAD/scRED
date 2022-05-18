#scRED-v1.0
import os,sys,re
import argparse
import time
import configparser
import subprocess

def parse_options():
    parser = argparse.ArgumentParser(description='scRED 1.0')
    parser.add_argument('-b', '--rnabam', help='The bam file to be analyzed')
    parser.add_argument('-o', '--outdir', help='The output directory')
    parser.add_argument('-r', '--reference', help='The reference FASTA file')
    parser.add_argument('-rg', '--refGene', help='The reference gene file')
    parser.add_argument('-sr', '--simpleRepeat', help='The simple repeat file')
    parser.add_argument('-a', '--alufile', help='The ALU file')
    parser.add_argument('-ds', '--dbsnpfile', help="The snp file")
    parser.add_argument('-an', '--annovardir', help="The annovar directory")
    parser.add_argument('-dt', '--deltempfile', help="Whether delete temp files(y/n)?")
    args = parser.parse_args()
    return args
        
def get_time():
    return "["+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())+"]"

def run_subprocess(cmd,log_out,log_err):
    p = subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False)
    stdout,stderr = p.communicate()
    stdout=stdout.decode()
    stderr=stderr.decode()
    print(stdout,end='')
    if stderr!='':
        print(stderr,end='')
    log_out+=stdout
    log_err+=stderr
    return log_out, log_err

def analyze(args):
    log_out=str()
    log_err=str()

    cmd=['echo', get_time(), " Get all parameters ...",]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    
    bindir = os.path.abspath(os.path.dirname(__file__))
    rnabam=args.rnabam
    outdir=args.outdir
    ref=args.reference
    refGene=args.refGene
    simpleRepeat=args.simpleRepeat
    alu=args.alufile
    dbsnp=args.dbsnpfile
    annovar=args.annovardir
    deltf=args.deltempfile
    
    cmd=['echo', get_time(), " Calculating mismatch ratio ...",]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    cmd=[bindir+"/MismatchStat", '-i', rnabam, '-x', '1000000', '-u', '-o', outdir+"/rna_mismatch_stat.txt"]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    
    cmd=['echo', get_time(), " Variation sites detecting ..."]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    cmd=[bindir+"/MutDet",'-i',rnabam,'-r',ref,'-q','50','-v',outdir+"/rna_mismatch_stat.txt",'-u','-o',outdir+"/variation.sites.txt.gz"]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    
    cmd=['echo', get_time(), " Binomial,variation reads,editing frequency,strandbias,reads end,simple repeat and homopolymer filtering...\n",]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    cmd=['perl',bindir+"/binom.reads.fre.end.strand.mism.poly.rep.filter.pl",'--in',outdir+"/variation.sites.txt.gz",'--out',outdir+"/binom.sp.hp.filter.sites",\
    '--simpleRepeat',simpleRepeat,'--reference',ref,'--endfre','0.9','--endp','0.05','--strandfre','0.9','--strandp','0.005','--homopolymer','5']
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    
    cmd=['echo', get_time(), " Removing the sites which in snp database..."]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    cmd=['python',bindir+"/filt.dbsnp.py",dbsnp, outdir+"/binom.sp.hp.filter.sites",  outdir+"/binom.sp.hp.filter.sites.dbsnp.txt"]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    
    cmd=['echo', get_time(), " RNA editing sites sorting..."]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    os.system('sort -k1.4n -k2n %s/binom.sp.hp.filter.sites.dbsnp.txt > %s/binom.sp.hp.filter.sites.dbsnp.sort.txt'\
    % (outdir, outdir))
    os.system('mv %s/binom.sp.hp.filter.sites.dbsnp.sort.txt %s/binom.sp.hp.filter.sites.dbsnp.txt'\
    % (outdir,outdir))
    
    cmd=['echo', get_time(), " Candidate RNA editing sites annotation ..."]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    cmd=['perl',bindir+"/annovar.pl",outdir+"/binom.sp.hp.filter.sites.dbsnp.txt",annovar,annovar+"/mousedb_mm10/mm10_refGene.txt",alu,outdir+"/RNA_editing.sites.annotation.txt"]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    
    cmd=['echo', get_time(), " Splicing region filtering ..."]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    cmd=['perl',bindir+"/alu.splicing.filter.pl",outdir+"/RNA_editing.sites.annotation.txt",outdir+"/RNA_editing.sites.annotation.splicing.filter.txt"]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    
    cmd=['echo', get_time(), " Basic information counting ..."]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    cmd=['perl',bindir+"/stat.pl",outdir+"/RNA_editing.sites.annotation.splicing.filter.txt",outdir+"/TypeDist.stat.txt",outdir+"/GenomeDist.stat.txt"]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)

    if deltf== 'y':
        cmd=['echo', get_time(), " Delete temp files ..."]
        log_out,log_err=run_subprocess(cmd,log_out,log_err)
        for root,dirs,files in os.walk(outdir):
            for f in files:
                if re.search("binom",f):
                    abs_f=os.path.join(root,f)
                    os.remove(abs_f)
    
    cmd=['echo', get_time(), " Completed!"]
    log_out,log_err=run_subprocess(cmd,log_out,log_err)
    
    log_out_f=open(outdir+"/log.out","w")
    log_out_f.write(log_out)
    log_err_f=open(outdir+"/log.err","w")
    log_err_f.write(log_err)

if __name__ == '__main__':
    args = parse_options()
    analyze(args)