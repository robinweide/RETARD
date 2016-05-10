#ranges in /rrna-w500.bed are 500, 370 and 156 (ignore the remaining 34bp ranges)
for l in 500 370 156
do
	bedtools random -l $l -n 1000000 -g /net/NGSanalysis/ref/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai >> random.bed
done

bedtools shuffle -noOverlapping  -g /net/NGSanalysis/ref/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai -i random.bed | bedtools sort > random-uniq.bed
bedtools nuc -fi /net/NGSanalysis/ref/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -bed random-uniq.bed > random-uniq-gc.bed


echo <<'EOF'
read.delim("random-uniq-gc.bed", header=T)->d
d$bin <- cut(d$X8_pct_gc, breaks=seq(0.5, 0.9, 0.01))
#create 3 random bed files
for(i in 1:3) {
	sr <- tapply(1:nrow(d), d$bin, function(r) if (length(r) > 1000) sample(r, 1000) else r)
  write.table(Reduce(rbind, lapply(sr, function(r) d[r,c(1,2,3)])), file=paste0("gcbins-",i,".bed"), sep="\t", row.names=F, col.names=F, quote=F)
}

#downsample bam to 350M(1), 100M(2), 50M(4) 20M(8) 10M(10) 5M(10)
#write the downsample numbers for both bam files:

n <- c(350,rep(100,2), rep(50,4), rep(20,8), rep(10,10), rep(5,10))
m <- c(0,0:1, 0:3, 0:7,0:9, 0:9)
s <- paste0(n,"M-",m,".bam")

write.table(cbind("DRR01683x.bam", paste0("DRR01683x",s), n*1e6 / 375550749,sample(2^31,length(n))), file="seed-DRR01683x.txt", sep="\n", quote=F, row.names=F, col.names=F)
write.table(cbind("DRR01685x.bam", paste0("DRR01685x",s), n*1e6 / 387978997,sample(2^31,length(n))), file="seed-DRR01685x.txt", sep="\n", quote=F, row.names=F, col.names=F)
EOF

#determine gc content in random reference regions
parallel '/net/NGSanalysis/apps/bedtools/bedtools2-2.22.0/bin/bedtools nuc -fi /net/NGSanalysis/ref/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -bed {} > {.}-gc.txt' ::: gcbins-1.bed gcbins-2.bed gcbins-3.bed

parallel -a seed-DRR01683x.txt -N4  -j8 '/opt/java/jdk1.7.0_79/bin/java -jar /net/NGSanalysis/apps/picard/picard-tools-1.131/picard.jar DownsampleSam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE INPUT={1} OUTPUT={2} P={3} R={4}'
parallel -a seed-DRR01685x.txt -N4  -j8 '/opt/java/jdk1.7.0_79/bin/java -jar /net/NGSanalysis/apps/picard/picard-tools-1.131/picard.jar DownsampleSam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE INPUT={1} OUTPUT={2} P={3} R={4}'

#subset previous rRNA alignment based on readnames in downsampled data
ls -1 DRR01683x*M-?.bam |parallel -j6 'perl subsetbyname.pl ../DRR01683x.bam {} |samtools view -Sbu - | samtools sort -m 1000000000 - {.}-rRNA'
ls -1 DRR01685x*M-?.bam |parallel -j12 'perl subsetbyname.pl ../DRR01685x.bam {} |samtools view -Sbu - | samtools sort -m 1000000000 - {.}-rRNA'

#quantify random regions in downsampled bams
parallel -j8 '/net/NGSanalysis/apps/bedtools/bedtools2-2.22.0/bin/bedtools coverage -d -abam {2} -b {1} |gzip -c > {2.}-{1.}.txt.gz' ::: gcbins-1.bed gcbins-2.bed gcbins-3.bed ::: `ls -1 DRR0168[35]x*M-?.bam |xargs`

#determine gc content in rRNA regions
/net/NGSanalysis/apps/bedtools/bedtools2-2.22.0/bin/bedtools nuc -fi ../../rDNA-Hs.fasta -bed rrna-w500.bed > rrna-w500-gc.bed
#quantify rRNA regions in downsamples rRNA alignment
ls -1 *-?-rRNA.bam |parallel -j12 '/net/NGSanalysis/apps/bedtools/bedtools2-2.22.0/bin/bedtools coverage -d -abam {} -b rrna-w500.bed |gzip -c > {.}-rrna-w500.txt.gz'
