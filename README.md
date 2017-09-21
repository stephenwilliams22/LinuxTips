# LinuxTips
A place to store lots of Tips and Tricks in Linux

```#Super fast searching with silver searcher (ag)
for i in `cat list.txt`; do ag -w $i results.txt >> out.txt; done

#grep with color
grep --color=always 'string' file

#grep lines that lack a string especially good for '#' and pipe to other comands
grep -v

#much tidy'r less
less -S

#time a command and write it to a file (ls example)
{ time ls  ; } 2> time.txt

#kill all background processes
killall -u stephen.williams

#find all the files in a dir and sub dir and copy them to a new dir
find grape_DS_0.5_HFCT2BCXY -name "short*" > shorts
cp `cat shorts` BUSCO_PHASE_SUM/

#read in dir of .bam dirs. use samtools to get the depth at each base. only get chr17, ash the BRCA1 locus, and print to a file with the LENA ID (taken from the dir).

for i in `cat dirlist`; do samtools depth $i | awk -v f="$i" 'BEGIN{split(f,fs,"/")} $1=="chr17" && $2 >=41196312 && $2 <=41277500 {print > fs[8]}' ; done 
 
#download without verbose wget. the nohup output can get pretty big. 
nohup wget --no-verbose  https://….. &


#delete everything before ‘>SL3.0ch01’ and write to file. Want to delete chr00 in the original file
sed -ne '/>SL3.0ch01/,$ p' Ref_S_lycopersicum_chromosomes.3.00.fa > Ref_S_lycopersicum_chromosomes.3.00.1-12.fa

#remove contigs less than a specific length (called removesmalls.pl) in Apps folder
#used like $ perl removesmalls.pl 200 contigs.fasta > contigs-l200.fasta

## removesmalls.pl
#!/usr/bin/perl
use strict;
use warnings;

my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";
}


#find all .fasta files and gzip them. Could modify these to gunzip
find . -type f -name '*.fasta' -exec gzip "{}" \;
find grape_HFCT2BCXY -type f -name '*.fasta' -exec gzip "{}" \;

#get the size of symlinked file
stat -Lc %s symlink.file

#make a little .bam from a big one and keep the header
samtools view -h big_possorted_bam.bam | head -n 50 | samtools view - -bhS > little.test.bam

#grep every file in a dir
grep “>” *

#Convert whitespace to tabs
awk -v OFS="\t" '$1=$1' example.vcf > example.tabs.vcf

#pause a process and restart in the background
kill -20 PID #pause 
kill -18 PID #restart in the background

#number of reads in a fastq.gz file
parallel "echo {} && gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'" ::: *.gz


#gzip all the .fasta files in a dir
parallel gzip ::: *.fasta

#find your jobs that have completed within the last 3 hrs and look at memory used etc.
qacct 3h | grep -A 36 -B 5 "stephen.williams"

#split a .bed file by car
awk '{close(f);f=$1}{print > f".bed"}' file.bed

#copy multiple files at once to a dir. be in the original file directory	
cp -t /path/to/destination/ file1 file2 file3

#run an R script on the terminal
Rscript -e ‘your R script’

#view a .tsv file and keep the spacing inline with the header
column -t 13846.hla.bubbles.pls.scores.out | head


#extract a specific part of a fast. can put regions in too like chr2:1:2000
samtools faidx human_genome.fa
samtools faidx human_genome.fa chr1

#get lengths of a record in a fast
samtools faidx sample.fa 
cut -f1-2 sample.fa.fai


#for a bunch of paths to .vcf files, make symlinks to them and extract the lenaID for naming
for i in `cat vcf.dir`; do ln -s $i `awk 'END{ var=FILENAME; split (var,a,/\//); print a[8]}' $i`; done
```


