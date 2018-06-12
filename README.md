# LinuxTips
A place to store lots of Tips and Tricks in Linux

```
#Super fast searching with silver searcher (ag)
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


#for a bunch of paths to .bam and .bai files, make symlinks to them and extract the lenaID for naming
for i in `cat bams.dir`; do ln -s $i `awk 'END{ var="'$i'"; split (var,a,/\//); print a[8]}' /dev/null`.bam; done
for i in `cat bais.dir`; do ln -s $i `awk 'END{ var="'$i'"; split (var,a,/\//); print a[8]}' /dev/null`.bai; done

# for a bunch of leana ids create symlinks
for i in {48868..48914}; do ln -s `find_marsoc_paths.py $i`/outs/phased_variants.vcf.gz.tbi ./$i.vcf.gz.tbi; done

#upload a file to an ftp server via curl (will be useful for SRA uploads)
nohup curl -T my_file.txt ftp://ftp.place.to.go --user username:password &

#put command in sftp bash script. cd into the correct dir of the sftp server then put stuff there (recursive)
nano put.sh
      #!/bin/bash
      export SSHPASS=PassWord
      sshpass -e sftp -P 2200 'stephen.williams@10xgenomics.com'@bioshare.bioinformatics.ucdavis.edu << !
      cd s0ryf7dj2tr4tef/
      put -r /mnt/home/stephen/yard/chili_SN2.0_BUSCO3
      bye
      !
nohup ./put.sh &o


#replace characters with sed. replce foo with bar
sed -i 's/foo/bar/g' *txt

#convert whitespace to return.
sed -i 's/\s\+/\n/g' file.txt

#replace characters with sed when there are backslashes such as in file paths
sed -i 's#HipSTR#HipSTR/DBS_TS#g' make_plots.sh

#Find vcfs and cound the number of variants in parallel
find *vcf | parallel 'echo {};cat {} | grep -v "#" | wc -l'

#modified to cound the number of reads in a file (seems very fast)
find read-RA_si-GGTAACGC_lane-001-chunk-001.fastq.gz | parallel 'echo {};zcat {} | wc -l'

#start a bunch of qsub jobs at once by making a shell script
 #file structure submit.qsub.sh
 cd 45801 ; qsub < 45801.sh ; cd ..
 cd 45802 ; qsub < 45802.sh ; cd ..
 etc.
 etc.
./submit.qsub.sh

#switch between different longranger dpe
cd longranger
git checkout master (or whatever maybe wei/LONGRANGER-1759)

#delete a queued job on hydra
cd /mnt/hydra/hydra/spool/pending
rm YOURJOBNUMBER

#find .tsv files in all directories excluding ones that end with del and copy them to another directory
find ./45* -type f -name "*.tsv" ! -path "./*.del/*" -exec cp {} /mnt/home/stephen/public_html/Heidi/MacRehm_new \;

#pretty .bam/.sam files
prettysam -r hg19-2.2.0/fasta/genome.fa phased_possorted_bam.bam

#delete all the files in a directory except an *.sh file for a list of directories
shopt -s extglob 
for i in `cat dirs.txt`; do cd $i; rm -- !($i.sh); cd ..;  done

#remove spaces from a file but keep tabs (nice if you coppied sample info off LENA)
sed 's/ \+//g' file

# use awk to "grep" for something and print the second field
generate data | awk '/something/ {print $2}' 

# Ultra fast find things such as cell barcodes with ripgrep
sambamba view -t 5 possorted_genome_bam.bam | rg -j 5 --no-line-number -F -f test_bc.txt > test.rg.sam

# A tidy way to look at a vcf
zcat output.gvcf.gz|  tail -n +285 | column -t

# Add something to bash profile automatically 

echo -e '\n#Your note\nexport PATH=/path/to/the/thing:$PATH' >> ~/.bash_profile

# Send an email to yourself at the end of a bash script. Has subject "your script is done" and content "test"

mail -s "your script is done" youremail@email.com <<< "test"

# run R script on cluster
/mnt/opt/R/R-3.4.3/bin/Rscript --vanilla /path/to/R_code.R

```


