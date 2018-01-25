## linking genomes to my folder

cd ~/data/data_genomes
mkdir hg19 mm9 mm10 hg18
for genome in hg19 mm9 mm10 hg18; do
	ln -s /k/genomes/$genome/chrom.sizes /home/daria/data/data_genomes/$genome/.
done

## copying genome sizes and removing all unnecessary chromosomes

for genome in hg19 mm9 mm10 hg18; do
	cp /k/genomes/$genome/chrom.sizes /home/daria/data/data_genomes/$genome/chrom.sizes.good
done

## Edit chrom.sizes.good to remove junk chromosomes (random, Un, M, etc)

emacs mm9/chrom.sizes.good
emacs hg19/chrom.sizes.good
emacs hg18/chrom.sizes.good
emacs dm3/chrom.sizes.good

## linking all genomes fasta files

for genome in hg19 mm9 mm10 hg18; do
	mkdir -p /home/daria/data/data_genomes/$genome/fa/
	ln -s /k/genomes/$genome/fa/*.fa /home/daria/data/data_genomes/$genome/fa/
done

## bowtie indices for al chromosomes, but only canonical

for genome in hg19 mm9 mm10 hg18; do
	mkdir -p /home/daria/data/data_genomes/$genome/index/bowtie2_canonical
	ln -s /k/genomes/$genome/index/bowtie2_canonical/* /home/daria/data/data_genomes/$genome/index/bowtie2_canonical/
done

for genome in hg19 mm9 mm10 hg18; do
	mkdir -p /home/daria/data/data_genomes/$genome/index/bowtie_canonical
	ln -s /k/genomes/$genome/index/bowtie_canonical/* /home/daria/data/data_genomes/$genome/index/bowtie_canonical/
done

## Motifs

mkdir -p /NextGenSeqData/project-data/daria/data_genomes/motifs

# Jaspar

mkdir -p /NextGenSeqData/project-data/daria/data_genomes/motifs/jaspar
cd /NextGenSeqData/project-data/daria/data_genomes/motifs/jaspar
# downloaded May 2016: http://jaspar.genereg.net/html/DOWNLOAD/all_data/FlatFileDir/FlatFileDir.tar.gz manually; matrix file is from here: http://jaspar.genereg.net/html/DOWNLOAD/all_data/FlatFileDir/matrix_list.txt; did fle on matrix file

jaspar2meme -pfm -logodds FlatFileDir | awk '{if($1=="MOTIF"){gsub(/\./,"_",$2);gsub(/::/,"_",$3);$2=$2"_"$3;$3=""};print $0}' > jaspar_core_nonredundant_vertebrates.meme
#rm FlatFileDir.tar.gz
#rm FlatFileDir/*
#rmdir FlatFileDir/*
mkdir jaspar_core_nonredundant_vertebrates # to store motifs to read meme files for mast.

awk '{if($1=="MOTIF"){F=1;file="jaspar_core_nonredundant_vertebrates/"$2".meme";print file>file;print header>file;print $0>file}else if(F==1){print $0>file}else{if(!header){header=$0}else(header=header"\n"$0)}}' jaspar_core_nonredundant_vertebrates.meme

meme2images jaspar_core_nonredundant_vertebrates.meme logo/

# error: awk: cmd. line:1: (FILENAME=jaspar_core_nonredundant_vertebrates.meme FNR=11268) fatal: can't redirect to `jaspar_core_nonredundant_vertebrates/MA0045_1_HMG-I/Y.meme' (No such file or directory)
# It seems only one case like this: corrected name manually (HMG-I/Y --> HMG-I-Y), otherwise impossible to run the whole skript.
# despite the error: number of motifs in jaspar_core_nonredundant_vertebrates.meme: 2049; in logo folder: 4098 (png and esp files)

# mapping motifs to the genomes
screen -R mast
for genome in hg19 ; do
mkdir -p mast_p03/${genome}
for motif in `ls -1 jaspar_core_nonredundant_vertebrates/ | sed 's/\.meme//'`; do
    for chr in `cat /home/daria/data/data_genomes/${genome}/chrom.sizes.good | cut -f1 | xargs`; do
	mast -hit_list -mt 0.001 jaspar_core_nonredundant_vertebrates/${motif}.meme /home/daria/data/data_genomes/${genome}/fa/${chr}.fa | awk -vOFS='\t' -vM=$motif '($1!~/#/){if($2=="+1"){s="+"}else{s="-"}print $1,$3,$4,M,$6,s,$5}'
    done | gzip > mast_p03/${genome}/${motif}.txt.gz
done &
done

for genome in hg19; do
for motif in `ls -1 jaspar_core_nonredundant_vertebrates/ | sed 's/\.meme//'`; do
    echo $motif
    gunzip -c mast_p03/${genome}/${motif}.txt.gz | cut -f1 | uniq | head -2
done
done

################## stop for now

# Uniprobe Bulyk
# addtional programs, added to bash

mkdir -p /NextGenSeqData/project-data/daria/data_genomes/motifs/uniprobe
cd /NextGenSeqData/project-data/daria/data_genomes/motifs/uniprobe
mkdir -p SCI09
wget http://the_brain.bwh.harvard.edu/uniprobe/downloads/SCI09/104_pwm_all_separate.zip
unzip 104_pwm_all_separate.zip -d SCI09

mkdir -p Cell08
wget http://the_brain.bwh.harvard.edu/uniprobe/downloads/Cell08/homeo_pwm_all_separate.zip
unzip homeo_pwm_all_separate.zip -d Cell08

mkdir -p mmu
uniprobe2meme=~/utils/meme_4.10.2/scripts/uniprobe2meme # could not add to bashrc
meme2images=~/utils/meme_4.10.2/src/meme2images

(cat SCI09/*_pwm_primary.txt | $uniprobe2meme -logodds | awk '{if($1=="MOTIF"){$2=$3;$3=""};print $0}' | sed 's/-primary//'
cat Cell08/*_pwm.txt | $uniprobe2meme -logodds | awk '{if($1=="MOTIF"){gsub(/\./,"_",$2)};print $0}'
)> mmu.meme
awk '{if($1=="MOTIF"){F=1;file="mmu/"$2".meme";print file>file;print header>file;print $0>file}else if(F==1){print $0>file}else{if(!header){header=$0}else(header=header"\n"$0)}}' mmu.meme
$meme2images mmu.meme logo/

cd /NextGenSeqData/project-data/daria/data_genomes/motifs/uniprobe
mast=~/utils/meme_4.10.2/src/mast

#for genome in mm9 hg19 hg18 mm10; do # did not do for mm10 and hg18
#for genome in mm9; do
for genome in hg19; do
mkdir -p mast_p03/${genome}
for motif in `ls -1 mmu/ | sed 's/\.meme//'`; do
    for chr in `cat ~/data/data_genomes/$genome/chrom.sizes.good | cut -f1 | xargs`; do
	$mast -hit_list -mt 0.001 mmu/${motif}.meme ~/data/data_genomes/${genome}/fa/${chr}.fa | awk -vOFS='\t' -vM=$motif '($1!~/#/){if($2=="+1"){s="+"}else{s="-"}print $1,$3,$4,M,$6,s,$5}'
	done | gzip > mast_p03/${genome}/${motif}.txt.gz
done &
done

# Repeats
mkdir ~/data/data_genomes/annotations
cd ~/data/data_genomes/annotations

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz -O ~/temp_NOBACKUP/hg19_rmsk.txt.gz
gunzip -c ~/temp_NOBACKUP/hg19_rmsk.txt.gz | awk -vOFS='\t' '{print $6,$7,$8,$10,$11,$12,$13}' | sort -k1,1 -k2,2n | gzip > hg19_repeats.txt.gz

rm ~/temp_NOBACKUP/hg19_rmsk.txt.gz

# http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/simpleRepeat.txt.gz

for chr in `cat ~/data/data_genomes/mm9/chrom.sizes.good | cut -f1 | xargs`; do
    wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/${chr}_rmsk.txt.gz -O ~/temp_NOBACKUP/mm9_${chr}_rmsk.txt.gz
    gunzip -c ~/temp_NOBACKUP/mm9_${chr}_rmsk.txt.gz | awk -vOFS='\t' '{print $6,$7,$8,$10,$11,$12,$13}'
done | sort -k1,1 -k2,2n | gzip > mm9_repeats.txt.gz

rm ~/temp_NOBACKUP/mm9_*

# liftOver files and blacklists
cd ~/data/data_genomes/mm9
wget http://www.broadinstitute.org/~anshul/projects/mouse/blacklist/mm9-blacklist.bed.gz

cd ~/data/data_genomes/hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz

# liftOverfiles are from here: http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/
cd ~/data/data_genomes/
mkdir liftOver; cd liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz 

# lifOver black lists

cd ~/data/data_genomes/mm9
chain=~/data/data_genomes/liftOver/mm9ToMm10.over.chain.gz
zcat mm9-blacklist.bed.gz > ~/data/temp_NOBACKUP/mm9-blacklist.bed
liftOver ~/data/temp_NOBACKUP/mm9-blacklist.bed $chain ~/data/temp_NOBACKUP/mm10-blacklist.bed ~/data/temp_NOBACKUP/unmapped.bed
gzip ~/data/temp_NOBACKUP/mm10-blacklist.bed
mv ~/data/temp_NOBACKUP/mm10-blacklist.bed.gz ~/data/data_genomes/mm10
rm ~/data/temp_NOBACKUP/unmapped.bed ~/data/temp_NOBACKUP/mm9-blacklist.bed
# ~25 regions were not liftable

cd ~/data/data_genomes/hg19
mv wgEncodeDacMapabilityConsensusExcludable.bed.gz hg19-blacklist.bed.gz
chain=~/data/data_genomes/liftOver/hg19ToHg38.over.chain.gz
zcat hg19-blacklist.bed.gz > ~/data/temp_NOBACKUP/hg19-blacklist.bed
liftOver ~/data/temp_NOBACKUP/hg19-blacklist.bed $chain ~/data/temp_NOBACKUP/hg38-blacklist.bed ~/data/temp_NOBACKUP/unmapped.bed
gzip ~/data/temp_NOBACKUP/hg38-blacklist.bed
mv ~/data/temp_NOBACKUP/hg38-blacklist.bed.gz ~/data/data_genomes/hg38
rm ~/data/temp_NOBACKUP/hg19-blacklist.bed ~/data/temp_NOBACKUP/unmapped.bed
# 10 regions were not liftable


## Gene ontology
mkdir -p ~/data/data_genomes/geneOntology
cd ~/data/data_genomes/geneOntology

wget http://geneontology.org/ontology/go.obo
wget http://www.geneontology.org/gene-associations/gene_association.goa_human.gz
wget http://www.geneontology.org/gene-associations/gene_association.mgi.gz

gunzip -c gene_association.mgi.gz | awk 'BEGIN{while(getline<"go.obo"){if($1=="[Term]"){F=1}else if($1~/\[/){F=0}else if($1=="id:"&&F==1){go=$2}else if($1=="name:"&&F==1){N[go]=$2;for(i=3;i<=NF;i++){N[go]=N[go]"_"$i};gsub(/[^A-Za-z0-9_]/,"",N[go])}else if($1=="namespace:"&&F==1){NS[go]=$2}else if($1=="is_a:"){if(!P[go]){P[go]=$2}else{P[go]=P[go]";"$2}}}}($1!~"!"&&$4~"GO"){go=$4;g=$3;if(NS[go]=="biological_process"){print g"\t"go"\t"N[go];if(P[go]){list=P[go];while(list!=0){split(list,X,";");print g"\t"X[1]"\t"N[X[1]];if(P[X[1]]){list=list";"P[X[1]]};if(list~";"){sub(/GO:[0-9]*;/,"",list)}else{list=0}}}}}' | sort -k1,1 | uniq > GO_biological_process_mouse.txt

gunzip -c gene_association.goa_human.gz | awk 'BEGIN{while(getline<"go.obo"){if($1=="[Term]"){F=1}else if($1~/\[/){F=0}else if($1=="id:"&&F==1){go=$2}else if($1=="name:"&&F==1){N[go]=$2;for(i=3;i<=NF;i++){N[go]=N[go]"_"$i};gsub(/[^A-Za-z0-9_]/,"",N[go])}else if($1=="namespace:"&&F==1){NS[go]=$2}else if($1=="is_a:"){if(!P[go]){P[go]=$2}else{P[go]=P[go]";"$2}}}}($1!~"!"&&$4~"GO"){go=$4;g=$3;if(NS[go]=="biological_process"){print g"\t"go"\t"N[go];if(P[go]){list=P[go];while(list!=0){split(list,X,";");print g"\t"X[1]"\t"N[X[1]];if(P[X[1]]){list=list";"P[X[1]]};if(list~";"){sub(/GO:[0-9]*;/,"",list)}else{list=0}}}}}' | sort -k1,1 | uniq > GO_biological_process_human.txt

## Pathways # wiki pathways later!


## mappability files and blacklists from UCSCS

# Create reads covering each genome using a sliding window of N bp by step of 1
# toupper converts lowere case to upper 

screen -R mappability
cd ~/data/data_genomes
for genome in mm10 hg38; do
    for N in 50 36; do
        echo $genome $N
        for fasta in `ls -1 ~/data/data_genomes/$genome/fa/*.fa`; do
            cat $fasta
        done | awk -vN=$N '{if(/^>/){s=""}else{s=(substr(s,length(s)-(N-2)) toupper($1));for(i=1;i<=(length(s)-(N-1));i++){print substr(s,i,N)}}}' | gzip > $genome/${genome}_${N}bp_reads.gz
    done
done

cd ~/data/data_genomes
# from Anais: (mm9 36bp => 6504m23.200s). uses 5 cores!
# mm10
for genome in hg38; do
    for N in 50; do
        echo $genome $N
        gunzip -c ${genome}/${genome}_${N}bp_reads.gz | bowtie --chunkmbs 200 -p 5 -r --sam -m 1 -v 3 --best --strata ~/data/data_genomes/$genome/index/bowtie_canonical/$genome - > $genome/${genome}_${N}bp.sam
        samtools view -Sb $genome/${genome}_${N}bp.sam | samtools sort - $genome/${genome}_${N}bp
        samtools index $genome/${genome}_${N}bp.bam
        rm $genome/${genome}_${N}bp_reads.gz $genome/${genome}_${N}bp.sam
    done
done

# Percent mapped
cd ~/data/data_genomes
#genome=mm10
genome=hg38
    for N in 50; do
	echo -en $genome"\t"$N"\t"
	samtools idxstats $genome/${genome}_${N}bp.bam | awk '{m=m+$3;u=u+$4}END{print m"\t"m*100/(m+u)}'
    done >> $genome/mapping_summary.txt
cat $genome/mapping_summary.txt | aa

# Mappable regions
# Only very few positions have >2 reads
for genome in mm10; do
    for N in 36 50; do
	bamToBed -i $genome/${genome}_${N}bp.bam | mergeBed -i stdin | gzip > $genome/${genome}_${N}bp_mappable.bed.gz
	rm $genome/${genome}_${N}bp.bam $genome/${genome}_${N}bp.bam.bai
    done
done

# Coverage
genome=mm10
    for N in 36 50; do
	echo -en $genome"\t"$N"\t"
	t=$(awk '{t=t+$2}END{print t}' /home/daria/data/data_genomes/${genome}/chrom.sizes)
	gunzip -c $genome/${genome}_${N}bp_mappable.bed.gz | awk -vt=$t '{m=m+$3-$2}END{print m"\t"m*100/t}'
    done >> $genome/coverage_summary.txt
cat $genome/coverage_summary.txt | aa

####################### for hg19 mappability

screen -R mappability
cd ~/data/data_genomes
# cutting genome into 50 bp fragments
for genome in hg19; do
    for N in 50; do
        echo $genome $N
        for fasta in `ls -1 ~/data/data_genomes/$genome/fa/*.fa`; do
            cat $fasta
        done | awk -vN=$N '{if(/^>/){s=""}else{s=(substr(s,length(s)-(N-2)) toupper($1));for(i=1;i<=(length(s)-(N-1));i++){print substr(s,i,N)}}}' | gzip > $genome/${genome}_${N}bp_reads.gz
    done
done

cd ~/data/data_genomes
# from Anais: (mm9 36bp => 6504m23.200s). uses 5 cores!
for genome in hg19; do
    for N in 50; do
        echo $genome $N
        gunzip -c ${genome}/${genome}_${N}bp_reads.gz | bowtie --chunkmbs 200 -p 5 -r --sam -m 1 -v 3 --best --strata ~/data/data_genomes/$genome/index/bowtie_canonical/$genome - > $genome/${genome}_${N}bp.sam
        samtools view -Sb $genome/${genome}_${N}bp.sam | samtools sort - $genome/${genome}_${N}bp
        samtools index $genome/${genome}_${N}bp.bam
        rm $genome/${genome}_${N}bp_reads.gz $genome/${genome}_${N}bp.sam
    done
done

# Percent mapped
cd ~/data/data_genomes
genome=hg19
    for N in 50; do
	echo -en $genome"\t"$N"\t"
	samtools idxstats $genome/${genome}_${N}bp.bam | awk '{m=m+$3;u=u+$4}END{print m"\t"m*100/(m+u)}'
    done >> $genome/mapping_summary.txt
cat $genome/mapping_summary.txt | aa
#hg19  50  2656472821  84.6777

# Mappable regions
# Only very few positions have >2 reads
for genome in hg19; do
    for N in 50; do
	bamToBed -i $genome/${genome}_${N}bp.bam | mergeBed -i stdin | gzip > $genome/${genome}_${N}bp_mappable.bed.gz
	rm $genome/${genome}_${N}bp.bam $genome/${genome}_${N}bp.bam.bai
    done
done

# Coverage
genome=hg19
    for N in 50; do
	echo -en $genome"\t"$N"\t"
	t=$(awk '{t=t+$2}END{print t}' /home/daria/data/data_genomes/${genome}/chrom.sizes)
	gunzip -c $genome/${genome}_${N}bp_mappable.bed.gz | awk -vt=$t '{m=m+$3-$2}END{print m"\t"m*100/t}'
    done >> $genome/coverage_summary.txt
cat $genome/coverage_summary.txt | aa
# hg19  50  2742962275  87.4345

#########################################################genome annotations

cd ~/data/data_genomes
mkdir -p ~/data/data_genomes/annotations/raw
cd ~/data/data_genomes/annotations/raw

# hg19
# from ucsc table browser refFlat
# GRCh37.74 the lastest and the last release for hg19
wget ftp://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens/Homo_sapiens.GRCh37.74.gtf.gz

gunzip -c Homo_sapiens.GRCh37.74.gtf.gz | awk -vOFS='\t' '($2=="protein_coding"&&(($1>=1&&$1<=22)||$1=="Y"||$1=="X")){print "chr"$0}' | sort -k1,1 -k4,4n > Homo_sapiens.GRCh37.74.gtf
gunzip -c Homo_sapiens.GRCh37.74.gtf.gz | awk -vFS='\t' -vOFS='\t' -vC="/home/daria/data/data_genomes/hg19/chrom.sizes" 'BEGIN{while(getline<C){S[$1]=$2}}{if("chr"$1 in S && $2=="protein_coding"){gid=$9;sub(/gene_id "/,"",gid);sub(/";.*/,"",gid);tid=$9;sub(/.*transcript_id "/,"",tid);sub(/";.*/,"",tid);gn=$9;sub(/.*gene_name "/,"",gn);sub(/";.*/,"",gn);print "chr"$1,$4,$5,$7,$3,gid,tid,gn}}' | sort -k1,1 -k2,2n | gzip > hg19_annotation.txt.gz

# mm9
wget ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz
gunzip -c Mus_musculus.NCBIM37.67.gtf.gz | awk '($2=="protein_coding"&&$1!~"NT"){if($1=="MT"){$1="M"};print "chr"$0}' | sort -k1,1 -k4,4n > Mus_musculus.NCBIM37.67.gtf
awk -vFS='\t' -vOFS='\t' -vC="/home/daria/data/data_genomes/mm9/chrom.sizes" 'BEGIN{while(getline<C){S[$1]=$2}}{if($1 in S){gid=$9;sub(/gene_id "/,"",gid);sub(/";.*/,"",gid);tid=$9;sub(/.*transcript_id "/,"",tid);sub(/";.*/,"",tid);gn=$9;sub(/.*gene_name "/,"",gn);sub(/";.*/,"",gn);print $1,$4,$5,$7,$3,gid,tid,gn}}' Mus_musculus.NCBIM37.67.gtf | sort -k1,1 -k2,2n | gzip > mm9_annotation.txt.gz

# mm10
# GRCm38.76
wget ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/Mus_musculus.GRCm38.76.gtf.gz
gunzip -c Mus_musculus.GRCm38.76.gtf.gz | awk '($2=="protein_coding"&&$1!~"NT"){if($1=="MT"){$1="M"};print "chr"$0}' | sort -k1,1 -k4,4n > Mus_musculus.GRCm38.76.gtf
awk -vFS='\t' -vOFS='\t' -vC="/home/daria/data/data_genomes/mm10/chrom.sizes" 'BEGIN{while(getline<C){S[$1]=$2}}{if($1 in S){gid=$9;sub(/gene_id "/,"",gid);sub(/";.*/,"",gid);tid=$9;sub(/.*transcript_id "/,"",tid);sub(/";.*/,"",tid);gn=$9;sub(/.*gene_name "/,"",gn);sub(/";.*/,"",gn);print $1,$4,$5,$7,$3,gid,tid,gn}}' Mus_musculus.GRCm38.76.gtf | sort -k1,1 -k2,2n | gzip > mm10_annotation.txt.gz

# only for mm10! 22013 entries: ENST (col7 - still reads as one column) ids were extracted with gene_id and a spce in between. Remove them
zcat mm10_annotation.txt.gz | awk -vFS='\t' -vOFS='\t' '{if($7~/gene_id/){split($7,new,"\"");newcol=new[2]}else{newcol=$7}; print $1,$2,$3,$4,$5,$6,newcol,$8}' | sort -k1,1 -k2,2n | gzip > /tmp/mm10.txt.gz

mv /tmp/mm10.txt.gz mm10_annotation.txt.gz

rm *.gtf.gz *.gtf

# Create mRNA and TSS positions
for genome in hg19 mm10 mm9; do
    gunzip -c ${genome}_annotation.txt.gz | awk -vOFS='\t' '($5="exon"){ch[$8"\t"$7]=$1;sd[$8"\t"$7]=$4;if(!min[$8"\t"$7]||$2<min[$8"\t"$7]){min[$8"\t"$7]=$2};if(!max[$8"\t"$7]||$3>max[$8"\t"$7]){max[$8"\t"$7]=$3}}END{for(t in ch){print ch[t],min[t],max[t],"mRNA","0",sd[t],t}}' | sort -k1,1 -k2,2n > ${genome}_mrna.txt
    awk -vOFS='\t' '{if($6=="+"){print $1,$2,$2,"TSS",$5,$6,$7,$8}else if($6=="-"){print $1,$3,$3,"TSS",$5,$6,$7,$8}}' ${genome}_mrna.txt > ${genome}_tss.txt
done

# Create genomic intervals in a BED format (chr/start/end/name/score/strand + gene/transcript_id)
TF1=$(mktemp)
TF2=$(mktemp)
# only human genome for now
for genome in hg19 mm10 mm9; do
    CHROM=/home/daria/data/data_genomes/$genome/chrom.sizes
    (# CDS
    gunzip -c ${genome}_annotation.txt.gz | awk -vOFS='\t' '($5=="CDS"){print $1,$2,$3,"CDS"}'
    # 5'UTR (1st exon - 1st CDS) & 3'UTR (last exon - last CDS)
    subtractBed -a <(gunzip -c ${genome}_annotation.txt.gz | awk -vOFS='\t' '($5=="exon")') -b <(gunzip -c ${genome}_annotation.txt.gz | awk -vOFS='\t' '($5=="CDS")') | sort -k1,1 -k2,2n  | closestBed -D "a" -a stdin -b <(gunzip -c ${genome}_annotation.txt.gz | awk -vOFS='\t' '($5=="CDS")') | awk -vOFS='\t' '{if($4=="+"){if($NF>0){print $1,$2,$3,"5UTR"}else{print $1,$2,$3,"3UTR"}}else if($4=="-"){if($NF>0){print $1,$2,$3,"3UTR"}else{print $1,$2,$3,"5UTR"}}}' | uniq
    # 5'UTR (1st exon - 1st CDS)
    #gunzip -c ${genome}_annotation.txt.gz | awk -vOFS='\t' '{ch[$8"\t"$7]=$1;sd[$8"\t"$7]=$4;if($4=="+"){if($5=="exon"&&(!min[$8"\t"$7]||$2<min[$8"\t"$7])){min[$8"\t"$7]=$2};if($5=="CDS"&&(!max[$8"\t"$7]||($2-1)<max[$8"\t"$7])){max[$8"\t"$7]=$2-1}}else if($4=="-"){if($5=="exon"&&(!max[$8"\t"$7]||$3>max[$8"\t"$7])){max[$8"\t"$7]=$3};if($5=="CDS"&&(!min[$8"\t"$7]||($3+1)>min[$8"\t"$7])){min[$8"\t"$7]=$3+1}}}END{for(t in ch){if(min[t]<max[t]){print ch[t],min[t],max[t],"5UTR"}}}' | sort -k1,1 -k2,2n | mergeBed -i stdin | awk -vOFS='\t' '{print $0,"5UTR"}'
    # 3'UTR (last exon - last CDS)
    #gunzip -c ${genome}_annotation.txt.gz | awk -vOFS='\t' '{ch[$8"\t"$7]=$1;sd[$8"\t"$7]=$4;if($4=="-"){if($5=="exon"&&(!min[$8"\t"$7]||$2<min[$8"\t"$7])){min[$8"\t"$7]=$2};if($5=="CDS"&&(!max[$8"\t"$7]||($2-1)<max[$8"\t"$7])){max[$8"\t"$7]=$2-1}}else if($4=="+"){if($5=="exon"&&(!max[$8"\t"$7]||$3>max[$8"\t"$7])){max[$8"\t"$7]=$3};if($5=="CDS"&&(!min[$8"\t"$7]||($3+1)>min[$8"\t"$7])){min[$8"\t"$7]=$3+1}}}END{for(t in ch){if(min[t]<max[t]){print ch[t],min[t],max[t],"3UTR"}}}' | sort -k1,1 -k2,2n | mergeBed -i stdin | awk -vOFS='\t' '{print $0,"3UTR"}'
    ) | sort -k1,1 -k2,2n > $TF1
    (cat $TF1
    # Intronic regions (mRNA - CDS - 5'UTR - 3'UTR) (!!! change to 1-based)
    subtractBed -a ${genome}_mrna.txt -b $TF1 | awk -vOFS='\t' '{if($3-$2>1){print $1,$2+1,$3-1}}' | sort -k1,1 -k2,2n | mergeBed -i stdin | awk -vOFS='\t' '{print $0,"INTRON"}'
    # Promoter regions (2Kb) (stop at the next/previous mRNA)
    awk -vL=2000 -vC=$CHROM -vOFS='\t' 'BEGIN{while(getline<C){S[$1]=$2}}{if(os="-"){min=op+1;if($1==och&&$2<=(op+2000)){max=$2-1}else{max=op+L};if(max>S[och]){max=S[och]};if(min<max){print och,min,max}};if($6=="+"){max=$2-1;if($1==och&&($2-L)<=op){min=op+1}else{min=$2-L};if(min<1){min=1};if(min<max){print $1,min,max}};och=$1;op=$3;os=$6;og=$7"\t"$8}' ${genome}_mrna.txt | sort -k1,1 -k2,2n | awk -vOFS='\t' '{print $0,"P2000"}'
    ) | sort -k1,1 -k2,2n > $TF2
    # Intergenic regions
    (cat $TF2
    awk '{print $1"\t0\t"$2}' $CHROM | subtractBed -a stdin -b $TF2 | awk -vOFS='\t' '{if($3-$2>1){print $1,$2+1,$3-1,"INTER"}}'
    ) | sort -k1,1 -k2,2n | awk -vOFS='\t' '{print $4,$0}' | prec-overlap CDS 5UTR 3UTR INTRON P2000 INTER | awk -vOFS='\t' '{print $2,$3,$4,$1}' | sort -k1,1 -k2,2n > ${genome}_genomic_features.txt
done
rm $TF1 $TF2

################ cutting genome into smaller fragments to get better genomic distribution for random regions

#cd ~/data/data_genomes/annotations/raw
# median of regions: smallest CDS 122 bp, 3UTR: 363, etc. Median all: 240 bp
#awk '($4=="INTER"){print $3-$2}' hg19_genomic_features.txt | median.R -i - 

# cutting hg19 into smaller pieces
cd ~/data/data_genomes/annotations/raw
bedtools makewindows -b hg19_genomic_features.txt -w 200 -i src > /tmp/random_regions_window200.txt&
gzip /tmp/random_regions_window200.txt 
mv /tmp/random_regions_window200.txt.gz .

################ transcriptome index files for a specific annotation for tophat2
GTF=/k/genomes/hg19/ann/2015_01_19/derived/refGene_hg19_all_excl_multi.gtf
GENOME=~/data/data_genomes/hg19/index/bowtie2_canonical/hg19

# takes ages, but need to build only once
cd ~/data/data_genomes/hg19
tophat2 -G $GTF --transcriptome-index=transcriptome_data/known_all_excl_multi $GENOME

GTF=/k/genomes/hg19/ann/2015_01_19/derived/refGene_hg19_longest_excl_multi.gtf
GENOME=~/data/data_genomes/hg19/index/bowtie2_canonical/hg19
cd ~/data/data_genomes/hg19
tophat2 -G $GTF --transcriptome-index=transcriptome_data/known_long_excl_multi $GENOME

# on GTF data with no exclusions
GTF=/k/genomes/hg19/ann/cur/derived/refGene_hg19.gtf
GENOME=~/data/data_genomes/hg19/index/bowtie2_canonical/hg19

cd ~/data/data_genomes/hg19
tophat2 -G $GTF --transcriptome-index=transcriptome_data/known_refGene_hg19 $GENOME

##################### fa file for tophat: build from bowtie indexes and move

# using tophat to build fa file, then transfer to bowtie2_canonical folder

screen -R expression
cd ~/data/temp_NOBACKUP

head -n 4000000 ~/data/published_data/expression/Rathert_2015/fastq/SRR1688372_1.trimmed.fq > test1
head -n 4000000 ~/data/published_data/expression/Rathert_2015/fastq/SRR1688372_1.trimmed.fq > test2

INDEX=/home/daria/data/data_genomes/hg19/index/bowtie2_canonical/hg19
TRANSCRIPTOME=~/data/data_genomes/hg19/transcriptome_data/known_refGene_hg19
LIBRARY=fr-firststrand
mkdir -p tophat; outFolder=tophat

tophat2 -p 5 --b2-very-fast -o $outFolder --transcriptome-index=$TRANSCRIPTOME --no-coverage-search --library-type=$LIBRARY $INDEX test1 test2

# moving hg19.fa to /home/daria/data/data_genomes/hg19/index/bowtie2_canonical/

mv /home/daria/data/data_genomes/hg19/tophat_out/tmp/hg19.fa /home/daria/data/data_genomes/hg19/index/bowtie2_canonical/hg19.fa
samtools faidx /home/daria/data/data_genomes/hg19/index/bowtie2_canonical/hg19.fa

######################
# fa file hg19

twoBitToFa /k/genomes/hg19/2bit/hg19.2bit /home/daria/data/data_genomes/hg19/hg19.fa
samtools faidx /home/daria/data/data_genomes/hg19/hg19.fa

######################

########################################### CpG islands: download and calculations
mkdir -p ~/data/data_genomes/annotations/raw
cd ~/data/data_genomes/annotations

#table cpgIslandExt
#"Describes the CpG Islands (includes observed/expected ratio)"
 #  (
  # string chrom;        "Chromosome or FPC contig"
   #uint   chromStart;   "Start position in chromosome"
   #uint   chromEnd;     "End position in chromosome"
   #string name;         "CpG Island"
   #uint   length;       "Island Length"
   #uint   cpgNum;       "Number of CpGs in island"
   #uint   gcNum;        "Number of C and G in island"
   #uint   perCpg;       "Percentage of island that is CpG"
   #uint   perGc;        "Percentage of island that is C or G"
   #float  obsExp;	"Ratio of observed(cpgNum) to expected(numC*numG/length) 
		#	 CpG in island"
   #)

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/cpgIslandExt.txt.gz -O raw/mm9_cpgIslandExt.txt.gz
gunzip -c raw/mm9_cpgIslandExt.txt.gz | awk -vOFS='\t' '($2!~"_"){print $2,$3,$4}' | sort -k1,1 -k2,2n > mm9_cgi.txt

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz -O raw/hg19_cpgIslandExt.txt.gz
gunzip -c raw/hg19_cpgIslandExt.txt.gz | awk -vOFS='\t' '($2!~"_"){print $2,$3,$4}' | sort -k1,1 -k2,2n > hg19_cgi.txt

cd ~/data/data_genomes/annotations

# vmatchPattern searches on both strands so CpGs are found twice and looking for Cs also reports Gs
R #(BSgenome.Mmusculus.UCSC.mm9 is not loaded in R-3?--> need to load genomes first see below)
source("http://bioconductor.org/biocLite.R")
# check which genomes are loaded
biocLite('BSgenome.Hsapiens.UCSC.hg19') # once it is there no need to install
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm9) #this genome was already there --> do the same for mm10
x=vmatchPattern("CG",Mmusculus)
y=vmatchPattern("C",Mmusculus)
for(chr in c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY")){
write.table(as.data.frame(x[seqnames(x)==chr&strand(x)=="+"])[,1:3],"mm9_cpg.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)
write.table(as.data.frame(y[seqnames(y)==chr&strand(y)=="+"])[,1:3],"mm9_c.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)
write.table(as.data.frame(y[seqnames(y)==chr&strand(y)=="-"])[,1:3],"mm9_g.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)
}
#######################
# mm10 --> not done yet 
biocLite('BSgenome.Mmusculus.UCSC.mm10')
library(BSgenome.Mmusculus.UCSC.mm10) #this genome was already there --> do the same for mm10
x=vmatchPattern("CG",Mmusculus)
y=vmatchPattern("C",Mmusculus)
for(chr in c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY")){
write.table(as.data.frame(x[seqnames(x)==chr&strand(x)=="+"])[,1:3],"mm10_cpg.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)
write.table(as.data.frame(y[seqnames(y)==chr&strand(y)=="+"])[,1:3],"mm10_c.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)
write.table(as.data.frame(y[seqnames(y)==chr&strand(y)=="-"])[,1:3],"mm10_g.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)
}
######################
library(BSgenome.Hsapiens.UCSC.hg19)
x=vmatchPattern("CG",Hsapiens)
y=vmatchPattern("C",Hsapiens)
for(chr in c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrM","chrX","chrY")){
write.table(as.data.frame(x[seqnames(x)==chr&strand(x)=="+"])[,1:3],"hg19_cpg.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)
write.table(as.data.frame(y[seqnames(y)==chr&strand(y)=="+"])[,1:3],"hg19_c.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)
write.table(as.data.frame(y[seqnames(y)==chr&strand(y)=="-"])[,1:3],"hg19_g.txt",quote=F,row.names=F,col.names=F,sep="\t",append=T)
}
#?#
q()
n
# sort -m : merge already sorted files
#for genome in mm10; do
for genome in mm9 hg19; do
    sort -m -k1,1 -k2,2n <(awk -vOFS='\t' '{print $1,$2,$2+1,"C"}' ${genome}_c.txt) <(awk -vOFS='\t' '{print $1,$2,$2+1,"G"}' ${genome}_g.txt) | intersectBed -sorted -wao -a stdin -b ${genome}_cpg.txt | awk -vOFS='\t' '{print $1,$2,$3,$4,$NF}' | gzip > ${genome}_gc_cpg.txt.gz
    rm ${genome}_c.txt ${genome}_g.txt ${genome}_cpg.txt
done

# CpG density = nb_bp * nb_CpG / (nb_C * nb_G)
# high-, medium- and low- CpG density regions (>=4 CpG/100bp => high ; <2 CpG/100bp => low)
# min=0; 33%=3; 66%=4; max=107

# included the sort function
screen -R overlap
for genome in mm9 hg19 mm10; do   
    intersectBed -sorted -c -a <(sort -k1,1 -k2,2n /home/daria/data/data_genomes/$genome/chrom.sizes | awk -vOFS='\t' '{for(i=1;i<=$2-500;i+=100){print $1,i,i+500}}') -b <(gunzip -c ${genome}_gc_cpg.txt.gz | awk '($5==1)' | sort -k1,1 -k2,2n) | awk -vOFS='\t' '{$4=$4*100/500;print $0}' | awk -vOFS='\t' 'BEGIN{x="16:8:4:2:0";split(x,X,":")}{for(i=1;i<=length(X);i++){if($4>=X[i]){print "BIN"X[i],$0;break}}}' | prec-overlap BIN16 BIN8 BIN4 BIN2 BIN0 | awk -vOFS='\t' '{print $2,$3,$4+1,$1}' | gzip > ${genome}_cpg_regions.txt.gz
done

genome=mm10
for bin in BIN16 BIN8 BIN4 BIN2 BIN0; do
    echo -e $bin"\t"$(gunzip -c ~/data/data_genomes/annotations/${genome}_cpg_regions.txt.gz | awk -vB=$bin '($4==B)' | wc -l)"\t"$(gunzip -c ${genome}_cpg_regions.txt.gz | awk -vB=$bin '($4==B){print $3-$2}' | median.R -i -)"\t"$(gunzip -c ${genome}_cpg_regions.txt.gz | awk -vB=$bin '($4==B){print $3-$2}' | sort -nr | head -1)
done

# bin; N of such bins; median length of bins; max length
# mm10
#BIN16   809     601     1701
#BIN8    12225   801     4501
#BIN4    59273   601     13101
#BIN2    376639  701     8701
#BIN0    329842  2299    3110000

#mm9
#BIN16   781     601     1801
#BIN8    12151   801     4401
#BIN4    58150   601     13101
#BIN2    370387  701     8801
#BIN0    324647  2299    58682401

#hg19
#BIN16   1237    601     3001
#BIN8    21773   901     28201
#BIN4    119037  601     26101
#BIN2    578807  701     71801
#BIN0    481183  1799    30000399

##################

# HOCOMOCO v11

mkdir ~/data/data_genomes/motifs/HOCOMOCOv11; cd ~/data/data_genomes/motifs/HOCOMOCOv11

# meme format includes only letter probabilities but not logg odds
#Prior to release 4.7.0 of the MEME Suite MAST required a log-odds matrix section to be specified, however the current version of the MEME Suite is capable of translating the letter-probability matrix into the log-odds matrix and vice-versa.
# mast -version --> 4.11.0

wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme


mkdir -p HOCOMOCOv11_core # to store motifs to read meme files for mast.

awk '{if($1=="MOTIF"){F=1;file="HOCOMOCOv11_core/"$2".meme";print file>file;print header>file;print $0>file}else if(F==1){print $0>file}else{if(!header){header=$0}else(header=header"\n"$0)}}' HOCOMOCOv11_core_HUMAN_mono_meme_format.meme 

mkdir -p logo; mkdir -p mast_p03
meme2images HOCOMOCOv11_core_HUMAN_mono_meme_format.meme logo/

screen -R mast

for genome in hg19 ; do
mkdir -p mast_p03/${genome}
for motif in `ls -1 HOCOMOCOv11_core/ | sed 's/\.meme//'`; do
    echo -en $motif
    for chr in `cat /home/daria/data/data_genomes/${genome}/chrom.sizes.good | cut -f1 | xargs`; do
    mast -hit_list -mt 0.001 HOCOMOCOv11_core/${motif}.meme /home/daria/data/data_genomes/${genome}/fa/${chr}.fa | awk -vOFS='\t' -vM=$motif '($1!~/#/){if($2=="+1"){s="+"}else{s="-"}print $1,$3,$4,M,$6,s,$5}'
    done | gzip > mast_p03/${genome}/${motif}.txt.gz
done &
done


# HOCOMOCO for mouse

mkdir -p HOCOMOCOv11_core_mouse # to store motifs to read meme files for mast.

wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/core/MOUSE/mono/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme

awk '{if($1=="MOTIF"){F=1;file="HOCOMOCOv11_core_mouse/"$2".meme";print file>file;print header>file;print $0>file}else if(F==1){print $0>file}else{if(!header){header=$0}else(header=header"\n"$0)}}' HOCOMOCOv11_core_MOUSE_mono_meme_format.meme

mkdir -p logo; mkdir -p mast_p03
meme2images HOCOMOCOv11_core_MOUSE_mono_meme_format.meme logo/

screen -R mast

for genome in mm9 ; do
mkdir -p mast_p03/${genome}
for motif in `ls -1 HOCOMOCOv11_core_mouse/ | sed 's/\.meme//'`; do
    echo -en $motif
    for chr in `cat /home/daria/data/data_genomes/${genome}/chrom.sizes.good | cut -f1 | xargs`; do
    mast -hit_list -mt 0.001 HOCOMOCOv11_core_mouse/${motif}.meme /home/daria/data/data_genomes/${genome}/fa/${chr}.fa | awk -vOFS='\t' -vM=$motif '($1!~/#/){if($2=="+1"){s="+"}else{s="-"}print $1,$3,$4,M,$6,s,$5}'
    done | gzip > mast_p03/${genome}/${motif}.txt.gz
done &
done





