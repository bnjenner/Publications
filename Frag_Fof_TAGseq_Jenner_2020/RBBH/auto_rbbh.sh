#!/bin/bash

blast_function () {

        Q=$1
	DB=$2
        OUT=$3
	
	blastn -query $Q -db $DB \
	       -outfmt '6 qseqid sseqid pident mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qlen slen length' \
	       -out $OUT

}

rbbh () {
	
	blast_1=$1
	fasta_1=$2
	blast_2=$3
	fasta_2=$4
	out=$5

	python2 rbbh.py $blast_1 $fasta_1 $blast_2 $fasta_2 0.001 0.80 > $out	

}

chrom_key=$1

head=false

while IFS= read -r line
do
	if [[ $head == false ]]
	then
		header=$line
		head=true
	else
		
		sample_1=`echo $header | awk '{ print $1 }'`
		sample_2=`echo $header | awk '{ print $2 }'`
		sample_3=`echo $header | awk '{ print $3 }'`
		sample_4=`echo $header | awk '{ print $4 }'`
		sample_5=`echo $header | awk '{ print $5 }'`
	
		chrom_1=`echo $line | awk '{ print $1 }'`
		chrom_2=`echo $line | awk '{ print $2 }'`
		chrom_3=`echo $line | awk '{ print $3 }'`
		chrom_4=`echo $line | awk '{ print $4 }'`
		chrom_5=`echo $line | awk '{ print $5 }'`

		file_1=`ls genes_by_chromosome/${sample_1}/chromosome_${chrom_1}_genes.${sample_1}.fasta`
		file_2=`ls genes_by_chromosome/${sample_2}/chromosome_${chrom_2}_genes.${sample_2}.fasta`
		file_3=`ls genes_by_chromosome/${sample_3}/chromosome_${chrom_3}_genes.${sample_3}.fasta`
		file_4=`ls genes_by_chromosome/${sample_4}/chromosome_${chrom_4}_genes.${sample_4}.fasta`
		file_5=`ls genes_by_chromosome/${sample_5}/chromosome_${chrom_5}_genes.${sample_5}.fasta`

		[[ -d blast_results/chrom_${chrom_1} ]] || mkdir blast_results/chrom_${chrom_1}
		[[ -d rbbh_results/chrom_${chrom_1} ]] || mkdir rbbh_results/chrom_${chrom_1}		
		
		# pivot vs others
		blast_function $file_1 databases/${sample_2} blast_results/chrom_${chrom_1}/chromosome_${chrom_1}_${sample_1}_vs_${sample_2}.txt
		blast_function $file_1 databases/${sample_3} blast_results/chrom_${chrom_1}/chromosome_${chrom_1}_${sample_1}_vs_${sample_3}.txt
		blast_function $file_1 databases/${sample_4} blast_results/chrom_${chrom_1}/chromosome_${chrom_1}_${sample_1}_vs_${sample_4}.txt
		blast_function $file_1 databases/${sample_5} blast_results/chrom_${chrom_1}/chromosome_${chrom_1}_${sample_1}_vs_${sample_5}.txt


		# others vs pivot 
		blast_function $file_2 databases/${sample_1} blast_results/chrom_${chrom_1}/chromosome_${chrom_2}_${sample_2}_vs_${sample_1}.txt
                blast_function $file_3 databases/${sample_1} blast_results/chrom_${chrom_1}/chromosome_${chrom_3}_${sample_3}_vs_${sample_1}.txt
                blast_function $file_4 databases/${sample_1} blast_results/chrom_${chrom_1}/chromosome_${chrom_4}_${sample_4}_vs_${sample_1}.txt
		blast_function $file_5 databases/${sample_1} blast_results/chrom_${chrom_1}/chromosome_${chrom_5}_${sample_5}_vs_${sample_1}.txt

		rbbh blast_results/chrom_${chrom_1}/chromosome_${chrom_1}_${sample_1}_vs_${sample_2}.txt $file_1 \
		     blast_results/chrom_${chrom_1}/chromosome_${chrom_2}_${sample_2}_vs_${sample_1}.txt $file_2 \
		     rbbh_results/chrom_${chrom_1}/chrom_${chrom_1}_${sample_1}_vs_chrom_${chrom_2}_${sample_2}.txt

		rbbh blast_results/chrom_${chrom_1}/chromosome_${chrom_1}_${sample_1}_vs_${sample_3}.txt $file_1 \
                     blast_results/chrom_${chrom_1}/chromosome_${chrom_3}_${sample_3}_vs_${sample_1}.txt $file_3 \
                     rbbh_results/chrom_${chrom_1}/chrom_${chrom_1}_${sample_1}_vs_chrom_${chrom_3}_${sample_3}.txt

		rbbh blast_results/chrom_${chrom_1}/chromosome_${chrom_1}_${sample_1}_vs_${sample_4}.txt $file_1 \
                     blast_results/chrom_${chrom_1}/chromosome_${chrom_4}_${sample_4}_vs_${sample_1}.txt $file_4 \
                     rbbh_results/chrom_${chrom_1}/chrom_${chrom_1}_${sample_1}_vs_chrom_${chrom_4}_${sample_4}.txt

		rbbh blast_results/chrom_${chrom_1}/chromosome_${chrom_1}_${sample_1}_vs_${sample_5}.txt $file_1 \
                     blast_results/chrom_${chrom_1}/chromosome_${chrom_5}_${sample_5}_vs_${sample_1}.txt $file_5 \
                     rbbh_results/chrom_${chrom_1}/chrom_${chrom_1}_${sample_1}_vs_chrom_${chrom_5}_${sample_5}.txt

	fi
done < "$chrom_key"
