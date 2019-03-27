#!/bin/bash
FMETA=$1
CT=$(wc -l $FMETA | awk '{print $1}')

mkdir stitched_reads aligned_reads &
wait

x=2
while [ $x -le $CT ]
do

	string="sed -n ${x}p $FMETA"

	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1, $2, $3, $4, $5, $6}')
	set -- $var
		c1=$1
		c2=$2
		c3=$3
		
		if [ "${c2}" == "Kevin" ]; then
			HDR="\n@SQ\tSN:KH\tLN:200\n@RG\tID:Kevin\tPL:ILLUMINA\tSM:Kevin"
			REFERENCE=KH_corrected_15apr2016.fa
		fi

		if [ "${c2}" == "Olivier" ]; then
			HDR="\n@SQ\tSN:OH\tLN:232\n@RG\tID:Olivier\tPL:ILLUMINA\tSM:Olivier"
			REFERENCE=OH_corrected_13jan2015.fa
		fi
		
		if [ "${c2}" == "Olivier_HbG" ]; then
			HDR="\n@SQ\tSN:OH_HbG\tLN:213\n@RG\tID:Olivier\tPL:ILLUMINA\tSM:Olivier"
			REFERENCE=OH_HbG_16sep2016.fa
		fi
		
		# Chris MegaTal
		if [ "${c2}" == "Chris_CCR5" ]; then
			HDR="\n@SQ\tSN:CP_CCR5\tLN:220\n@RG\tID:Chris\tPL:ILLUMINA\tSM:Chris"
			REFERENCE=CP_CCR5_11jan2017.fa
		fi

		# Chris ZFN
		if [ "${c2}" == "Chris_CCR5-2" ]; then
			HDR="\n@SQ\tSN:CP_CCR5-2\tLN:166\n@RG\tID:Chris\tPL:ILLUMINA\tSM:Chris"
			REFERENCE=CP_CCR5-2_14feb2017.fa
		fi

		### Copy Reference
		cp ~/Resources/$REFERENCE ./

		###STEP ONE: Run PEAR to stitch reads together
		~/pear/bin/pear -f ./fastq_files/${c1}_L001_R1_001.fastq \
			-r ./fastq_files/${c1}_L001_R2_001.fastq \
			-o ./stitched_reads/${c1}

		###STEP TWO: Filter Stitched Reads Based On Average Quality Score
		~/Scripts/Python/quality_filter_fastq.py ./stitched_reads/${c1}.assembled.fastq ./stitched_reads/qc-filtered_${c1}.assembled.fastq

		###STEP THREE: Align reads to reference
		needle -asequence $REFERENCE \
			-bsequence ./stitched_reads/qc-filtered_${c1}.assembled.fastq \
			-gapopen 10.0 \
			-gapextend 0.5 \
			-aformat3 sam \
			-outfile ./aligned_reads/${c1}.sam

		###STEP FIVE: Convert SAM to BAM, and add custom header
		FIXED_ARGS='VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING'
		mv ./aligned_reads/${c1}.sam ./aligned_reads/TMP
		sed "s/VN:1.3/VN1.3${HDR}/" ./aligned_reads/TMP > ./aligned_reads/${c1}.sam
		rm -f ./aligned_reads/TMP
		IP="INPUT=./aligned_reads/${c1}.sam"
		OP="OUTPUT=./aligned_reads/${c1}.bam"
		LG="./aligned_reads/${c1}_SamFormatConverter.log"
		picard-tools SamFormatConverter $FIXED_ARGS $IP $OP >& $LG

		###STEP SIX: Run find_indels.R
		

	echo ${c1}
	x=$(( $x + 1 ))

done
