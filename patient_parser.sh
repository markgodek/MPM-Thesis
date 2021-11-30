INPUT=$1
SCRIPT_OUT=$2

# debug variables
#INPUT=raw_data_path/SAM_patient_link.csv
#INPUT=raw_data_path/debug_link.csv
#SCRIPT_OUT=user_dir/example_familial/scripts/complete_pipeline/debug; mkdir -p ${SCRIPT_OUT}

# clean up a previous run

rm ${SCRIPT_OUT}/PATIENT_LIST.tmp ${SCRIPT_OUT}/NORMAL_SAMPLES.tmp 2> /dev/null
rm ${SCRIPT_OUT}/ALL_SAMPLES.tmp ${SCRIPT_OUT}/TUMOR_ONLY_LINK.csv ${SCRIPT_OUT}/TUMOR_PATIENTS.tmp 2> /dev/null

SAMPLE_COUNT=$(($(wc -l < ${INPUT}) - 1))
INDEX=1

# create header for tumor only link file
head -n 1 ${INPUT} > ${SCRIPT_OUT}/TUMOR_ONLY_LINK.csv

for ((INDEX=1; INDEX<=${SAMPLE_COUNT}; INDEX++))
do
	row=(`head -n $((${INDEX} +1)) ${INPUT} | tail -n 1`)

	OIFS=$IFS;
	IFS=",";

	patient=(`head -n $((${INDEX} +1)) ${INPUT} | tail -n 1`)

	IFS=$OIFS;
	
	# trim whitespace and store in variables
	code="$(echo -e "${patient[0]}" | tr -d '[:space:]')"
	tumor_sam="$(echo -e "${patient[1]}" | tr -d '[:space:]')"
	normal_sam="$(echo -e "${patient[2]}" | tr -d '[:space:]')"

	if test ! -z "${code}"
	then	
		echo ${code} >> ${SCRIPT_OUT}/PATIENT_LIST.tmp
	fi	

	if test ! -z "${tumor_sam}"
                then
			echo ${tumor_sam} >> ${SCRIPT_OUT}/ALL_SAMPLES.tmp
			echo ${code} >> ${SCRIPT_OUT}/TUMOR_PATIENTS.tmp
                        echo ${row} >> ${SCRIPT_OUT}/TUMOR_ONLY_LINK.csv
		fi
	if test ! -z "${normal_sam}"
                then
			echo ${normal_sam} >> ${SCRIPT_OUT}/ALL_SAMPLES.tmp
			echo ${normal_sam} >> ${SCRIPT_OUT}/NORMAL_SAMPLES.tmp
		fi
done
