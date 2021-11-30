#input=all_VCFs.txt

#while IFS= read -r line
#do
#	echo "$line"
#done<${input}

DEBUG="False"
RESELECT="True"



if [[ "${DEBUG}" == "False" ]] && [[ "${RESELECT}" == "False" ]]
then
	echo "Main program"
	echo "DEBUG: ${DEBUG}"
	echo "RESELECT: ${RESELECT}"
fi

if [[ "${RESELECT}" == "True" ]]
then

	echo "RESELECTING"
	echo "RESELECT: ${RESELECT}"
fi
