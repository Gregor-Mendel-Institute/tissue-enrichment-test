for j in 4MM perfect u_perfect
do
	echo $j
	for i in $(ls alignment_summary_$j)
	do
		file=alignment_summary_$j/$i
		echo $file $(sed.exe '6q;d' $file | sed.exe 's/Number of input reads |//') \
		$(sed.exe '9q;d' $file | sed.exe 's/Uniquely mapped reads number |//') \
		$(sed.exe '26q;d' $file | sed.exe 's/Number of reads mapped to too many loci |//')
	done
done