#
# Randomly "deletes" (changes to 'x') 5% of the SNP values in the input file
#

outfile=$1
percent=$2

function d100 {
       d=$[ ( $RANDOM % 100  ) + 1 ]
}

if [ -a $outfile ]
then rm $outfile
fi

echo "Perturbing input file... this may take awhile."

done=0

read ch
while [ $done -eq 0 ]

do
	d100
	case $ch in
	'1')	if [ $d -lt $percent ]
		then
			echo -n 'x '>> $outfile
		else
			echo -n '1 '>> $outfile
		fi
		;;
	'0')	if [ $d -lt $percent ]
		then
			echo -n 'x '>> $outfile
		else
			echo -n '0 '>> $outfile
		fi
		;;
	'x')	echo -n 'x '>> $outfile
		;;
	'')	echo '' >> $outfile
		done=1
		;;
	esac
	read ch
done
