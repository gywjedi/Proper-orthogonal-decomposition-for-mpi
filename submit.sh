st=40
end=60

for (( x=$st ; x<$end ; x++ ))
do
    cidx=`printf %03i $x`

    if [ $x == $st ] ; then
        previd=`sbatch q_batch.$cidx.sh | awk '{print $4}'`
    else
        previd=`sbatch -d afterok:$previd q_batch.$cidx.sh  | awk '{print $4}'`
    fi


done
