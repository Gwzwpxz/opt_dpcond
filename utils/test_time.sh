for file in /Users/gaowenzhi/Desktop/public/opt_dpcond/datasets/suitesparse/sdp/*-M.dat-s
do
echo $file
./DSDP6 $file > $file.log
done
