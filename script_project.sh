cont=0
cd SNAPSHOTS/

#for N in Dwarf1_snapshot_100

for N in ./*.hdf5
do
let cont=cont+1
DataGas=*"$cont"Gas.csv 
DataDark=*"$cont"Dark.csv 

echo "galaxia: "$N
for k in 0.05 0.1 0.2 0.4
do
echo "dr = "$k
python ../script_project2.py $DataGas $DataDark $k $N  >> pendiente.txt
done
done

