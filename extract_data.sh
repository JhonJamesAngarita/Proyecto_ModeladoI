etiqueta=0
cd SNAPSHOTS/
for N in ./*.hdf5
do
let etiqueta=etiqueta+1
python ../extract.py $N $etiqueta
done

