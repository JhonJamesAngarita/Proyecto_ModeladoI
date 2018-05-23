import gadget as gd
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as units

filename='SNAPSHOTS/Dwarf1_snapshot_100.hdf5'
sn=gd.Snapshot(filename)

distancia=np.sqrt((sn.part1.pos**2).sum(axis=1))    # calculo de la distancia de la particula
cp=[]                                               # cantidad de particulas por radio r
r=np.arange(0.1,4.1,0.1)
for i in r:
    cp.append((distancia <= i).sum())  

particulas=np.array(cp)                        
masa=(particulas*sn.masses[1])
plt.plot(r,masa*units.M_sun*1e10,'*',color="g")
plt.show()
