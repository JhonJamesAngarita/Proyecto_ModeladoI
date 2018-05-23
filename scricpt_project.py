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

# Grafica de la densidad
difmasa=[]                                  #diferencial de masa_real
for i in range(np.shape(masa)[0]-1):
    difmasa.append(masa[i+1]-masa[i])
difmasa=np.array(difmasa)
difmasa=np.concatenate((masa[:1],difmasa))
dr=0.1
rho=(difmasa*units.M_sun*1e10/dr*units.kpc)/(4*np.pi*(r*units.kpc)**2)
plt.plot(np.log(np.array(r)),np.log(np.array(rho)),'.')
plt.show()

# Plot de velocidad
Radio=np.arange(0.1,4,0.1)
velcir=np.sqrt(const.G*masa*units.M_sun*1e10/(r*units.kpc)).to('km/s')
plt.plot(r,velcir,'.',color="r")
plt.show()

# Velocidades tangencial
x=[]
y=[]
condz = ((sn.part0.pos[:,2] <= 1) & (sn.part0.vel[:,2] >= -1) ) 
x=sn.part0.pos[condz,0]   # Eliminando la restriccion 
y=sn.part0.pos[condz,1]   # Eliminando la restriccion 
plt.plot(x,y,'.')
plt.show()

vx=sn.part0.vel[condz,0]   # Eliminando la restriccion 
vy=sn.part0.vel[condz,1]   # Eliminando la restriccion 
unx=-y/(np.sqrt(x**2+y**2))
uny=x/(np.sqrt(x**2+y**2))
velcirc=[]   
for i in r:
    condr = ((x**2+y**2 <= i) & (x**2+y**2 >i-0.1) ) 
    velcirc.append(((vx*unx+vy*uny)*condr).sum()/(condr).sum())
velcirc=velcirc*units.km/units.s
plt.plot(r,np.abs(velcirc),'*',r,velcir,'.')
plt.show()

masa1=[]     #calculo de masa1

masa1=(((velcirc**2)*r*units.kpc)/const.G).to(units.M_sun)
masa1=np.array(masa1)
difmasa1=[]
for i in range(np.shape(masa1)[0]-1):
    difmasa1.append(masa1[i+1]-masa1[i])
difmasa1=np.array(difmasa1)
difmasa1=np.concatenate((masa1[:1],difmasa1))

rho1=(difmasa1/dr)/(4*np.pi*r**2)
plt.plot(np.log(np.array(r)),np.log(np.array(rho1)),'.',np.log(np.array(r)),np.log(np.array(rho)),'*')
plt.show()
