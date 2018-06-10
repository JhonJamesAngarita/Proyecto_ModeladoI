import gadget as gd
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy import constants as const
from astropy import units as units
import scipy.optimize
import matplotlib.style as style


def Dif_Masa(masa):
    # input: masa, output: diferencial de masa  
    difmasa=[]
    for i in range(np.shape(masa)[0]-1):
        difmasa.append(masa[i+1]-masa[i])
    difmasa=np.array(difmasa)
    difmasa=np.concatenate((masa[:1],difmasa))
    return difmasa

def Part_Radio(r): 
    #input: radio, output: cantidad de particulas por radio r
    distancia=np.sqrt((sn.part1.pos**2).sum(axis=1))   
    cp=[]                                               
    for i in r:
        cp.append((distancia <= i).sum())  
    return cp

Gas = np.genfromtxt(sys.argv[1],delimiter=",")  #datos de Gas   columns: 0-->pos_x, 1-->pos_y, 2-->pos_z, 
                                                #               columns: 3-->vel_x, 5-->vel_y, 5-->vel_z, 6 -->masses
Dark = np.genfromtxt(sys.argv[2],delimiter=",") #datos de Dark  columns: 0-->pos_x, 1 -->pos_y, 2 -->pos_z, 3 -->masses 
dr  = float(sys.argv[3])
gal  = str(sys.argv[4])

distancia=np.sqrt((Dark[:,0]**2+Dark[:,1]**2+Dark[:,2]**2))  
r=np.arange(dr,4+dr,dr)
cp=[]                # cantidad de particulas por radio r
for i in r:
    cp.append((distancia <= i).sum())  

particulas=np.array(cp) 
masa=(particulas*Dark[0,3])  # Masa real

difmasa=Dif_Masa(masa)
rho=(difmasa*units.M_sun*1e10/dr*units.kpc)/(4*np.pi*(r*units.kpc)**2)  # Densidad real


velcir=np.sqrt(const.G*masa*units.M_sun*1e10/(r*units.kpc)).to('km/s')  # Vel. real

######################################################################
##### Datos de Gas, metodo observacional  ############################
######################################################################

x=[]
y=[]
condz = ((Gas[:,2] <= 1) & (Gas[:,2] >= -1))  # restriccion del plano z entre -1 y 1

x=Gas[condz,0]   # Eliminando la restriccion del plano z
y=Gas[condz,1]   # Eliminando la restriccion del plano z

vx=Gas[condz,3]   # Eliminando la restriccion del plano z
vy=Gas[condz,4]   # Eliminando la restriccion del plano z
unx=-y/(np.sqrt(x**2+y**2))   # vectores unitarios
uny=x/(np.sqrt(x**2+y**2))

##################################################
########## Vel calculada #########################
velcirc=[]   
for i in r:
    condr = ((x**2+y**2 <= i**2) & (x**2+y**2 >(i-dr)**2)) 
    velcirc.append(((vx*unx+vy*uny)*condr).sum()/(condr).sum())
velcirc=velcirc*units.km/units.s     # Vel circular observacional


##################################################
########## masa calculada ########################
masa1=[]    
masa1=(((velcirc**2)*r*units.kpc)/const.G).to(units.M_sun)
masa1=np.array(masa1)

##################################################
########## densidad calculada #########################
difmasa1=Dif_Masa(masa1)
rho1=(difmasa1/dr)/(4*np.pi*r**2)

##################################################
##### restriccion de los primeros 200 pc #########
cond_pc = (r > 0.2) 
rho1_pc=rho1[cond_pc]
rho_pc=rho[cond_pc]
r_pc=r[cond_pc]

##################################################
##### plot: rho real vs rho observacional ########
fig=plt.figure(figsize=(8,6))
graph=fig.add_subplot(111)
plt.plot(np.log(np.array(r_pc)),np.log(np.array(rho1_pc)),'.',np.log(np.array(r_pc)),np.log(np.array(rho_pc)),'*')
graph.legend(('Densidad observada', 'Densidad real'))
graph.set_xlabel('Radius [Kpc]')
plt.savefig(str(gal[2:4]+gal[-8:-5]+'_rho_dr'+ sys.argv[3])+'.svg')

x=np.log(r_pc[:3])
y=np.log(np.array(rho1_pc[:3]))
yr=np.log(np.array(rho_pc[:3]))

p0 = [0,1]
errfunc         = lambda p: np.ravel(p[0]*x+p[1]-y)
fitparams       = scipy.optimize.leastsq(errfunc,p0,full_output=1)[0]
fitparams       = list(fitparams)
m = fitparams[0]
b = fitparams[1]
print(m)
