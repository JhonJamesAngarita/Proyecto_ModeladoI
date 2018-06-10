import gadget as gd
import numpy as np
import sys

sn = gd.Snapshot(sys.argv[1])  # le ponemos el primer alrgumento de bash, la cual es la ruta

Gas=[];
Dark=[];
Gas.append(sn.part0.pos[:,0]);
Gas.append(sn.part0.pos[:,1]);
Gas.append(sn.part0.pos[:,2]);
Gas.append(sn.part0.vel[:,0]);
Gas.append(sn.part0.vel[:,1]);
Gas.append(sn.part0.vel[:,2]);
Gas.append(sn.part0.mass);
Dark.append(sn.part1.pos[:,0]);
Dark.append(sn.part1.pos[:,1]);
Dark.append(sn.part1.pos[:,2]);
Dark.append(np.ones(np.shape(Dark)[1])*sn.masses[1]);
Gas=np.array(Gas).T
Dark=np.array(Dark).T
np.savetxt(str('SNAPSHOT'+sys.argv[2])+'Gas.csv', Gas, delimiter=',', fmt='%s')
np.savetxt(str('SNAPSHOT'+sys.argv[2])+'Dark.csv', Dark, delimiter=',', fmt='%s')

