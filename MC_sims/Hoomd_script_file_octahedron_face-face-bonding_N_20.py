###############
### Modules ### # hoomd v3.11
################

import hoomd
from hoomd import hpmc
import os
import numpy as np
import random
import gsd.hoomd
import datetime

# The script file to perform a simulation of face-to-face bonding octahedron with a chain length of 20

#   Find the path of the simulation file
path = os.path.dirname(os.path.abspath(__file__))
path = path + "/"
####################################################    
eq_Bond_L = 0.1     # Define eqilibrium bond length       
N_length = 20           # Define polymer chain length   

# calculate octahedron edge length       
cube_half_L = 0.5    
R_octa = 1.0/np.sqrt(2)                 
dia_dis = np.sqrt(2)

# dump samples every 1700 steps, collect 2000 samples totally
ndump = 1700                  
total_run = ndump*2000   
# define simulation box length
Lbox = (dia_dis)*N_length*2 

# define simulation snapshot
pts_init = []
typeid_init = []
charges_init = []
types_str = []
orientation = [1,0,0,0]
for i in range(1,N_length+1):
    if eq_Bond_L < dia_dis:
        x_temp = -Lbox/4+i*(dia_dis + eq_Bond_L)
    else:
        x_temp = -Lbox/4+i*eq_Bond_L
    y_temp = 0
    z_temp = 0
    xyz_add = [x_temp,y_temp,z_temp]
    pts_init.append(xyz_add)
    typeid_init.append(i-1)
    types_str.append(str(i-1))
    charges_init.append(i-1)
snapshot = gsd.hoomd.Snapshot()
snapshot.particles.N = len(pts_init)
snapshot.particles.position = pts_init
snapshot.particles.typeid = typeid_init      
snapshot.configuration.box = [Lbox, Lbox, Lbox, 0, 0, 0]
snapshot.particles.types = types_str
snapshot.particles.orientation = N_length*orientation
snapshot.particles.charge = charges_init


r_cut_dis = 5.0 # cut-off distance for user-defined potential
#************************* user-defined bonding potential (harmonic bond) *********************************
harmonic_bond = ''' float diff_ij = abs(charge_i - charge_j);
                    float diff_typeid_ij = abs(float(type_i) - float(type_j));
                    float rsq = dot(r_ij, r_ij);
                    if (int(diff_ij) == 1 && int(diff_typeid_ij) == 1) {
                        float bondenergy = 3.0f/2.0f/0.01f*rsq;  //  equilibrium bond length 
                        return bondenergy;
                    }
                    else {
                        return 0.0f; 
                    }
                '''
patch = hoomd.hpmc.pair.user.CPPPotentialUnion(
    r_cut_constituent=r_cut_dis,                     # cut-off
    r_cut_isotropic=0.0,
    code_constituent=harmonic_bond,
    code_isotropic='',
    param_array_constituent=[],
    param_array_isotropic=[],
)
# define the particle and bonding sites as patchy particle, using charges to define the bonding sequence
for i in range(0,N_length):
    patch.positions[str(i)] = [(-R_octa/3.0, -R_octa/3.0, -R_octa/3.0), 
                                (R_octa/3.0, R_octa/3.0, R_octa/3.0)]
    patch.orientations[str(i)] = [(1, 0, 0, 0), (1, 0, 0, 0)]
    patch.typeids[str(i)] = [i, i]
    patch.diameters[str(i)] = [0, 0]
    patch.charges[str(i)] = [i*2+1, i*2+2]


# define particle shape in the simulation
mc = hoomd.hpmc.integrate.ConvexPolyhedron(default_d=eq_Bond_L, default_a=0.5)
octa_verts =[(0.0, 0.0, R_octa), (R_octa, 0.0, 0.0), (0.0, R_octa, 0.0),
                (-R_octa, 0.0, 0.0), (0.0, -R_octa, 0.0), (0.0, 0.0, -R_octa)]                               
for i in range(0,N_length):
    mc.shape[str(i)] = dict(vertices=octa_verts)
mc.pair_potential = patch

# create the simulation
cpu = hoomd.device.CPU()
sim = hoomd.Simulation(device=cpu, seed=2501)
sim.operations.integrator = mc
sim.create_state_from_snapshot(snapshot)
#output initial cfg
hoomd.write.GSD.write(state=sim.state, mode='wb', filename=path + 'initial20.gsd')
#*****************************************************************


### Define logger for tps printing ###
logger = hoomd.logging.Logger(categories=['scalar', 'string'])
logger.add(sim, quantities=['timestep', 'tps'])
class Status():
    def __init__(self, sim):
        self.sim = sim
    @property
    def seconds_remaining(self):
        try:
            return (self.sim.final_timestep - self.sim.timestep) / self.sim.tps
        except ZeroDivisionError:
            return 0
    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining))
status = Status(sim)
table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=ndump),
                          logger=logger)
sim.operations.writers.append(table)
# Shape logging
logger_shape = hoomd.logging.Logger()
logger_shape.add(mc,quantities=['type_shapes'])
######################################



# Dump samples to gsd files
file_name = "eqBL_" + str(eq_Bond_L) + "_N_" + str(N_length) + "_ndump_" + str(ndump) + '.gsd'
gsd_writer = hoomd.write.GSD(filename=path + file_name,                 
                             trigger=hoomd.trigger.Periodic(ndump),
                             mode='wb',
                             logger=logger_shape)
sim.operations.writers.append(gsd_writer)
# run the simulation
sim.run(total_run)          
