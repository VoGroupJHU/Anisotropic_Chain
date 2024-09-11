###############
### Modules ###
###############

'''
Paper: "Elucidating the Interplay Between Entropy-Driven and Patch-Mediated 
Bonding in Directing Nanoscale Assemblies"
Sample Code for MD run
Authors: Kireeti Akkunuri, Xiangyu Zhang, Thi Vo
Dated: September 2024
'''

from hoomd import *
from hoomd import md
from hoomd import dem
from hoomd import deprecated
import os
import numpy as np
import random
import glob

###############
###############

#####################################
### Read in initial configuration ###
#####################################

# Read in file
file_name = 'init_bond.txt'
f = open(file_name,'r')
message = f.readlines()
bond_config = []
bond_type_config = []
for line in message:

  # Read file
  tmp_line = line
  bond_type,bond_i,bond_j = tmp_line.split()
  bond_tmp = [int(bond_i),int(bond_j)]  

  # Add
  bond_config.append(bond_tmp)
  bond_type_config.append(bond_type)

# Close file
f.close()

# Read in file - monomer core
file_name = 'verts.txt'
f = open(file_name,'r')
message = f.readlines()
vertices_A = []
for line in message:

  # Read file
  tmp_line = line
  x_tmp,y_tmp,z_tmp,v_tmp,Ixx_A,Iyy_A,Izz_A = tmp_line.split()
  pts_tmp = [float(x_tmp),float(y_tmp),float(z_tmp)]  

  # Add
  vertices_A.append(pts_tmp)  

# Close file
f.close()
vertices_A=vertices_A-np.mean(vertices_A,axis=0)

# Volume
v_A = float(v_tmp)

# Moment of inertia
Ixx_A = float(Ixx_A)
Iyy_A = float(Iyy_A)
Izz_A = float(Izz_A)

# Read in file - ghost monomers
file_name = 'verts_rigid.txt'
f = open(file_name,'r')
message = f.readlines()
pts_rigid_A = []
q_rigid_A = []
type_rigid_A = []
for line in message:

  # Read file
  tmp_line = line
  type_tmp,x_tmp,y_tmp,z_tmp,qs_tmp,qx_tmp,qy_tmp,qz_tmp = tmp_line.split()
  pts_tmp = [float(x_tmp),float(y_tmp),float(z_tmp)]  
  q_tmp = [float(qs_tmp),float(qx_tmp),float(qy_tmp),float(qz_tmp)]  

  # Add
  pts_rigid_A.append(pts_tmp)  
  q_rigid_A.append(q_tmp)
  type_rigid_A.append('B')
  
# Close file
f.close()

# Read in file
file_name = 'init.txt'
f = open(file_name,'r')
message = f.readlines()
pts_config = []
q_config = []
type_config = []
body_config = []
moment_config = []
mass_config = []
for line in message:

  # Read file
  tmp_line = line
  type_tmp,body_tmp,x_tmp,y_tmp,z_tmp,qs,qx,qy,qz,sigma_A,sigma_Aavg,sigma_B,Lx,Ly,Lz = tmp_line.split()
  pts_tmp = [float(x_tmp),float(y_tmp),float(z_tmp)]    
  q_tmp = [float(qs),float(qx),float(qy),float(qz)]

  # Add
  pts_config.append(pts_tmp)
  type_config.append(int(type_tmp))
  body_config.append(int(body_tmp))
  q_config.append(q_tmp)
  if int(type_tmp) == 0:
    moment_factor = 3.0     # to make certain types more rigid to rotation i.e. 'gentler' rotation
    monomer_add = [moment_factor*Ixx_A,moment_factor*Iyy_A,moment_factor*Izz_A]
    moment_config.append(monomer_add)
    mass_config.append(1.0)
  else:
    Isph = 0.4*(float(sigma_B)/2.0)**2.0
    moment_add = [Isph,Isph,Isph]
    moment_add = [1.0,1.0,1.0]
    moment_config.append(moment_add)
    v_ghost = (4.0/3.0)*np.pi*(float(sigma_B)/2.0)**3.0
    mass_config.append(1.0)

# Close file
f.close()

# Box size
Lx = float(Lx)
Ly = float(Ly)
Lz = float(Lz)
Lx = Lz
Ly = Lz

# Particle size
sigma_A = float(sigma_A)

# Effective particle size
sigma_Aavg = float(sigma_Aavg)

# Particle size
sigma_ghost = float(sigma_B)

# Count A
indx_A = np.where(np.array(type_config) == 0)
N_A = len(indx_A[0])

# Count B
indx_B = np.where(np.array(type_config) == 1)
N_B = len(indx_B[0])
 
#####################################
#####################################

###################
### Set up sims ###
###################

# Box factor
box_factor = 1

# Load HOOMD
context.initialize("--mode=cpu");

# Domain decomposition
## manually change each split. Better for the system to be even
comm.decomposition(nx=2,ny=2,nz=2)

# Make snapshot
# Note #
# A: shape A
# B: ghost linker
########
# Read existing GSD files
file_list_1 = glob.glob('sample_single_part*.gsd')
if not file_list_1:
    gsd_already_exists = False
    gsd_output = 'sample_single_part1.gsd'
else:
    gsd_already_exists = True
    file_list_1.sort()
    file_found = os.path.basename(file_list_1[-1])
    gsd_output = 'sample_single_part'+ str(len(file_list_1)+1) + '.gsd'

########

if gsd_already_exists is False:
    current_step = 0
    snapshot = data.make_snapshot(N=len(type_config),
                box=data.boxdim(Lx=box_factor*Lx, Ly=box_factor*Ly, Lz=box_factor*Lz),
                particle_types=['A','B'],
                bond_types=['ghost-ghost'])

    # Define particle positions
    snapshot.particles.resize(len(pts_config));
    snapshot.particles.position[:] = pts_config
    snapshot.particles.orientation[:] = q_config
    snapshot.particles.moment_inertia[:] = moment_config
    snapshot.particles.typeid[:] = type_config
    snapshot.particles.mass[:] = mass_config
    snapshot.particles.body[:] = body_config

    # Bonds
    snapshot.bonds.resize(len(bond_config));
    snapshot.bonds.group[:] = bond_config
    snapshot.bonds.typeid[:] = bond_type_config

    # Read in system
    system = init.read_snapshot(snapshot);

    # Replicate
    nrep = 1
    system.replicate(nx=nrep,ny=nrep,nz=nrep)


    # Calculate number of chains
    snapshot_now = system.take_snapshot()
    N_final = snapshot_now.particles.N
    N_init = len(pts_config)
    N_chains = int(N_final/N_init)
    N= int(len(body_config)/3) # number of monomers per chain

    # Store chain index into charge
    charge_config = []
    for i in range(0,N_chains):
        for j in range(0,N):
            for k in range(0,3):
                charge_config.append(i)

    snapshot_now.particles.charge[:] = charge_config
    system.restore_snapshot(snapshot_now)
else:
    system = init.read_gsd(file_found,frame=-1)
    snapshot_now = system.take_snapshot()
    current_step = get_step()
    


# Check rigid body
rigid = md.constrain.rigid()
rigid.set_param('A',types=type_rigid_A, positions=pts_rigid_A); # prototype ghost vertices
rigid.validate_bodies()

# Define current z direction
Lz_current = system.box.Lz
Lx_current = system.box.Lx
Ly_current = system.box.Ly

# Compute volumes
v_ghost = (4.0/3.0)*np.pi*(sigma_ghost/2.0)**(3.0)
V_chain = v_A*N_A + v_ghost*N_B     # Single chain

# Define final monomer density
rho_target = 0.85   # Volume packing fraction
rho_target_sigma = rho_target/(sigma_Aavg**3.0)

# Effective number of monomer - volume average
Navg = N_A*(v_A*N_A/V_chain) + N_B*(v_ghost*0.5*N_B/V_chain)

# Box size
Lfinal = (Navg/(rho_target_sigma))**(1.0/3.0)

 
###################
###################

###############################
### Define potential params ###
###############################

# Create groups
all = group.all();

# Neighbor list
nl = md.nlist.cell()

# Create groups
all = group.all();

# Neighbor list
nl = md.nlist.cell()

# LJ parameters
kbt = 0.9
eps_rep = 1.0
eps_att = 1.0

# Distances ghost
sigma_Aghost = 0.5*(sigma_A+sigma_ghost)

# Distance cutoff repulsive
cutoff_rep = 2.0**(1.0/6.0)
rcut_AA_rep = cutoff_rep*sigma_A*(3.0)

# Distace cutoff repulsive ghost
rcut_Aghost_rep = cutoff_rep*sigma_Aghost
rcut_ghost_rep = cutoff_rep*sigma_ghost

# Distance cutoff attractive
cutoff_rep = 3.0
rcut_AA_att = cutoff_rep*sigma_A

# Distace cutoff attractive ghost
rcut_Aghost_att = cutoff_rep*sigma_Aghost
rcut_ghost_att = cutoff_rep*sigma_ghost

# Shape Potential
buffer_factor = 1.15
lj = md.pair.alj(r_cut=buffer_factor*rcut_AA_rep,nlist=nl)
lj.set_params(mode='shift')
lj.shape['A'] = {'vertices': list(vertices_A)}
lj.shape['B'] = {'rounding_radii': [sigma_ghost*0.5]}

# Set interactions
# Pure
lj.pair_coeff.set('A','A', epsilon=eps_rep,sigma_i=sigma_A,sigma_j=sigma_A,alpha=0, r_cut=buffer_factor*rcut_AA_rep)
lj.pair_coeff.set('B','B', epsilon=eps_rep,sigma_i=sigma_ghost,sigma_j=sigma_ghost,alpha=0, r_cut=buffer_factor*rcut_ghost_rep)
# Ghost
lj.pair_coeff.set('A','B', epsilon=eps_rep,sigma_i=sigma_A,sigma_j=sigma_ghost,alpha=0, r_cut=buffer_factor*rcut_AA_rep)

# Define bond potential
kspring_gg = 30.0
ro = 1.5*1.15

harmonic = md.bond.fene()
harmonic.bond_coeff.set('ghost-ghost', k=kspring_gg, r0=ro*sigma_ghost, sigma=sigma_ghost, epsilon= eps_rep)

###############################
###############################

####################################
### Set up integration and dumps ###
####################################

# Group
full = group.all()
groupA = group.type(type='A')
groupB = group.type(type='B')

# Set up the system
## Set integration step size using 'dt'
md.integrate.mode_standard(dt=2.5E-4)
integrator = md.integrate.nvt(group=groupA, kT=kbt, tau=0.1)
## Change the 'seed' during each run to randomize initial velocities 
integrator.randomize_velocities(seed=2355)
## total simulation run time
t_final = 10E7

# Dumps
gsd_name = 'sample.gsd'
gsd_dump = dump.gsd(
    filename=gsd_output,
    period=5E3,
    phase=0,
    group=full,
    overwrite=True,
    dynamic=['attribute','property','momentum','topology'])
gsd_dump.dump_shape(lj)
util.quiet_status()
steps_spring = int(1E4)
ro_space_array = np.linspace(25.0*sigma_ghost,sigma_ghost,int(501))

# First, slowly relax the overstretched bonds to the desired bond length 
if int(current_step/steps_spring) < len(ro_space_array):
    for i in range(int(current_step/steps_spring),len(ro_space_array)):
    	print('Run: ',i)
    	harmonic.bond_coeff.set('ghost-ghost', k=kspring_gg, r0=ro*ro_space_array[i], sigma=sigma_ghost, epsilon= eps_rep)
    	run(steps_spring,quiet=True)
    util.unquiet_status()

# Then, run the MD simulation to relax the whole chain
run_upto(t_final)


####################################
####################################
