import gsd.hoomd
import numpy as np
import pytransform3d.plot_utils as ppu
from distance3d import gjk, random, geometry, plotting, colliders
from ahrs import Quaternion 
import pytransform3d.transformations as pt
import pytransform3d.rotations as pr

# closest points distribution calculation for octahedron particles by using distance3d package

filename = "octa_vv_eqBL_1.0_N_80_ndump_35000"
path_temp = '/scratch4/tvo12/xzhan357/1_octahedron/v-v/BL=1.0/'

file_temp = path_temp + filename + '.gsd'
f = gsd.hoomd.open(name = file_temp, mode='r')
N_frame = len(f)        # number of frames
f.close()
N_start = 20        # discard first 20 samples

size = np.array([1.0,1.0,1.0])
box_L = np.array([0.0, 0.0, 0.0])

dist_list = []
point_location = []
q_initial = Quaternion([1.0, 0.0, 0.0, 0.0])
R_octa = 1.0/np.sqrt(2)
octa_verts =[(0.0, 0.0, R_octa), (R_octa, 0.0, 0.0), (0.0, R_octa, 0.0),
                (-R_octa, 0.0, 0.0), (0.0, -R_octa, 0.0), (0.0, 0.0, -R_octa)]    
for i in range (N_start,N_frame): #
    f = gsd.hoomd.open(name = file_temp, mode='r')
    frame = f[i]
    f.close()
    box_L[0] = frame.configuration.box[0]
    box_L[1] = frame.configuration.box[1]
    box_L[2] = frame.configuration.box[2]  
    N_length = frame.particles.N
    posi = frame.particles.position
    orien = frame.particles.orientation
    for j in range(0,N_length-1):#
        xyz_temp = posi[j]
        xyz_temp1 = posi[j+1]
        q_temp = orien[j]
        q_temp1 = orien[j+1]
        for k in range(3):       # unwrap coordinates
            while xyz_temp[k] - xyz_temp1[k] >= box_L[k]/2.0:
                xyz_temp1[k] = xyz_temp1[k] + box_L[k]
            while xyz_temp1[k] - xyz_temp[k] >= box_L[k]/2.0:
                xyz_temp[k] = xyz_temp[k] + box_L[k]
        
        octa_verts_0 = []
        octa_verts_1 = []
        for k in range(0, 6):
            temp_0 = pr.q_prod_vector(q_temp, octa_verts[k])
            temp_1 = pr.q_prod_vector(q_temp1, octa_verts[k])
            octa_verts_0.append(temp_0)
            octa_verts_1.append(temp_1)
         
        octa_verts_0 = np.array(octa_verts_0) + xyz_temp
        octa_verts_1 = np.array(octa_verts_1) + xyz_temp1
        octa_collider = colliders.ConvexHullVertices(octa_verts_0)
        octa_collider1 = colliders.ConvexHullVertices(octa_verts_1)
        dist, closest_point_box, closest_point_box1, _ = gjk.gjk(
        octa_collider, octa_collider1)        
        point = closest_point_box - xyz_temp        # calculate the vector relative to the center of the particle    
        point1 = closest_point_box1 - xyz_temp1
        q_inverse = Quaternion(q_temp).inverse
        q_1_inverse = Quaternion(q_temp1).inverse
        q_rotation = q_initial.product(q_inverse)
        q_1_rotation = q_initial.product(q_1_inverse)
        point_after_r = pr.q_prod_vector(q_rotation, point)         # rotate to the normal position 
        point_after_r1 = pr.q_prod_vector(q_1_rotation, point1)
        
        dis2center = np.sqrt(point_after_r[0]**2 + point_after_r[1]**2 + point_after_r[2]**2)       # distance calculation
        dis2center1 = np.sqrt(point_after_r1[0]**2 + point_after_r1[1]**2 + point_after_r1[2]**2)
        temp_array = np.insert(point_after_r, 0, dis2center)
        temp_array = [float(element) for element in temp_array]
        temp_array1 = np.insert(point_after_r1, 0, dis2center1)
        temp_array1 = [float(element) for element in temp_array1]
        dist_list.append(temp_array)
        dist_list.append(temp_array1)
    print(i)

filetemp =  path_temp + 'octa_vv_' + filename + "_closest.txt"       
output_file = open(filetemp, 'w')
for i in range(len(dist_list)):
    output_file.write(f"{dist_list[i][1]:.12f} {dist_list[i][2]:.12f} {dist_list[i][3]:.12f}\n")
       
        