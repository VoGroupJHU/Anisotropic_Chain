import gsd.hoomd
import numpy as np
import math

# The file to calculate chain end-to-end distance

filename = "eqBL_0.1_N_20_ndump_1700"
path_temp = '/scratch4/tvo12/xzhan357/1_octahedron/f-f/BL=0.1/'

file_temp = path_temp + filename + '.gsd'
f = gsd.hoomd.open(name = file_temp, mode='rb')
N_frame = len(f)        # number of frames
N_start = 10        # discard first 10 samples

step = [0 for y in range(N_frame)]
Re_array = []
Re2_array = []
step = []

check = 0
counting = 0
for i in range (N_start,N_frame, 1):
    frame_temp = f[i]
    Lx = frame_temp.configuration.box[0]
    Ly = frame_temp.configuration.box[1]
    Lz = frame_temp.configuration.box[2]
    N_length = frame_temp.particles.N
    head = frame_temp.particles.position[0]
    tail = frame_temp.particles.position[N_length-1]
    Re = tail - head
    Re[0] = Re[0] - Lx * round(Re[0]/Lx, 0)
    Re[1] = Re[1] - Ly * round(Re[1]/Ly, 0) 
    Re[2] = Re[2] - Lz * round(Re[2]/Lz, 0) 
    norm = Re[0]**2 + Re[1]**2 + Re[2]**2
    Re_array.append(float(np.sqrt(norm)))
    Re2_array.append(float(norm))
    step.append(float(frame_temp.configuration.step))
    
    check = check + Re2_array[counting]
    counting = counting + 1

print(check/(counting), N_start, N_frame, counting)

file_temp = path_temp + filename + '_Re.txt'
output_file = open(file_temp, 'w') 
output_file.write('Re2 average and std are:' + str(np.average(Re2_array)) + ' ' + str(np.std(Re2_array)) + '\n')
output_file.write('Re average and std are:' + str(np.average(Re_array)) + ' ' + str(np.std(Re_array)) + '\n')
for i in range(len(Re_array)):
    output_file.write(str(step[i]) + ' ' + str(Re_array[i]) + '\n')

output_file.close()    
f.close()
