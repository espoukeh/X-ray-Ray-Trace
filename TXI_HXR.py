import os
import sys
import subprocess
import numpy as np

#define propagation direction
source_angle = 0.007*4 #rad note 2 mirrors with 2theta angle
px = np.sin(source_angle) #easy way to get a prop_vector
pz = np.cos(source_angle) 

#Using BTM3L1 as source
source_x = -1.0023862 # source position [m]
source_z = 749.6411940 # source position [m]
source_size = 15.97*1.0E-3 / 2.0 # source diameter [m]
source = "./make_source -p [{},0,{}] -v [{},0,{}] -d [{},0,0] -s [{},0,0] -n 10000 -r uniform".format(source_x,source_z,px,pz,source_size, 0.00015)
subprocess.call(source + " -o in.dat", shell=True)
subprocess.call("./print_rays -i in.dat -o in.csv", shell=True)

#BTM4L1 as first collimator
c1_x = -0.8282077 # [m]
c1_z = 755.8602290 # [m]
c1_radius = 11.49*1E-3/2.0 # m
c1 = "./circ_aperture --C=[{},0,{}] --R={} --N=[{},0.0,{}]".format(c1_x,c1_z,c1_radius,px,pz)
subprocess.call(c1 + " -i in.dat -o c1.dat", shell=True)
subprocess.call("./print_rays -i c1.dat -o c1.csv", shell=True)

#BTM5L1 as second collimator
c2_x = -0.4864914 # [m]
c2_z = 768.0611940 # [m]
c2_radius = 11.49*1E-3/2.0 # m
c2 = "./circ_aperture --C=[{},0,{}] --R={} --N=[{},0.0,{}]".format(c2_x,c2_z,c2_radius,px,pz)
subprocess.call(c2 + " -i c1.dat -o c2.dat", shell=True)
subprocess.call("./print_rays -i c2.dat -o c2.csv", shell=True)

#MR3L1 Horizontal KB Mirror
M_z = 771.665 #mirror Z center position [m]
M_x = -0.385558472 # mirror X center position [m]
M_theta = source_angle + 0.007 # radians
M_Lx = 1.0 # mirror length [m]
M_Ly = 0.05  # mirror width [m] a bit arbutary for 2D but needed
#calculate bounding points of mirror		
P0_x = M_x - ((M_Lx/2.0)*np.sin(M_theta))
P0_y = (M_Ly/2.0)
P0_z = M_z - ((M_Lx/2.0)*np.cos(M_theta))

P1_x = M_x + ((M_Lx/2.0)*np.sin(M_theta))
P1_y = (M_Ly/2.0)
P1_z = M_z + ((M_Lx/2.0)*np.cos(M_theta))

P2_x = M_x - ((M_Lx/2.0)*np.sin(M_theta))
P2_y = -1.0*(M_Ly/2.0)
P2_z = M_z - ((M_Lx/2.0)*np.cos(M_theta))

mirror = "./flat_mirror --P0=[{},{},{}] --P1=[{},{},{}] --P2=[{},{},{}]".format(P0_x,P0_y,P0_z,P1_x,P1_y,P1_z,P2_x,P2_y,P2_z)
subprocess.call(mirror + " -i c2.dat -o m1.dat", shell=True)
subprocess.call("./print_rays -i m1.dat -o m1.csv", shell=True)

#beamdump location
dump_z = 778.667000  #1 m downstream of photon terminator
dump = "./propagate_to --P0=[0,0,{}] --P1=[1,0,{}] --P2=[0,1,{}] ".format(dump_z,dump_z,dump_z)
subprocess.call(dump + " -i m1.dat -o d.dat", shell=True)
subprocess.call("./print_rays -i d.dat -o d.csv", shell=True)

