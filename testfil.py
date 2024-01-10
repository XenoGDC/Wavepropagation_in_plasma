import SimMoviemaker as mkr
import numpy as np
import os
import Wave_Sim_reader as rd
import WavePropagationSim as sim
from matplotlib import pyplot as plt



# I,J = 600,600


# Matrix,hundred = rd.ReadBigSim(data_name)

# I,J,T = np.shape(Matrix)

# pMatrix = rd.gimidensity(I,J,0,'Blob',x0=int(I/3),y0=int(J/2),signy=41,peak=5e8)

# rd.plotplasmadens(I,J,pMatrix,'Blob','X-mode',B0=0.5)

mat1 = np.load('ne_1.npy')
mat2 = np.load('ne_2.npy')
mat3 = np.load('ne_3.npy')
x_list = np.load('x_list.npy')
y_list = np.load('y_list.npy')



print('Mat dim {dim}'.format(dim=np.shape(mat1)))
print('x_list dim {dim}'.format(dim=np.shape(x_list)))
print('y_list dim {dim}'.format(dim=np.shape(y_list)))


I,J = np.shape(mat1)


x_list = x_list[:I]
y_list = y_list[:J]

# fig,ax = plt.subplots()

# map1 = ax.pcolormesh(y_list,x_list,mat2)
# bar = plt.colorbar(map1)
# plt.show()
# plt.clf()


# fig,ax = plt.subplots()
# map2 = ax.pcolormesh(y_list,x_list,mat3)
# bar = plt.colorbar(map2)
# plt.show()
# plt.clf()

sigma = 0.02/rd.dx

I += 1
J += 1

# Write the directory for where you want to store your files here
data_collection_location = 'C:/Users/augus/Onedrive/skrivebord/DTU/Fagprojekt (bølgeudbredelse i plasme)'
os.chdir(data_collection_location)
angle_X = np.pi/2
name1x = 'Sim_for_ne_1_X'
name2x = 'Sim_for_ne_2_X'
name3x = 'Sim_for_ne_3_X'
name1o = 'Sim_for_ne_1_O'
name2o = 'Sim_for_ne_2_O'
name3o = 'Sim_for_ne_3_O'

# sim.WaveSim(200,I,J,sigma,field='Ey',CustomPMatrix=mat2,B0=[0,0,0.5],wave_polarity_angle=angle_X,CustomName=name2x)
# sim.WaveSim(200,I,J,sigma,field='Ey',CustomPMatrix=mat3,B0=[0,0,0.5],wave_polarity_angle=angle_X,CustomName=name3x)

# sim.WaveSim(200,I,J,sigma,field='Ez',CustomPMatrix=mat1,B0=[0,0,0.5],wave_polarity_angle=0,CustomName=name1o)
# sim.WaveSim(200,I,J,sigma,field='Ez',CustomPMatrix=mat2,B0=[0,0,0.5],wave_polarity_angle=0,CustomName=name2o)
# sim.WaveSim(200,I,J,sigma,field='Ez',CustomPMatrix=mat3,B0=[0,0,0.5],wave_polarity_angle=0,CustomName=name3o)


# name1x = 'Sim_for_ne_1_X'
# sim.WaveSim(300,I,J,sigma,field='Ey',CustomPMatrix=mat,B0=[0,0,0.5],wave_polarity_angle=0,CustomName=name)

# # newmat,hundred = rd.ReadBigSim(name)

# mkr.BigMovieMaker(name1x,mat1)
# mkr.BigMovieMaker(name2x,mat2)
# mkr.BigMovieMaker(name3x,mat3)
# mkr.BigMovieMaker(name1o,mat1)
# mkr.BigMovieMaker(name2o,mat2)
# mkr.BigMovieMaker(name3o,mat3)

mat2,hundred = rd.ReadBigSim(name3o,1)
# rd.GaussAnalyser(mat2,50,sig0=sigma,Vacuum=False)
mat1 = rd.SimReader('Misc Data/Ez_vacuum_dmpl50_dmp_i_7')
rd.GaussAnalyser(mat1,150)