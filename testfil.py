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

mat = np.load('ne_1.npy')
x_list = np.load('x_list.npy')
y_list = np.load('y_list.npy')



print('Mat dim {dim}'.format(dim=np.shape(mat)))
print('x_list dim {dim}'.format(dim=np.shape(x_list)))
print('y_list dim {dim}'.format(dim=np.shape(y_list)))

I,J = np.shape(mat)
x_list = x_list[:I]
y_list = y_list[:J]

# fig,ax = plt.subplots()

# map1 = ax.pcolormesh(y_list,x_list,mat)
# bar = plt.colorbar(map1)
# plt.show()

sigma = 0.02/rd.dx
I += 1
J += 1


data_collection_location = 'C:/Users/augus/Onedrive/skrivebord/DTU/Fagprojekt (bølgeudbredelse i plasme)'
os.chdir(data_collection_location)
# data_folder_name = 'Blob_n_5e18_width_var'
# data_name = 'X-mode_blob_test_sig_5e+1841_n_'
# os.chdir(data_folder_name)
# angle = np.pi/2
name = 'Sim_for_ne_1_X'
# sim.WaveSim(400,I,J,sigma,field='Ey',CustomPMatrix=mat,B0=[0,0,0.5],wave_polarity_angle=angle,CustomName=name)

# newmat,hundred = rd.ReadBigSim(name)

mkr.BigMovieMaker(name)