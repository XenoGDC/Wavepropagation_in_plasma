import numpy as np
import matplotlib.pyplot as plt
import time
import os
# import csv

# Constants
c = 2.99792458e8
mu = 4e-7 * np.pi
epsilon = 8.8541878e-12 # vacuum permebility
e_m = 9.1093837015e-31
e_q = 1.6022e-19

# Sets our parameters for angle frequency, dx, dy and dt.
# O = 2 * np.pi * 2.45e9
O = 2 * np.pi * 28e9
#x,y, and t is dx, dy, and dt
x = np.pi * c / (32 * O)
y = np.pi * c / (32 * O)
t = x / c / 2

# The function takes [the amount of time intervals, the x-size, y-size, dt, dx, dy, 
# angle frequency, time interval step size,inputwave width,analyse gaussian spread,plasma density type]
def WaveSim(TimeSize: int, Xsize: int, Ysize: int, Gauss_beam_FWHM: int, angle_frequency = O, 
            dt = t, dx = x, dy = y, TimeInterval = 10, 
            Plasma_type = 'Vacuum',wave_polarity_angle = 0, B0 = [0,0,0.1], Linear_angle = 0, dampening_length:int = 50
            ,dampening_intensity:float = 7,blobintensity:float = 0.9e20,blobradius: float = 40,
            blobposition:list = [0,0], CustomName = None,field:str = 'Ez',Give_Pdensity: bool = False,
            Linear_cutoff: int = 0,CustomPMatrix = None):

    B0x = B0[0]
    B0y = B0[1]
    B0z = B0[2]
    T = TimeSize
    I = Xsize
    J = Ysize
    omega = angle_frequency
    intval = TimeInterval
    sig0 = Gauss_beam_FWHM
    Pmode = Plasma_type
    angle = wave_polarity_angle

    fallback = os.getcwd()


    
    if not type(CustomName) == str:
            mapname = str(field + '_matrix_for_' + Pmode)
    else:
            mapname = str(CustomName)

    os.mkdir(mapname)
    os.chdir(mapname)

    #Sets the lowest
    if intval <= 1:
        intval = 2
    
    # sig = J / 100 * dy
    sig = sig0*dy
    n = omega**2/(e_q**2/e_m/epsilon)
    lam = c/omega
    try:
        if CustomPMatrix == None:
            if Pmode == 'Vacuum':
                scale = 0

            if Pmode == 'Constant':
                n_e = 0.9e17
                scale = n_e
            
            elif Pmode == 'Linear':
                # The scale is used as a value to multiply the plasma density into, 
                # in case the plasma doesn't have a constant density.
                scale = np.zeros([I-1,J-1])
                # n_e = 1.5*np.linspace(0,n,J-1)
                # lintheta = np.pi/4
                lintheta = Linear_angle
                if Linear_cutoff == 0:
                    Linear_cutoff = I/2

                # n0 = 0.9e19/(400*np.cos(lintheta)+400*np.sin(lintheta))

                def lindensityincrease(lintheta,x0,y0,n0,matrix):
                    n = n0/(x0*np.cos(lintheta)+y0*np.sin(lintheta))
                    scaletemp = np.ones(matrix.shape)
                    for xx in range(len(matrix[:,0])):
                        for yy in range(len(matrix[0,:])):
                            scaletemp[xx,yy] = n*(xx*np.cos(lintheta)+yy*np.sin(lintheta))
                    return scaletemp
                
                scale += lindensityincrease(lintheta,Linear_cutoff,J/2,0.97e19,scale)
                
                # for ii in range(I-1):
                #     scale[ii,:] = scale[ii,:]*n_e[ii]
            

            elif Pmode == 'Blob':
                # n0 = 0.9e20
                # x0 = J/2
                # y0 = I/2
                # signy = 40
                n0 = blobintensity
                scale = np.zeros([I-1,J-1])
                x = np.linspace(0,I,I)
                y = np.linspace(0,J,J)
                def blobmaker(x0,y0,signy,peak):
                    scaletemp = np.zeros([I-1,J-1])
                    for i in range(I-1):
                        for j in range(J-1):
                            scaletemp[i,j] = peak*np.exp(-(x[i]-x0)**2/signy**2-(y[j]-y0)**2/signy**2)
                    return scaletemp
                

                scale += blobmaker(blobposition[0],blobposition[1],blobradius,blobintensity)


            if Pmode == 'Blob' or Pmode == 'Linear':
                map1 = plt.pcolormesh(scale[:,:],cmap='seismic')
                cbar = plt.colorbar(map1)
                cbar.set_label('Plasma density')
                plt.draw()
                # plt.show()
    except:
        scale = np.zeros([I-1,J-1])

        scale += CustomPMatrix[:I,:J]
    
    
    # The plasma electron frequency
    Ope = (scale*e_q**2/e_m/epsilon)**0.5
    
    # Setting up the fields as empty arrays to hold values.
    Ex = np.zeros([intval + 1, I + 1, J + 1])
    Ey = np.zeros([intval + 1, I + 1, J + 1])
    Ez = np.zeros([intval + 1, I + 1, J + 1])

    Bx = np.zeros([intval + 1, I + 1, J + 1])
    By = np.zeros([intval + 1, I + 1, J + 1])
    Bz = np.zeros([intval + 1, I + 1, J + 1])

    Jx = np.zeros([intval + 1, I + 1, J + 1])
    Jy = np.zeros([intval + 1, I + 1, J + 1])
    Jz = np.zeros([intval + 1, I + 1, J + 1])

    #Parameters for counting when to take the gaussian spreading.
    tdot = 0
    count = False
    count2 = False
    start0 = True
    twave = 0

    # Stores the time
    t0 = time.time()
    t000 = time.time()
    
    # Dampening boundary conditions:
    # dmp_l sets the distance from the boundary that the dampening begins
    dmp_l = dampening_length
    dmp_mx = np.ones([1-dmp_l+I,dmp_l])*np.linspace(1,dmp_l,dmp_l)
    dmp_my = np.flip(np.transpose(np.ones([1-2*dmp_l+J,dmp_l])*np.linspace(1,dmp_l,dmp_l)))
    dmp_mx_t = np.transpose(dmp_mx[:dmp_l,:])
    dist_m = (dmp_mx[:dmp_l,:]**2 + dmp_mx_t**2)**0.5
    
    dmp_cn = np.zeros([dmp_l,dmp_l])
    for ii in range(dmp_l):
        for jj in range(dmp_l):
            dmp_cn[ii,jj] = min([dmp_mx[0,ii],dmp_my[jj,0]])
    dmp_cnr = np.flip(dmp_cn,0)
    dmp_cnl = np.flip(dmp_cnr,1)
    
    # dmp_i sets the intensity of the dampening, with higher values dampening over longer periods.
    dmp_i = dmp_l / dampening_intensity

    def dmp_model(dmp_mat,dmp_i):
        return (1-np.exp((1-dmp_mat) / dmp_i))

    Constant_save = list([dx,dy,I,J,dt,Ope,sig0])
    
    y_ax = np.linspace(0,J,J)

    T_big = int(round(T/100))

        # Just a little something if we want to save what kind of plasma density we have
    if Give_Pdensity:
        place = np.ones([I-1,J-1])
        Pdensity_matrix = place*scale
        Pname = str('Pdensity_' + txtname)
        np.save(Pname,Pdensity_matrix)


    for kk in range(T_big):
    # For every 100 steps. We define the matrix like this, so we don't use all our DRAM when we are storing the larger matrix we use later for analysis
        for k in range(100):
            # Remembers which time we are in the simulation for use in the sinus function
            t00 = k * intval + 100*kk
            

            # Simulates the wave equations for the next interval
            for t in range(intval):
                # Sums the earlier known time placement with the loops additional time placement to determine the current time placement.
                b = int(t00 + t)

                # Electric field
                Ex[t + 1, 1:-1, 1:-1] = Ex[t, 1:-1, 1:-1] + c ** 2 * dt * (
                            (Bz[t, 1:-1, 1:-1] - Bz[t, 1:-1, :-2]) / dy - mu * Jx[t, 1:-1, 1:-1])
                Ey[t + 1, 1:-1, 1:-1] = Ey[t, 1:-1, 1:-1] - c ** 2 * dt * (
                            (Bz[t, 1:-1, 1:-1] - Bz[t, :-2, 1:-1]) / dx + mu * Jy[t, 1:-1, 1:-1])
                Ez[t + 1, 1:-1, 1:-1] = Ez[t, 1:-1, 1:-1] + c ** 2 * dt * ((By[t, 1:-1, 1:-1] - By[t, :-2, 1:-1]) / dx - (
                            Bx[t, 1:-1, 1:-1] - Bx[t, 1:-1, :-2]) / dy - mu * Jz[t, 1:-1, 1:-1])

                # Gaussian distribution of sin-wave across y=0 boundary
                sinf = np.sin(omega * dt * (b))
                
                Ez[t + 1, 1, :-1] = np.cos(angle)*np.exp((-((y_ax - J / 2)*dy / sig) ** 2)) * sinf
                Ey[t + 1, 1, :-1] = np.sin(angle)*np.exp((-((y_ax - J / 2)*dy / sig) ** 2)) * sinf



                # Boundary conditions for electric field
                
                # #Reflecting boundary condition
                # Ey[t+1,-1,:] = 0
                # Ey[t+1,:,-1] = 0
                # Ey[t+1,:, 0] = 0
                # Ez[t+1,-1,:] = 0
                # Ez[t+1,:,-1] = 0
                # Ez[t+1,:, 0] = 0

                
                # For the sides
                Ex[t + 1, -dmp_l:, dmp_l:-dmp_l] = dmp_model(dmp_my,dmp_i) * Ex[t + 1, -dmp_l:, dmp_l:-dmp_l]
                Ey[t + 1, -dmp_l:, dmp_l:-dmp_l] = dmp_model(dmp_my,dmp_i) * Ey[t + 1, -dmp_l:, dmp_l:-dmp_l]
                Ez[t + 1, -dmp_l:, dmp_l:-dmp_l] = dmp_model(dmp_my,dmp_i) * Ez[t + 1, -dmp_l:, dmp_l:-dmp_l]

                Ex[t + 1, :-dmp_l, :dmp_l] = dmp_model(dmp_mx,dmp_i) * Ex[t + 1, :-dmp_l, :dmp_l]
                Ey[t + 1, :-dmp_l, :dmp_l] = dmp_model(dmp_mx,dmp_i) * Ey[t + 1, :-dmp_l, :dmp_l]
                Ez[t + 1, :-dmp_l, :dmp_l] = dmp_model(dmp_mx,dmp_i) * Ez[t + 1, :-dmp_l, :dmp_l]

                Ex[t + 1, :-dmp_l, -dmp_l:] = dmp_model(np.flip(dmp_mx,1),dmp_i) * Ex[t + 1, :-dmp_l, -dmp_l:]
                Ey[t + 1, :-dmp_l, -dmp_l:] = dmp_model(np.flip(dmp_mx,1),dmp_i) * Ey[t + 1, :-dmp_l, -dmp_l:]
                Ez[t + 1, :-dmp_l, -dmp_l:] = dmp_model(np.flip(dmp_mx,1),dmp_i) * Ez[t + 1, :-dmp_l, -dmp_l:]
                
                # For the corners
                Ex[t + 1, -dmp_l:, -dmp_l:] = dmp_model(dmp_cnr,dmp_i) * Ex[t + 1, -dmp_l:, -dmp_l:]
                Ey[t + 1, -dmp_l:, -dmp_l:] = dmp_model(dmp_cnr,dmp_i) * Ey[t + 1, -dmp_l:, -dmp_l:]
                Ez[t + 1, -dmp_l:, -dmp_l:] = dmp_model(dmp_cnr,dmp_i) * Ez[t + 1, -dmp_l:, -dmp_l:]

                Ex[t + 1, -dmp_l:, :dmp_l] = dmp_model(dmp_cnl,dmp_i) * Ex[t + 1, -dmp_l:, :dmp_l]
                Ey[t + 1, -dmp_l:, :dmp_l] = dmp_model(dmp_cnl,dmp_i) * Ey[t + 1, -dmp_l:, :dmp_l]
                Ez[t + 1, -dmp_l:, :dmp_l] = dmp_model(dmp_cnl,dmp_i) * Ez[t + 1, -dmp_l:, :dmp_l]


                #Check if the dampening is correct



                # Magnetic field
                Bx[t + 1, 1:-1, 1:-1] = Bx[t, 1:-1, 1:-1] - dt / dy * (Ez[t + 1, 1:-1, 2:] - Ez[t + 1, 1:-1, 1:-1])
                By[t + 1, 1:-1, 1:-1] = By[t, 1:-1, 1:-1] + dt / dx * (Ez[t + 1, 2:, 1:-1] - Ez[t + 1, 1:-1, 1:-1])
                Bz[t + 1, 1:-1, 1:-1] = Bz[t, 1:-1, 1:-1] - dt / dx * (
                            Ey[t + 1, 2:, 1:-1] - Ey[t + 1, 1:-1, 1:-1]) + dt / dy * (
                                                    Ex[t + 1, 1:-1, 2:] - Ex[t + 1, 1:-1, 1:-1])

                # Boundary conditions for the magnetic field

                # Dampen B-field at boundary
                # For the sides
                Bx[t + 1, -dmp_l:, dmp_l:-dmp_l] = dmp_model(dmp_my,dmp_i) * Bx[t + 1, -dmp_l:, dmp_l:-dmp_l]
                By[t + 1, -dmp_l:, dmp_l:-dmp_l] = dmp_model(dmp_my,dmp_i) * By[t + 1, -dmp_l:, dmp_l:-dmp_l]
                Bz[t + 1, -dmp_l:, dmp_l:-dmp_l] = dmp_model(dmp_my,dmp_i) * Bz[t + 1, -dmp_l:, dmp_l:-dmp_l]

                Bx[t + 1, :-dmp_l, :dmp_l] = dmp_model(dmp_mx,dmp_i) * Bx[t + 1, :-dmp_l, :dmp_l]
                By[t + 1, :-dmp_l, :dmp_l] = dmp_model(dmp_mx,dmp_i) * By[t + 1, :-dmp_l, :dmp_l]
                Bz[t + 1, :-dmp_l, :dmp_l] = dmp_model(dmp_mx,dmp_i) * Bz[t + 1, :-dmp_l, :dmp_l]

                Bx[t + 1, :-dmp_l, -dmp_l:] = dmp_model(np.flip(dmp_mx,1),dmp_i) * Bx[t + 1, :-dmp_l, -dmp_l:]
                By[t + 1, :-dmp_l, -dmp_l:] = dmp_model(np.flip(dmp_mx,1),dmp_i) * By[t + 1, :-dmp_l, -dmp_l:]
                Bz[t + 1, :-dmp_l, -dmp_l:] = dmp_model(np.flip(dmp_mx,1),dmp_i) * Bz[t + 1, :-dmp_l, -dmp_l:]
                
                # For the corners
                Bx[t + 1, -dmp_l:, -dmp_l:] = dmp_model(dmp_cnr,dmp_i) * Bx[t + 1, -dmp_l:, -dmp_l:]
                By[t + 1, -dmp_l:, -dmp_l:] = dmp_model(dmp_cnr,dmp_i) * By[t + 1, -dmp_l:, -dmp_l:]
                Bz[t + 1, -dmp_l:, -dmp_l:] = dmp_model(dmp_cnr,dmp_i) * Bz[t + 1, -dmp_l:, -dmp_l:]

                Bx[t + 1, -dmp_l:, :dmp_l] = dmp_model(dmp_cnl,dmp_i) * Bx[t + 1, -dmp_l:, :dmp_l]
                By[t + 1, -dmp_l:, :dmp_l] = dmp_model(dmp_cnl,dmp_i) * By[t + 1, -dmp_l:, :dmp_l]
                Bz[t + 1, -dmp_l:, :dmp_l] = dmp_model(dmp_cnl,dmp_i) * Bz[t + 1, -dmp_l:, :dmp_l]

                # #Reflecting boundary condition
                # Bx[t+1, 1:,-1] = 0
                # Bx[t+1,1:,1] = 0
                # Bx[t+1,-1,1:] = 0
                # By[t+1, 1:,-1] = 0
                # By[t+1,1:,1] = 0
                # By[t+1,-1,1:] = 0
                # Bz[t+1, 1:,-1] = 0
                # Bz[t+1,1:,1] = 0
                # Bz[t+1,-1,1:] = 0


                # New Current density
                Jx[t + 1, 1:-1, 1: -1] = (Jx[t, 1:-1, 1:-1] + dt * epsilon * Ope ** 2 * Ex[t + 1, 1:-1, 1:-1] 
                                            - (e_q * dt) / (e_m) * (((Jy[t, 1:-1, 1:-1] + Jy[t, 1: -1, :-2] + Jy[t,2:,1:-1] + Jy[t,2:,:-2]) * 2 * B0z) / 8 
                                                                                - ((Jz[t, 1:-1, 1:-1] + Jz[t,2:,1:-1]) * B0y / 2)))
                
                Jy[t + 1, 1:-1, 1:-1] = (Jy[t, 1:-1, 1:-1] + dt * epsilon * Ope ** 2 * Ey[t + 1, 1:-1, 1:-1] 
                                            - (e_q * dt) / (e_m) * (((Jz[t, 1:-1, 2:] + Jz[t, 1:-1,1:-1]) * B0x) / 2 
                                                                                - ((Jx[t, 1:-1,1:-1] + Jx[t, :-2,1:-1] + Jx[t,1:-1,2:] + Jx[t,:-2,2:]) * 2 * B0z) / 8))

                Jz[t + 1, 1:-1, 1:-1] = (Jz[t, 1:-1, 1:-1] + dt * epsilon * Ope ** 2 * Ez[t + 1, 1:-1, 1:-1] 
                                            - (e_q * dt) / (e_m) * (((Jx[t, 1:-1, 1:-1] + Jx[t, :-2,1:-1]) * 2 * B0y) / 4 
                                                                                - ((Jy[t, 1:-1, 1:-1] + Jy[t,1:-1,:-2]) * 2 * B0x) / 4))

                # Dampen J-field at boundary
                # For the sides
                Jx[t + 1, -dmp_l:, dmp_l:-dmp_l] = dmp_model(dmp_my,dmp_i) * Jx[t + 1, -dmp_l:, dmp_l:-dmp_l]
                Jy[t + 1, -dmp_l:, dmp_l:-dmp_l] = dmp_model(dmp_my,dmp_i) * Jy[t + 1, -dmp_l:, dmp_l:-dmp_l]
                Jz[t + 1, -dmp_l:, dmp_l:-dmp_l] = dmp_model(dmp_my,dmp_i) * Jz[t + 1, -dmp_l:, dmp_l:-dmp_l]

                Jx[t + 1, :-dmp_l, :dmp_l] = dmp_model(dmp_mx,dmp_i) * Jx[t + 1, :-dmp_l, :dmp_l]
                Jy[t + 1, :-dmp_l, :dmp_l] = dmp_model(dmp_mx,dmp_i) * Jy[t + 1, :-dmp_l, :dmp_l]
                Jz[t + 1, :-dmp_l, :dmp_l] = dmp_model(dmp_mx,dmp_i) * Jz[t + 1, :-dmp_l, :dmp_l]

                Jx[t + 1, :-dmp_l, -dmp_l:] = dmp_model(np.flip(dmp_mx,1),dmp_i) * Jx[t + 1, :-dmp_l, -dmp_l:]
                Jy[t + 1, :-dmp_l, -dmp_l:] = dmp_model(np.flip(dmp_mx,1),dmp_i) * Jy[t + 1, :-dmp_l, -dmp_l:]
                Jz[t + 1, :-dmp_l, -dmp_l:] = dmp_model(np.flip(dmp_mx,1),dmp_i) * Jz[t + 1, :-dmp_l, -dmp_l:]
                
                # For the corners
                Jx[t + 1, -dmp_l:, -dmp_l:] = dmp_model(dmp_cnr,dmp_i) * Jx[t + 1, -dmp_l:, -dmp_l:]
                Jy[t + 1, -dmp_l:, -dmp_l:] = dmp_model(dmp_cnr,dmp_i) * Jy[t + 1, -dmp_l:, -dmp_l:]
                Jz[t + 1, -dmp_l:, -dmp_l:] = dmp_model(dmp_cnr,dmp_i) * Jz[t + 1, -dmp_l:, -dmp_l:]

                Jx[t + 1, -dmp_l:, :dmp_l] = dmp_model(dmp_cnl,dmp_i) * Jx[t + 1, -dmp_l:, :dmp_l]
                Jy[t + 1, -dmp_l:, :dmp_l] = dmp_model(dmp_cnl,dmp_i) * Jy[t + 1, -dmp_l:, :dmp_l]
                Jz[t + 1, -dmp_l:, :dmp_l] = dmp_model(dmp_cnl,dmp_i) * Jz[t + 1, -dmp_l:, :dmp_l]

                # #Reflecting boundary condition
                # Jx[t+1, 1:,-1] = 0
                # Jx[t+1,1:,1] = 0
                # Jx[t+1,-1,1:] = 0
                # Jy[t+1, 1:,-1] = 0
                # Jy[t+1,1:,1] = 0
                # Jy[t+1,-1,1:] = 0
                # Jz[t+1, 1:,-1] = 0
                # Jz[t+1,1:,1] = 0
                # Jz[t+1,-1,1:] = 0

                if start0:
                    if field == 'Ez':
                        SaveMat = Ez[t,:,:]
                    elif field == 'Ey':
                        SaveMat = Ey[t,:,:]
                    start0 = False
                
            if field == 'Ez':
                SaveMat = np.dstack((SaveMat,Ez[intval,:,:]))
            if field == 'Ey':
                SaveMat = np.dstack((SaveMat,Ey[intval,:,:]))
                

            # Replaces the first row in the matrices, for use in the next interval
            Ex[0, :, :] = Ex[intval, :, :]
            Ey[0, :, :] = Ey[intval, :, :]
            Ez[0, :, :] = Ez[intval, :, :]
            Bx[0, :, :] = Bx[intval, :, :]
            By[0, :, :] = By[intval, :, :]
            Bz[0, :, :] = Bz[intval, :, :]
            Jx[0, :, :] = Jx[intval, :, :]
            Jy[0, :, :] = Jy[intval, :, :]
            Jz[0, :, :] = Jz[intval, :, :]

            #Predict computation time remaining:
            t1 = time.time()
            t_remain = (t1-t000) / ((k+1+kk*100)/T) - (t1-t000)
            t_remain_min = int(t_remain/60)
            t_remain_sec = int((t_remain/60 - t_remain_min) * 60)
            
            if t_remain_min >= 1:
                print(str(int(round((100*kk+k)/T*100,3))) + '% of simulation done. Est time remaining: ' + str(t_remain_min) + 'min ' + str(t_remain_sec) + 'sec ', end= '\r')
            else:
                print(str(int(round((100*kk+k)/T*100,3))) + '% of simulation done. Est time remaining: ' + str(t_remain_sec) + 'sec    ', end= '\r')


        timestamp = 't' + str(kk*100) + '-' + str(kk*100+99) 
        
        print('\n')

        #Saves the matrix to a text file:    
        print('t = '+ str(b+100*kk))
        print('k = '+ str(k+100*kk))
        t2 = time.time()
        print('Computation time for timestamp ' + timestamp + ': ' + str(t2-t000) + 's')
        print(SaveMat[:-1,:,:].shape)
        if not type(CustomName) == str:
            txtname = str(mapname + '_to_' + timestamp + '.npy')
        else:
            txtname = str(CustomName + '_to_' + timestamp + '.npy')

        np.save(txtname, SaveMat[:,:,:-1])

        
        print('Writing time: '+ str(time.time()-t2) + 's')

        SaveMat = SaveMat[:,:,-1]

    os.chdir(fallback)