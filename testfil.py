import SimMoviemaker as mkr
import numpy as np
import os
import Wave_Sim_reader as rd
import WavePropagationSim as sim
from matplotlib import pyplot as plt
import BlobDensityDispersionAnalyser as blb
from playsound import playsound


try:
    # I,J = 600,600
    dfont = 'DejaVu Sans'

    # Get the current path
    path0 = os.getcwd()


    # Matrix,hundred = rd.ReadBigSim(data_name)

    # I,J,T = np.shape(Matrix)

    # pMatrix = rd.gimidensity(I,J,0,'Blob',x0=int(I/3),y0=int(J/2),signy=41,peak=5e8)

    # rd.plotplasmadens(I,J,pMatrix,'Blob','X-mode',B0=0.5)

    mat1 = np.load(os.path.join(path0,'ne_1_scaled.npy'))
    # mat2 = np.load(os.path.join(path0,'ne_2_scaled.npy'))
    # mat3 = np.load(os.path.join(path0,'ne_3_scaled.npy'))
    x_list = np.load(os.path.join(path0,'x_list.npy'))
    y_list = np.load(os.path.join(path0,'y_list.npy'))

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
    

    # name_vacuum = 'CustomP_vacuum_comparison'
    # sim.WaveSim(200,I,J,sigma,field='Ez',B0=[0,0,0.5],wave_polarity_angle=0,CustomName=name_vacuum)
    
    # Name of folder to hold the densities and the simulated data
    directory_name = 'Real_blobs_over_time'
    os.chdir(directory_name)
    files = os.listdir()

    # Make list over plasma profiles
    blobs = [file for file in files if file.endswith('.npy')]
    blobs = sorted(blobs)
    
    # Plot the plasma profiles
    if True:
        for i in range(len(blobs)):
            blob_matrix = np.load(blobs[i])
            plt.pcolormesh(y_list,x_list[:I-1],blob_matrix)
            plt.show()

    # Simulate all blobs
    if True:
        print('Simulating for the blobs')
        B0 = [0,0,0.5]
        for blob in blobs:
            # Clear the screen
            os.system('cls')
            filename = '{blob}_simulation'.format(blob=blob[:-4])
            p_matrix = np.load(blob)
            I,J = np.shape(p_matrix)
            I += 1
            J += 1

            if True:
                print('Simulating X-mode for {blob}'.format(blob=blob))
                filenameX = '{default}_X-mode'.format(default=filename)
                # X-mode first
                sim.WaveSim(200,I,J,sigma,B0=B0,CustomName=filenameX,CustomPMatrix=p_matrix)

            if False:
                print('Simulating O-mode for {blob}'.format(blob=blob))
                filenameO = '{default}_O-mode'.format(default=filename)
                # O-mode after
                sim.WaveSim(200,I,J,sigma,B0=B0,Linear_angle=angle_X,field='Ey',CustomName=filenameO,CustomPMatrix=p_matrix)
            
            
    # Make sorted lists of simulations
    data_sim = [file for file in files if not '.' in file]
    Xmodes = [file for file in data_sim if file.endswith('X-mode')]
    Omodes = [file for file in data_sim if file.endswith('O-mode')]
    Xmodes = sorted(Xmodes)
    Omodes = sorted(Omodes)

    # If there should be videos
    if False:
        # Making videos for the X-mode simulations
        for i in range(len(Xmodes)):
            mkr.BigMovieMaker(Xmodes[i],Pmatrix=blob[i])

        # Making videos for the O-mode simulations
        for i in range(len(Omodes)):
            mkr.BigMovieMaker(Omodes[i],Pmatrix=blob[i])

    # Wave width over time
    if True:
        print('Plotting induced widths over time:')
        position = int(0.06/rd.dy)
        Widths = np.array([0.04555448373672055])
        dWidths = np.array([1.4987636159051634e-05])
        times = np.array([0])
        t = 0.7
        for file, blob in zip(Xmodes,blobs):
            matrix,hundred = rd.ReadBigSim(file,1)

            sig,amp,dsig,damp = blb.GaussSingleAnalysis(matrix,50,position=position,sig0=sigma)
            Widths = np.append(Widths,sig*2*rd.dx)
            dWidths = np.append(dWidths,dsig*2*rd.dx)
            times = np.append(times,t)
            t += 0.7
        
        plt.figure()
        
        plt.scatter(times,Widths)
        plt.errorbar(times,Widths,dWidths)
        plt.title('Widths at the plasma, 6cm away from start width of 4cm')
        plt.xlabel('Timepoint [$\mu$s]')
        plt.ylabel('Width at plasma [m]')
        plt.savefig('Beam_widths_over_time.png')
        plt.show()


    # name1x = 'Sim_for_ne_1_X'
    # sim.WaveSim(300,I,J,sigma,field='Ey',CustomPMatrix=mat,B0=[0,0,0.5],wave_polarity_angle=0,CustomName=name)

    # # newmat,hundred = rd.ReadBigSim(name)

    
    # mat2,hundred = rd.ReadBigSim(name3o,1)
    # rd.GaussAnalyser(mat2,50,sig0=sigma,Vacuum=False)
    # mat1 = rd.SimReader('Misc Data/Ez_vacuum_dmpl50_dmp_i_7')
    # rd.GaussAnalyser(mat1,150)

    # mat11, hun1 = rd.ReadBigSim(name_vacuum,1)
    # mat11, hun1 = rd.ReadBigSim(name1o,1)
    # mat21, hun2 = rd.ReadBigSim(name2o,1)
    # mat31, hun3 = rd.ReadBigSim(name3o,1)

    # name1 = 'Wave_1_O_capt field plot'  
    # name2 = 'Wave_2_O_capt field plot'
    # name2 = 'Wave_3_O_capt field plot'

    # # rd.PlotTimePoint(mat11,name1,60,hun1,field='Ez',Pmatrix=mat1)
    # # rd.PlotTimePoint(mat21,name2,60,hun2,field='Ez',Pmatrix=mat2)
    # # rd.PlotTimePoint(mat31,name2,60,hun2,field='Ez',Pmatrix=mat3)

    # # d_matrix = rd.gimidensity(2400,2400,0,'Linear',2400*3/4)
    # # rd.plotplasmadens(2400,2400,d_matrix,'Linear','Both',B0=0.25,cutoffp=int(2400*3/4))

    # blb.Blobdispersion('Blob_n_5e18_width_var',250,plotname = 'Blob_Dispersion_sig_var_n5e18',plottitle = 'Blob Dispersion with varying blob width'.format(n=round(5e18,3)),varmode = 'Width')
    # blb.Blobdispersion('Blob_n_1e18-1e19',250,plotname = 'Blob_Dispersion_n_var_sig30',plottitle = 'Blob Dispersion with varying plasma density'.format(n=round(5e18,3)))
    
    
    # position = int(0.06/rd.dy)
    # amp1,sig1,damp1,dsig1 = blb.GaussSingleAnalysis(mat11,50,position=position,sig0 = sigma)
    # amp2,sig2,damp2,dsig2 = blb.GaussSingleAnalysis(mat21,50,position=position,sig0 = sigma)
    # amp3,sig3,damp3,dsig3 = blb.GaussSingleAnalysis(mat31,50,position=position,sig0 = sigma)

    # sigs = np.array([2*sig1,2*sig2,2*sig3])*rd.dy
    # dsigs = np.array([2*dsig1,2*dsig2,2*dsig3])*rd.dy

    # print('The widths are {sigs} $\plusmn$ {dsigs}'.format(sigs=sig1*2*rd.dy,dsigs=dsig1*2*rd.dy))
    # input()
    



except:
    import traceback
    traceback.print_exc()
    input()