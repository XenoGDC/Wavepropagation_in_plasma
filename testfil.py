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

    mat1 = np.load(os.path.join(path0, 'ne_1.npy'))
    x_list = np.load(os.path.join(path0, 'x_list.npy'))
    y_list = np.load(os.path.join(path0, 'y_list.npy'))

    I, J = np.shape(mat1)

    x_list = x_list[:I]
    y_list = y_list[:J]

    sigma = 0.02/rd.dx

    I += 1
    J += 1

    # Write the directory for where you want to store your files here. It should contain the folder which holds the plasma profiles.
    data_collection_location = 'C:/Users/augus/Onedrive/skrivebord/DTU/Fagprojekt (bølgeudbredelse i plasme)'
    os.chdir(data_collection_location)

    # Name of folder to hold the densities and the simulated data
    directory_name = 'Real_blobs_over_time'
    os.chdir(directory_name)
    files = os.listdir()
    angle_X = np.pi/2

    # name_vacuum = 'CustomP_vacuum_comparison'
    # sim.WaveSim(200,I,J,sigma,field='Ez',B0=[0,0,0.5],wave_polarity_angle=0,CustomName=name_vacuum)

    # Make list over plasma profiles
    blobs = [file for file in files if file.endswith('.npy')]
    blobs = sorted(blobs)

    # Plot the plasma profiles
    if False:
        for i in range(len(blobs)):
            blob_matrix = np.load(blobs[i])
            plt.pcolormesh(y_list, x_list[:I-1], blob_matrix)
            plt.show()

    # Simulate all blobs
    if False:
        print('Simulating for the blobs')
        B0 = [0, 0, 0.5]
        for blob in blobs:
            # Clear the screen
            os.system('cls')
            filename = '{blob}_simulation'.format(blob=blob[:-4])
            p_matrix = np.load(blob)
            I, J = np.shape(p_matrix)
            I += 1
            J += 1

            if True:
                print('Simulating X-mode for {blob}'.format(blob=blob))
                filenameX = '{default}_X-mode'.format(default=filename)
                # X-mode first
                sim.WaveSim(200, I, J, sigma, B0=B0,
                            CustomName=filenameX, CustomPMatrix=p_matrix)

            if False:
                print('Simulating O-mode for {blob}'.format(blob=blob))
                filenameO = '{default}_O-mode'.format(default=filename)
                # O-mode after
                sim.WaveSim(200, I, J, sigma, B0=B0, Linear_angle=angle_X,
                            field='Ey', CustomName=filenameO, CustomPMatrix=p_matrix)

    # Make sorted lists of simulations
    data_sim = [file for file in files if not '.' in file]
    Xmodes = [file for file in data_sim if file.endswith('X-mode')]
    Omodes = [file for file in data_sim if file.endswith('O-mode')]
    Xmodes = sorted(Xmodes)
    Omodes = sorted(Omodes)

    # If there should be videos
    if False:
        if True:
            # Making videos for the X-mode simulations
            for Xmode, blob in zip(Xmodes[1:], blobs[1:]):
                # for i in range(len(Xmodes)):
                p_matrix = np.load(blob)
                mkr.BigMovieMaker(Xmode, Pmatrix=p_matrix)

        if False:
            # Making videos for the O-mode simulations
            for Omode, blob in zip(Omodes, blobs):
                # for i in range(len(Omodes)):
                p_matrix = np.load(blob)
                mkr.BigMovieMaker(Omode, Pmatrix=p_matrix)

    # Wave width over time
    if False:
        time_diff = 1.46e-6
        print('Plotting induced widths over time:')
        position = int(0.06/rd.dy)
        Widths = np.array([0.04555448373672055])
        dWidths = np.array([1.4987636159051634e-05])
        times = np.array([0])
        t = time_diff*10**6
        for file, blob in zip(Xmodes, blobs):
            matrix, hundred = rd.ReadBigSim(file, 1)

            sig, amp, dsig, damp = blb.GaussSingleAnalysis(
                matrix, 50, position=position, sig0=sigma)
            Widths = np.append(Widths, sig*2*rd.dx)
            dWidths = np.append(dWidths, dsig*2*rd.dx)
            times = np.append(times, t)
            t += time_diff*10**6

        # Plot interesting plasma profiles
        if True:
            blob1 = np.load(blobs[0])
            t1 = time_diff*1e6
            blob2 = np.load(blobs[3])
            t2 = time_diff*4e6
            matrix_ex, hundred = rd.ReadBigSim(Xmodes[3], 1)
            mat_capt = matrix_ex[:, :, 50]

            blob3 = np.load(blobs[7])
            t3 = time_diff*8e6
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

            x_list = np.arange(I)*rd.dx
            y_list = np.arange(J)*rd.dy

            ymarkline = np.array([position*rd.dx, position*rd.dx])
            xmarkline = np.array([0, (J-1)*rd.dx])

            fig.suptitle(
                'Widths at the plasma, 6cm away from start width of 4cm')

            ax1.set_title('Measured widths')
            ax1.scatter(times, Widths)
            ax1.errorbar(times, Widths, dWidths)
            ax1.set_xlabel('Timepoint [$\mu$s]')
            ax1.set_ylabel('Width at plasma [m]')

            ax2.pcolormesh(y_list, x_list, blob1)
            ax2.set_title('t = {t1} [$\mu$s]'.format(t1=t1))
            ax2.plot(xmarkline, ymarkline, label='Measured point')
            ax2.legend()
            ax2.set_xlabel('y [m]')
            ax2.set_ylabel('x [m]')

            ax3.pcolormesh(y_list[:], x_list[:],
                           mat_capt[:-1, :-1], cmap='seismic')
            ax3.contour(y_list[:-1], x_list[:-1], blob2)
            ax3.set_title('t = {t2} [$\mu$s]'.format(t2=t2))
            ax3.set_xlabel('y [m]')
            ax3.set_ylabel('x [m]')

            cmap = ax4.pcolormesh(y_list, x_list, blob3)
            cbar = fig.colorbar(cmap)
            cbar.set_label('Plasma density')
            ax4.set_title('t = {t3} [$\mu$s]'.format(t3=t3))
            ax4.set_xlabel('y [m]')
            ax4.set_ylabel('x [m]')

            fig.tight_layout(pad=1.06)

            fig.savefig('Beam_widths_over_time_with_blobs.png')
            fig.show()
            input()

        if False:
            plt.figure()

            plt.scatter(times, Widths)
            plt.errorbar(times, Widths, dWidths)
            plt.title('Widths at the plasma, 6cm away from start width of 4cm')
            plt.xlabel('Timepoint [$\mu$s]')
            plt.ylabel('Width at plasma [m]')
            plt.savefig('Beam_widths_over_time.png')
            plt.show()

        # State the increase in width
        increase_width = np.max(Widths)/np.min(Widths)-1
        print('The beam gets {inc}% wider at the worst case'.format(
            inc=round(100*increase_width, 3)))
        print(Widths)

    # Plot Gaussian fitting
    if False:
        os.chdir(os.path.join(data_collection_location, 'Misc data'))
        vac_mat = rd.SimReader('Ez_vacuum_dmpl50_dmp_i_7.npy')
        blb.GaussSingleAnalysis(vac_mat, 150, 500, plotcurve=True)

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

    # Plot the density for the linear increase
    if False:
        d_matrix = rd.gimidensity(2400, 2400, 0, 'Linear', 2400*3/4)
        rd.plotplasmadens(2400, 2400, d_matrix, 'Linear', 'Both',
                          B0=0.25, cutoffp=int(2400*3/4), cmap='plasma')

    # Plot different captures of the fields:
    if False:
        # for the blob:
        if False:

            os.chdir(os.path.join(
                data_collection_location, 'Blob_n_5e18_width_var'))
            d_matrix = rd.gimidensity(2400, 2400, 0, 'Blob', x0=int(
                2400/3), y0=1200, signy=41, peak=5e18)
            rd.plotplasmadens(2400, 2400, d_matrix, 'Blob', cmap='plasma')

            mat, hundred = rd.ReadBigSim('X-mode_blob_test_sig_5e+1841_n_', 2)

            rd.PlotTimePoint(mat, 'Blob', 40, hundred, 'Ey')

        # For the vacuum
        if False:
            os.chdir(data_collection_location)
            os.chdir('Misc data')
            mat = rd.SimReader('Ez_vacuum_dmpl50_dmp_i_7')
            rd.PlotTimePoint(mat, 'Vacuum', 150, 0, 'Ez')

        # For the linear
        if False:
            print('Plotting the linear simulation')
            os.chdir(data_collection_location)
            mat, hundred = rd.ReadBigSim('Ez_Linear_V_big', 5)
            rd.PlotTimePoint(mat, 'Linear', 50, hundred, 'Ez')

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

    rd.PlotDampenedBoundary(100, 100, 50, 7)


except:
    import traceback
    traceback.print_exc()
    input()
