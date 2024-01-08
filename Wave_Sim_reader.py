import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import WavePropagationSim as sim


# Constants
c = 2.99792458e8
mu = 4e-7 * np.pi
epsilon = 8.8541878e-12 # vacuum permebility
e_m = 9.1093837015e-31
e_q = 1.6022e-19

O = 2 * np.pi * 28e9
#x,y, and t is dx, dy, and dt
dx = np.pi * c / (32 * O)
dy = np.pi * c / (32 * O)
dt = dx / c / 2

lambdt = 2*np.pi/O
lam = 2*np.pi*c/O
sig0 = 60
sig = sig0*dx
path0 = os.getcwd()

# class Sim():
#     def __init__(self) -> None:
#         pass

def SimReader(Customname: str = None):
    if Customname == None:

        files = os.listdir(path0)
        data_files = [file for file in files
                      if file.endswith('.npy')]
        print('Select file with the number')
        for i in range(len(data_files)):
            print(str(i) + '. ' + str(data_files[i]))
        
        Input = True
        while Input:
            try:
                type0 = input()
                chosen = str(data_files[int(type0)])
                arr = np.load(chosen)
                Input = False
            except:
                print('Please choose a valid file')
    else:
        
        if Customname.endswith('.npy'):
            arr = np.load(Customname)
        else:
            arr = np.load(Customname + '.npy')
    
    return arr
        
def ReadBigSim(filename: str = None,hundred: int = None):
    fallback = os.getcwd()
    os.chdir(filename)

    
    filelist = os.listdir()
    datalist = [dat for dat in filelist if dat.endswith('.npy')]

    if hundred == None:
        for i in range(len(datalist)):
            print(str(i) + '. ' + datalist[i])

        choice = int(input('Choose file: '))
    else:
        choice = hundred
    choicedat = datalist[int(choice)]

    Mat = SimReader(Customname=choicedat)
    os.chdir(fallback)
    return Mat,choice    

def GaussAnalyser(Matrix, TimePosition: int, ShowStar = False,sig0 = 60):
    Mat = Matrix
    I,J,T = np.shape(Mat)
    print(I,J,T)

    sig = sig0*dx
    wait = True
    ad = 0
    X = len(Mat[:,0,0])
    Y = len(Mat[0,:,0])
    # w_count = int(input('How many slices?\n'))
    w_count = 100

    #Model for the gaussian curve
    def model_f(x,a,b):
        return b*np.exp(-((x - Y / 2)**2 / a**2))
    
    # while wait:
    

    # Defines X-positions where we analyse the curve along the y-axis
    sig_pos = np.array([1+ad])
    for ii in range(1,w_count):
        sig_pos = np.append(sig_pos,int(ii*X/w_count+ad))
    

    Matcapt = np.matrix(Mat[:,:,TimePosition])
                
    SigMatBig = np.abs(np.array(Mat[sig_pos[0],:,TimePosition]))
            
    for ii in range(1,w_count):
        absMat = np.abs(Mat[sig_pos[ii],:,TimePosition])
        SigMatBig = np.vstack((SigMatBig,absMat))
        
    for tt in range(TimePosition+1,int(TimePosition+2*lambdt/dt/10)):
        SigMatBigTemp = np.abs(np.array(Mat[sig_pos[0],:,tt]))
            
        for ii in range(1,w_count):
            SigMatBigTemp = np.vstack((SigMatBigTemp,np.abs(Mat[sig_pos[ii],:,tt])))
        SigMatBig = np.dstack((SigMatBig,SigMatBigTemp))
    

    #Takes the averages
    n_maybe = len(SigMatBig[0,0,:])
    SigMatBig_sum = np.zeros(SigMatBigTemp.shape)
    for ii in range(n_maybe):
        SigMatBig_sum += SigMatBig[:,:,ii]

    SigMat = (SigMatBig_sum/n_maybe)

    x_ting = np.linspace(0,Y,Y)
    
    popt, pcov = curve_fit(model_f,x_ting,SigMat[0,:],p0=[sig0,1])
    perr = np.sqrt(np.diag(pcov))[0]
    if ShowStar:
        plt.title(str('Y-aksen plottet på x = ' + str(sig_pos[0])))
        # plt.plot(model_f(x_ting,popt))
        plt.plot(x_ting,SigMat[0,:])
        plt.draw()
        plt.pause(0.01)
        plt.clf()


    
    for ii in range(1,w_count):
        #Makes a new popt and pcov we call poptx and pcovx for each sample.
        poptx, pcovx = curve_fit(model_f,x_ting,SigMat[ii,:],p0=[np.max(popt),1],maxfev=5000)
        #Stacks the individual popt and pcov on top of each other, even though the pcov isn't going to be used because the dimensions do not fit.
        popt, pcov = np.vstack((popt,poptx)), np.vstack((pcov,pcovx))
        
        #Makes an individual error for the sigma to append into a list of the errors
        perrx = np.sqrt(np.diag(pcovx))[0]
        perr = np.append(perr,perrx)
        
        #Plots each fittet gaussian curve to give an idea of what is going on during the process
        # plt.plot(model_f(x_ting,popt[ii,0],popt[ii,1]))
        if ShowStar:
            plt.plot(x_ting,SigMat[ii,:],'r')
            plt.title(str('Y-aksen plottet på x = ' + str(sig_pos[ii])))
            plt.draw()
            plt.pause(0.01)
            plt.clf()




    #Plots the final graphs
    #____________________________________
    fig, (ax1,ax2) = plt.subplots(1,2)
    zr = np.pi*sig**2/lam
    linx = np.linspace(0,X+1)
    liny = sig * np.sqrt(1+(linx*dx/zr)**2)


    #In the errorbar we scale the value of sigma and it's error to the real size by multiplying it with dy
    perr = dy*perr
    ax2.scatter(sig_pos*dx,dy*popt[:,0],label='Simulation')
    ax2.plot(linx*dx, liny,label='Theory')
    ax2.errorbar(sig_pos*dx,dy*popt[:,0],yerr=perr,fmt='o')
    ax2.set_title('Gaussian beam dispersion')
    ax2.set_xlabel('Position [m]')
    ax2.set_ylabel('$\sigma$ [m]')
    ax2.legend()
    leg = ax2.get_legend()
    leg.legendHandles[0].set_color('orange')
    leg.legendHandles[1].set_color('blue')
    
    #____________________________________
    
    
    #On the left subplot the electric field in the Ez direction is plottet at the time the fitting takes place and at the lines we fit the gaussian curve to.
    # for ii in range(w_count):
    #     Mat[int(sig_pos[ii]), :,TimePosition] = 100

    xlist = np.linspace(0,I,I)*dx
    ylist = np.linspace(-J/2,J/2,J)*dy

    map1 = ax1.pcolormesh(ylist,xlist,Matcapt, cmap='seismic',vmin=-1, vmax=1)
    cbar = plt.colorbar(map1)
    cbar.set_label('Ez')
    ax1.set_title('Ez-field map')
    ax1.set_xlabel('y [m]')
    ax1.set_ylabel('x [m]')
    


    fig.tight_layout()
    fig.savefig('Sigma_val_testplot.png',format='png')
    # plt.show()
    plt.draw()

    plt.pause(0.01)

    print(liny[0])

    

    pass

def PlotTimePoint(Matrix,name,t:int,hundred:int,field):
    I,J,T = np.shape(Matrix)
    linex = np.arange(I) * dx
    liney = np.arange(J) * dy
    timepoint = hundred*100 + t
    
    plt.figure()
    map1 = plt.pcolormesh(linex,liney,Matrix[:,:,t], cmap='seismic',vmin=-1,vmax=1)
    cbar = plt.colorbar(map1)
    cbar.set_label('{field} intensity'.format(field=field))
    plt.title('{name} at t = {time}'.format(name = name,time=timepoint))
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.savefig('{name} field plot.png'.format(name=name))
    plt.show()

def Waveabsavg(Matrix,start):
    difference = 2*int(lambdt/dt/11)
    print(difference)
    MatAvg = np.zeros(Matrix[:,:,0].shape)
    for tt in range(start,start+difference):
        MatAvg += np.abs(Matrix[:,:,tt])
    
    MatAvg = MatAvg/difference
    # np.sum(np.abs(Matrix[:,:,start:start+int(lambdt/dt)]),keepdims=True)/(start+int(lambdt/dt))
    
    map1 = plt.pcolormesh(MatAvg, cmap='seismic',vmin=-1, vmax=1)
    cbar = plt.colorbar(map1)
    cbar.set_label('Ez intensity')
    # plt.set_title('Ez-field as average')
    # plt.set_xlabel('y')
    # plt.set_ylabel('x')

def CutoffAnalyser(Matrix,TimePosition: int(),choicedx: float = dx
                   ,customname: str = ''):
    Mat = np.abs(Matrix)
    I = int((len(Mat[:,0,0])-1)/2)
    

    Midline = Mat[:,I,TimePosition]
    Top = []
    Pos = []
    
    for i in range(1,len(Midline)-1):
        Before = Midline[i] - Midline[i-1]
        After = Midline[i] - Midline[i+1]
        if Before > 0 and After > 0:
            Top.append(Midline[i])
            Pos.append(i)
    
    WaveLength = []
    for i in range(1,len(Pos)):
        Wave = (Pos[i]-Pos[i-1])*2*choicedx
        WaveLength.append(Wave)
    
    # plt.scatter(np.array(Pos)[1:]*choicedx,WaveLength)
    plt.xlabel('Position')
    # plt.ylabel('Wavelength')
    plt.title(customname)
    # plt.savefig(customname + '.png')
    # plt.close()

    plt.scatter(Pos,Top)
    plt.show()

def CMAanalysis(Matrix,density_matrix,cutoff:int,B0,
                TimePosition: int = None,choicedx: float = dx,hundred: int = 0,
                field: str = 'E',omega: float = O,
                mode: str = 'O',
                title:str = 'CMA_testplot.png'):
    # Taking the absolute to the matrix to sort out the waves
    abMat = np.abs(Matrix)

    if len(Matrix[0,0,:]) > 100:
        TimePosition += hundred*100
        hundred = 0




    I = int((len(abMat[0,:,0])-1)/2)
    II,J,T = np.shape(Matrix)

    # Calculate omega_ce and omega_pe
    omega_ce = B0/e_m*np.abs(e_q)
    OPemidline = (density_matrix[:,I]*e_q**2/(e_m*epsilon))**0.5



    Midline = abMat[:,I,TimePosition]
    
    Top = []
    Pos = []
    Ope_point = []
    
    for i in range(1,len(Midline)-1):
        Before = Midline[i] - Midline[i-1]
        After = Midline[i] - Midline[i+1]
        if Before > 0 and After > 0:
            Top.append(Midline[i])
            Pos.append(i)
            Ope_point.append(OPemidline[i])

    
    WaveLength = []
    for i in range(1,len(Pos)):
        # We're taking half wavelengths when we subtract the top points of the absolute field, 
        # so we should probably time this difference by 2 to get the full wavelength right?
        Wave = (Pos[i]-Pos[i-1])*choicedx*2
        WaveLength.append(Wave)


    S_Top = []
    S_Pos = []
    S_Ope = []
    S_wavelength = []

    # S_wavelength.append(WaveLength[0])

    for t in range(1,len(Pos[1:])):
        if Top[t] > 0.3:
            S_Top.append(Top[t])
            S_Pos.append(Pos[t])
            S_Ope.append(Ope_point[t])
            S_wavelength.append(WaveLength[t-1])
        else:
            break
    
    S_k = 2*np.pi/np.array(S_wavelength)

    Oce_pr_O = omega_ce/omega

    if mode == 'O':
        k_f = (omega**2/c**2-np.array(OPemidline)**2/c**2)
        cutoffpoint = int(cutoff)
    elif mode == 'X':
        k_f = ((omega**2-omega*omega_ce-OPemidline**2)*(omega**2+omega*omega_ce-OPemidline**2)/
        (c**2*(omega**2-omega_ce**2-OPemidline**2)))

        cutoffpoint = int(cutoff*(1-Oce_pr_O/np.tan(np.pi/4)))

    
    
    omega2 = c*S_k

    Ope_pr_O_sq = (np.array(S_Ope)/omega)**2
    


    # Ting bliver meget store og jeg ved ikke helt hvorfor
    lindturn = np.linspace(0,1,10)
    Ilindturn = np.flip(lindturn)
    

    #Plots the final graphs
    #____________________________________

    xlist = np.linspace(-II/2,II/2,II)*dx
    ylist = np.linspace(0,J,J)*dy


    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    
    fig.set_figwidth(12*0.8)
    fig.set_figheight(9*0.8)
    fig.set_tight_layout(tight=True)
    
    ax2.scatter(Ope_pr_O_sq,Oce_pr_O*np.ones(len(Ope_pr_O_sq)))
    ax2.plot(lindturn,Ilindturn)

    # For S = 0
    x_ax = np.linspace(0,1,100)
    y_ax = (1-x_ax)**0.5
    ax2.plot(x_ax,y_ax)

    ax2.set_title('CMA-Diagram')
    ax2.set_xlabel('$\omega_{pe}^2/\omega^2$')
    ax2.set_ylabel('$\omega_{ce}/\omega$')
    #____________________________________
    
    
    #On the left subplot the electric field in the Ez direction is plottet at the time the fitting takes place and at the lines we fit the gaussian curve to.
    # for ii in range(w_count):
    #     Mat[int(sig_pos[ii]), :,TimePosition] = 100
    MatWithCutoff = Matrix[:,:,TimePosition]
    if len(MatWithCutoff[:,0]) >= 800:
        MatWithCutoff[cutoffpoint-4:cutoffpoint+4,:] = 1000
    else:
        MatWithCutoff[cutoffpoint,:] = 1000
    map1 = ax1.pcolormesh(xlist,ylist,MatWithCutoff, cmap='seismic',vmin=-1, vmax=1)
    cbar = plt.colorbar(map1)
    cbar.set_label(field + ' intensity')
    ax1.set_title(field + '-field at t = ' + str(np.round((TimePosition+100*hundred)*dt*10**9,3)/10**9) + 's')
    ax1.set_xlabel('y [m]')
    ax1.set_ylabel('x [m]')

    ax3.plot(np.linspace(0,len(Midline),len(Midline))*dx,Midline)
    ax3.set_title('Absolute Wave')
    ax3.set_xlabel('x [m]')
    ax3.set_ylabel('Amplitude')

    ax4.plot(np.array(OPemidline)**2,k_f)
    ax4.set_title('Dispersion on ' + mode + '-mode')
    ax4.set_xlabel('$\omega_{pe}^2$')
    ax4.set_ylabel('$k^2$',)
    ax4.set_ylim(-1e6,1e6)
    ax4.scatter(np.array(S_Ope)**2,S_k**2)
    # ax4.plot(S_k**2,np.array(S_Ope)**2)
    # # ax4.set_title('Wavelengths')
    # ax4.set_xlabel('k^2')
    # ax4.set_ylabel('Ope^2')
    

    fig.tight_layout()
    fig.savefig(title,format='png')
    plt.draw()


    plt.show()

def gimidensity(I,J,Linear_angle,cutoff: int = 0):
    scale = np.zeros([I-1,J-1])
    lintheta = Linear_angle
    if cutoff == 0:
        cutoff = I/2
    def lindensityincrease(lintheta,x0,y0,n0,matrix):
        n = n0/(x0*np.cos(lintheta)+y0*np.sin(lintheta))
        scaletemp = np.ones(matrix.shape)
        for xx in range(len(matrix[:,0])):
            for yy in range(len(matrix[0,:])):
                scaletemp[xx,yy] = n*(xx*np.cos(lintheta)+yy*np.sin(lintheta))
        return scaletemp

    scale += lindensityincrease(lintheta,cutoff,J/2,0.97e19,scale)
    return scale

def plotplasmadens(I,J:int,cutoffp,densitymatrix,B0:float,mode:str):
    fig, ax = plt.subplots()
    linex = np.arange(I)*dx
    liney = np.arange(J)*dy
    x1 = np.array([0,J-1])
    xmode = 1-B0/np.tan(np.pi/4)
    if mode == 'O':
        cutoff = cutoffp
    elif mode == 'X':
        cutoff = cutoffp * xmode
    y1 = np.array([dx*cutoff,dx*cutoff])
    ax.plot(x1*dx, y1)
    map1 = ax.pcolormesh(linex,liney,densitymatrix, cmap='seismic')
    cbar = plt.colorbar(map1)
    cbar.set_label('Plasmadensity')
    ax.set_title('Plasmadensity')
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    fig.savefig('densityfig' + mode + '.png')
    fig.show()

def AnalyseBlobDispersion(Matrix,Wavepoint,plotname:str='Blob_plot.png',timepoint:int=50):
    I,J,T = np.shape(Matrix)

    liney = np.arange(J)

    Abs_Matrix = np.abs(Matrix)

    Wavelength_time = lambdt/dt/10

    #Model for the expected 2 gaussian curves
    def model_b(x,a,b,lx0,c,d,rx0):
        left = b*np.exp(-((x - lx0)**2 / a**2))
        right = d*np.exp(-((x - rx0)**2 / c**2))
        return left + right
    
    Matrix_sum = Abs_Matrix[:,:,timepoint]
    timerange = np.linspace(timepoint+1,timepoint+int(Wavelength_time),int(Wavelength_time-1))

    
    if False:
        
        # For the bunches of wavewidths

        Area = np.linspace(I/2,I*4/5,int(I/24))

        Width = []

        guess_lx0,guess_rx0 = J/4,J*3/4
        guess_amp_l,guess_amp_r = 0.1,0.1
        guess_sig0l,guess_sig0r = 30,30
        for y in Area:
            try:
                for t in timerange:
                    Matrix_sum += Abs_Matrix[:,:,int(round(t))]
                Matrix_avg = Matrix_sum/Wavelength_time

                # line = Matrix_avg[int(I*2/3),:]
                line = Matrix_avg[int(y),:]

                # plt.figure()
                # plt.plot(liney,line)
                # plt.show()
                
                # x_ting = np.linspace(0,Y,Y)

                
                popt, pcov = curve_fit(model_b,liney,line,p0=[guess_sig0l,guess_amp_l,guess_lx0
                                                            ,guess_sig0r,guess_amp_r,guess_rx0],maxfev=5000)
                perr = np.sqrt(np.diag(pcov))[0]

                aa,bb,lx00,cc,dd,rx00 = popt
                guess_lx0,guess_rx0,guess_amp_l,guess_amp_r,guess_sig0l,guess_sig0r = popt

                Regressed_function = model_b(liney,aa,bb,lx00,cc,dd,rx00)

                # plt.figure()
                # plt.plot(liney,Regressed_function)
                # plt.show()

                Split_width = np.abs(lx00-rx00) + aa + cc

                Width.append(Split_width)
            except:
                Width.append(1)

        # print(Split_width)
        # print(Split_width*dy)

        plt.figure()
        plt.plot(Area*dx,np.array(Width)*dy)
        plt.show()



    # For single points:
    # ____________________________________________________________
    if True:
        print('Plotting for x = {timepoint}'.format(timepoint=timepoint))
        guess_lx0,guess_rx0 = J/4,J*3/4
        guess_amp_l,guess_amp_r = 0.1,0.1
        guess_sig0l,guess_sig0r = 50,50

        for t in timerange:
            Matrix_sum += Abs_Matrix[:,:,int(round(t))]
        Matrix_avg = Matrix_sum/Wavelength_time

        line = Matrix_avg[int(Wavepoint),:]

        plt.figure()
        plt.plot(liney,line)
        plt.show()
        
        popt, pcov = curve_fit(model_b,liney,line,p0=[guess_sig0l,guess_amp_l,guess_lx0
                                                    ,guess_sig0r,guess_amp_r,guess_rx0],maxfev=5000)
        perr = np.sqrt((pcov))[0]

        aa,bb,lx00,cc,dd,rx00 = popt
        
        Regressed_function = model_b(liney,aa,bb,lx00,cc,dd,rx00)

        plt.figure()
        plt.plot(liney*dy,Regressed_function)
        plt.show()

        Split_width = np.abs(lx00-rx00) + aa + cc



