import numpy as np
import Wave_Sim_reader as rd
from matplotlib import pyplot as plt
import os
from scipy.optimize import curve_fit
import time

this_sig0 = 0.025/rd.dx


def GaussSingleAnalysis(Mat,timepoint:int,position:int,sig0:int = this_sig0):

    I,J,T = np.shape(Mat)

    #Model for the gaussian curve
    def model_f(x,a,b,x0):
        return b*np.exp(-((x - x0)**2 / a**2))
                
    SigMatBig = np.abs(np.array(Mat[position,:,timepoint]))
    n = 1
                    
    for tt in range(timepoint+1,int(timepoint+2*rd.lambdt/rd.dt/10)):
        SigMatBigTemp = np.abs(np.array(Mat[position,:,tt]))
        SigMatBig += SigMatBigTemp
        n += 1
    
    x_ting = np.linspace(0,J,J)

    #Takes the averages
    SigMat = (SigMatBig/n)
    
    # plt.figure()
    # plt.plot(x_ting,SigMat)
    # plt.draw()
    # plt.pause(0.1)
    # time.sleep(3)
    # plt.clf()


    
    popt, pcov = curve_fit(model_f,x_ting,SigMat,p0=[sig0,1,J/2])
    perr = np.sqrt(np.diag(pcov))
    # print(popt)


    
    sig,amp,x0 = popt
    dsig,damp,dx0 = perr
        
    return sig,amp,dsig,damp

    
def Blobdispersion(File:str,timepoint:int,varmode:str = 'Density',sig0:int = this_sig0,plotname:str = 'Blob Dispersion',plottitle = ''):
    
    hundred = int(timepoint/100)
    timepoint_s = int(timepoint-hundred*100)
    fallback = os.getcwd()

    os.chdir(File)

    files = os.listdir()
    sims = [sim for sim in files if not sim.endswith('.png')]
    sims = sorted(sims)
    sigb = np.array([])
    dens = np.array([])
    width = np.array([])
    dsigm = np.array([])
    # print(files)

    for i in range(len(sims)):
        
        if varmode == 'Density':
            dens = np.append(dens,(i)*10**18)
        elif varmode == 'Width':
            if i == 0:
                sigb = np.append(sigb,0)
            else:
                sigb = np.append(sigb,1+10*i)

        matrix,hundreddddd = rd.ReadBigSim(sims[i],hundred)
        # print(np.shape(matrix))

        I,J,T = np.shape(matrix)

        sig,amp,dsig,damp = GaussSingleAnalysis(matrix,timepoint_s,position=int(I*2/3),sig0 = sig0)
        # print(sig)
        width = np.append(width,2*sig)
        dsigm = np.append(dsigm,2*dsig)
        
    
    if varmode == 'Density':
        variable = dens
    elif varmode == 'Width':
        variable = sigb*rd.dy*2
        varmode = varmode + ' [m]'
    # print(width)
    

    plt.figure()
    plt.plot(variable,np.abs(width*rd.dy))
    # plt.scatter(variable,width*rd.dy)
    plt.errorbar(variable,np.abs(width*rd.dy),dsigm*rd.dy)
    plt.title(plottitle)
    plt.xlabel('Plasma' + varmode)
    plt.ylabel('Apparent width at position ' + str(round(int(I*2/3)*rd.dx,3)) + ' [m]')
    plt.savefig(plotname + str('.png'))

    plt.show()
    