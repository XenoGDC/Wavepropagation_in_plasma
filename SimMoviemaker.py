from matplotlib import pyplot as plt # Needed to create the plots and save them as png's
import cv2 # Needed to read the images
import os # Needed to make a folder and to read said folder
import shutil #Needed to delete folder
import Wave_Sim_reader as rd
from time import time
import numpy as np
import traceback

dx,dy = rd.dx,rd.dy


def SimMatMoviemaker(Matrix,EndFrame: int, Moviename: str,StartFrame: int = 0, Pmatrix = None):
    t1 = time()
    #renaming some variables
    Mat = Matrix
    name = Moviename
    Start = StartFrame
    Stop = EndFrame

    J,I,T = np.shape(Matrix)
    
    print('\nMaking movie for ' + name + '\n')

    plasma_present = False

    try:
        if Pmatrix.all() != None:
            plasma_present = True
    except:
        pass

    # Sets up variables to remember the path and the file name for holding the images
    Path0 = os.getcwd()
    file = 'throwaway_img_folder'

    # Checks if the file is there
    files = os.listdir(os.getcwd())
    if not file in files:
            
        print('Making throwaway folder')
        os.mkdir(file)
    else:
        # Removes the folder used to store the images to save space
        shutil.rmtree(file)
        print('Making new throwaway folder')
        os.mkdir(file)

        
    t0 = time()

    xlist = np.linspace(-I/2,I/2,I)*dx
    ylist = np.linspace(0,J,J)*dy


    # Makes the renditions of the field in the matrix and makes a png for each time point
    for t in range(Start,Stop):
        plt.figure(figsize=[16/2,9/2])
        map1 = plt.pcolormesh(xlist,ylist,Mat[:,:,t],cmap='seismic',vmin=-1,vmax=1)
        cbar = plt.colorbar(map1)
        fieldname = str(name[0:2] + ' intensity')
        cbar.set_label(fieldname)

        # Adds a contour plot in case a plasma density has been specified
        if plasma_present:
            cmap = plt.contour(xlist[:-2],ylist[:-2],Pmatrix)
            conbar = plt.colorbar(cmap)
            conbar.set_label('Plasma density')
            # plt.tight_layout()

        plt.title(name)
        plt.xlabel('y [m]')
        plt.ylabel('x [m]')

        l1 = len(str(t))
        l2 = len(str(Stop))
        length = l2-l1 + 1
        imgname = 'img'
        for l in range(length):
            imgname = str(imgname + '0')

        imgname = str(imgname + str(t))
        plt.savefig(os.path.join(file, imgname))
        plt.close()

        t1 = time()
        diff = t1-t0
        minut = int((diff/((t+1)/Stop) - diff)/60)
        sec = round(((diff/((t+1)/Stop) - diff)/60-minut)*60)
        print(str(round(t/Stop*100)) + '% done with Etr: ' + str(minut) + 'min ' + str(sec) + 'sec',end='\r')


    # Makes a list for the images and only uses the png's
    imglist = sorted(os.listdir(os.path.join(Path0,file)))
    
    images = [img for img in imglist
              if img.endswith('.png')]

    # Reads the first frame to determine dimensions
    frame = cv2.imread(os.path.join(file,imglist[0]))
    height, width, layers = frame.shape

    # Sets up the video
    vidname = str(str(name) + '.avi')
    video = cv2.VideoWriter(vidname, 0 , 10, (width,height)) # Gives the video the determined name and sets the fps to 10

    # Writes the images into the video file
    for image in images:
        img = cv2.imread(os.path.join(file,image))
        video.write(img)
    
    # Frees up memory I think
    cv2.destroyAllWindows()
    # Releases the video as a file
    video.release()

    # Removes the folder used to store the images to save space
    shutil.rmtree(file)

    t2 = time() - t1
    print('Elapsed time: ' + str(int(round(t2))) + 'sec\n')

def concatenate_videos(new_video_path, videos):

    vid = cv2.VideoCapture(videos[0])
    while vid.isOpened():
        r,img = vid.read()
        break

    height,width,layers = img.shape

    video = cv2.VideoWriter(new_video_path, 0, 12, (width,height))

    for v in videos:
        curr_v = cv2.VideoCapture(v)
        while curr_v.isOpened():
            r, frame = curr_v.read()
            if not r:
                break
            video.write(frame)

    video.release()

    cv2.destroyAllWindows()

    pass

def BigMovieMaker(Filename:str,Pmatrix = None):
    Traceback = os.getcwd()
    os.chdir(Filename)
    plasma_present = False
    print(np.shape(Pmatrix))
    try:
        if Pmatrix.all() != None:
            print(np.shape(Pmatrix))
            plasma_present = True
    except:
        traceback.print_exc()
        pass

    filelist = os.listdir()
    files = [data for data in filelist if data.endswith('.npy')]

    
    # print(files)
    for i in range(len(files)):
        name = files[i]
        mat = rd.SimReader(name)
        videoname = name[:-4]
        SimMatMoviemaker(mat,len(mat[0,0,:]),videoname,Pmatrix=Pmatrix)

    filelist = os.listdir()
    vidlist = [vid for vid in filelist if vid.endswith('99.avi')]
    vidlist.sort

    try:
        concatenate_videos(str(Filename +'.avi'),vidlist)
    except:
        traceback.print_exc()
        input()
        pass

    for oldvid in vidlist:
        os.remove(oldvid)
    cv2.destroyAllWindows()

    os.chdir(Traceback)

