import os
import numpy as np

def loadlattice_LEIR():
        #os.system("tail -n+49 /eos/home-m/mzampeta/SWAN_projects/LEIR/IBS_Analytic/madx/output/twiss_length.tfs > ./lib/twiss_for_IBS.tfs")
        data = np.loadtxt('./lib/twiss_for_IBS.tfs' , usecols=(1,2,3,4,6,7,8,9,11,12))
        #column=name,s,l,betx,alfx,mux,x,px,dx,dpx,bety,alfy,muy,y,py,dy,dpy;
        Posi1  = np.array(data[:,0])
        dels1  = np.array(data[:,1])
        betx1  = np.array(data[:,2])
        alfx1  = np.array(data[:,3])
        etax1  = np.array(data[:,4])
        etadx1 = np.array(data[:,5])
        bety1  = np.array(data[:,6])
        alfy1  = np.array(data[:,7])
        etay1  = np.array(data[:,8])
        etady1 = np.array(data[:,9])

        Posi  = np.array(Posi1[0])
        betx  = np.array(betx1[0])
        alfx  = np.array(alfx1[0])
        etax  = np.array(etax1[0])
        etadx = np.array(etadx1[0])
        bety  = np.array(bety1[0])
        alfy  = np.array(alfy1[0])
        etay  = np.array(etay1[0])
        etady = np.array(etady1[0])
        dels  = np.array(dels1[0])

        DimT = len(Posi1)
        item = 1
        while item < DimT :
            if abs(Posi1[item] - Posi1[item-1] >= 0.005):
                Posi  = np.append(Posi , Posi1[item])
                betx  = np.append(betx , betx1[item])
                alfx  = np.append(alfx , alfx1[item])
                etax  = np.append(etax , etax1[item])
                etadx = np.append(etadx, etadx1[item])
                bety  = np.append(bety , bety1[item])
                alfy  = np.append(alfy , alfy1[item])
                etay  = np.append(etay , etay1[item])
                etady = np.append(etady, etady1[item])
                dels  = np.append(dels , dels1[item])
            item += 1
        DimT = len(Posi)
        print('Number of Elements:', DimT)
        return Posi, dels, betx, alfx, etax, etadx, bety, alfy, etay, etady
