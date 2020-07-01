import string
import numpy as np
from matplotlib import pylab as plt
from matplotlib.pylab import *
#import 
from astropy.io import fits
import scipy.optimize
from numpy import *

class data():
       
    def __init__(self, filename):
        """Define the data class"""
        """self.lines"""
        """self.elts"""
        """self.data"""
        """ """
        """self.lon"""
        """self.scans"""
        """self.Nchan"""
        """self.chan"""
        """self.velo"""
        """self.Ta[scan]  a dictionary"""
        #
        f = open(filename, 'r')
        self.lines=f.readlines()
        f.close()
        self.lines = list(map(str.strip,self.lines))
        self.elts = list(map(str.split,self.lines))
        for i in range(len(self.lines)):
            if (len(self.elts[i]) > 6) and ((self.elts[i][5]=='galactic') or (self.elts[i][5]=='CygEMN') or (self.elts[i][5]=='Sun') ) :
                self.lon=float(self.elts[i][6])
                self.lat=float(self.elts[i][7])

        Nlines=0
        for i in range(len(self.lines)):
            if ((self.elts[i][0][0:2])=='20'):
                self. deltaf=float(self.elts[i][6])
                Nlines=Nlines+1
                self.Nchan = int(self.elts[i][8])
                self.mydate=self.elts[i][0][:8]
                self.myt=self.elts[i][0][9:]
                self.Az=float(self.elts[i][1])
                self.El=float(self.elts[i][2])
                self.freq0=float(self.elts[i][5])
                self.Vlsr=float(self.elts[i][self.Nchan+10])

        self.data=np.zeros((Nlines,self.Nchan),float)
        Nlines=0
        self.Ta={}
        for i in range(len(self.lines)):
            if ((self.elts[i][0][0:2]) == '20'):
                for j in range(self.Nchan):
                    self.data[Nlines,j] = float(self.elts[i][9+j])
                #self.Ta["scan" + str(Nlines)]=self.data[Nlines]
                self.Ta[Nlines]=self.data[Nlines]
                Nlines=Nlines+1
        self.scans=self.Ta.keys()

        self.chan=range(1,self.Nchan+1)
        c=299792.458
        self.deltaV=-(self.deltaf*c/(self.freq0)) #+int(self.Nchan/2.)*self.deltaf)
        self.Nchan0=int(self.Nchan/2.)
        self.velo=np.zeros(self.Nchan,float)


        #Alex
        frest=1420.4
        self.Voffset=((frest-(self.freq0+self.Nchan0*self.deltaf))/frest)*c
        self.Voffset=0.

        for j in range(self.Nchan):
            V0=c*((self.deltaf*self.Nchan0)/self.freq0)
            self.velo[j]=V0+j*self.deltaV
            #self.velo[j]=-self.Nchan0*self.deltaV+j*self.deltaV-self.Vlsr
        

    def onpick(self,event):
        """Used for zooming in plotscans """
        thisline = event.artist
        xdata, ydata = thisline.get_data()
        ind = event.ind
        figi = figure(111)
        ax = figi.add_subplot(1,1,1)
        ax.plot(xdata[ind], ydata[ind])
        figi.show()
        return True

    def plotscans(self, scanList):
        """Plot all scans in individual subplots. Clic on one scan to popup a zoom. One need to close the popup window to zoom another scan"""
        #
        figsrc=plt.figure(100, figsize=(14, 11))
        figsrc.subplots_adjust(hspace=0.5)
        figsrc.subplots_adjust(wspace=0.3)
        Lscan=len(scanList)
        if (Lscan <= 16):
            Nx=4
        else:
            Nx=8
        Ny=int(Lscan/Nx)+1
        #print(Lscan, Nx, Ny)
        for scan in scanList:
            axsrc = figsrc.add_subplot(Nx,Ny, scan+1)
            axsrc.plot(self.velo, self.Ta[scan], picker=50000000)
            plt.xticks(fontsize=4)
            plt.yticks(fontsize=8)
            #plt.ylabel("Ta [K]")
            #plt.xlabel("velocity [km/s]")
            plt.title("Scan"+str(scan), fontsize=10)
            
        #figsrc.canvas.mpl_connect('pick_event', self.onpick)
        #plt.savefig('plotscans.ps')
        plt.show()
        return

    def plotallscans(self, scanList):
        """Plot all scan of the scan list together """
        #
        figsrc=plt.figure(100)
        for scan in scanList:
            plot(self.velo[8:self.Nchan-8], self.Ta[scan][8:self.Nchan-8])
            plt.ylabel("Ta [K]")
            plt.xlabel("velocity [km/s]")
            plt.title("Scan"+str(scanList))
     
        plt.show()
        #plt.savefig('./figures/plotallscans.ps')
        return

    def plotallscansChan(self, scanList):
        """Plot all scan of the scan list together """
        #
        for scan in scanList:
            plot(self.chan[8:self.Nchan-8], self.Ta[scan][8:self.Nchan-8])
            plt.ylabel("Ta [K]")
            plt.xlabel("Channels")
            plt.title("Scan"+str(scanList))
            
        #plt.savefig('./figures/plotallscansChan.ps')
        plt.show()
        return

    def baseline1w(self, scanList, w1, w2):
        """Remove a baseline for one scan between channels w1 and w2"""
        #
        for scan in scanList:
            x=w1+arange(len(self.Ta[scan][w1:w2])+len(self.Ta[scan][w3:w4]))
            y=self.Ta[scan][w1:w2]
            fitfunc = lambda q, x: q[0]*x+q[1]
            errfunc = lambda q, x, y: fitfunc(q,x)-y
            q = scipy.c_[0.5,mean(self.Ta[scan][w1:w2])]
            
            mean(self.Ta[scan][w1:w2])
            q1, success = scipy.optimize.leastsq(errfunc, q.copy()[0], args=(x,y))
            corrfit = fitfunc(q1, x)
            #
            x=arange(len(self.Ta[scan]))
            corrfit = fitfunc(q1, x)
            xw1=w1+np.zeros(10,float)
            xw2=w2+np.zeros(10,float)
            yw1=-50+40*arange(10)
            yw2=-50+40*arange(10)
            plt.figure(1)
            plt.plot(xw1,yw1)
            plt.plot(xw2,yw2)
            plt.plot(self.Ta[scan])
            plt.plot(corrfit, 'b-')
            self.Ta[scan]=self.Ta[scan]-corrfit
        #plt.show()


    def baseline2w(self, scanList, w1, w2, w3, w4):
        """Remove a baseline for one scan between channels w1 and w2"""
        #
        for scan in scanList:
            x=append(w1+arange(len(self.Ta[scan][w1:w2])), w3+arange(len(self.Ta[scan][w3:w4])))
            y=append(self.Ta[scan][w1:w2], self.Ta[scan][w3:w4])
            #print(x)
            #print(y)
            fitfunc = lambda q, x: q[0]*x+q[1]
            errfunc = lambda q, x, y: fitfunc(q,x)-y
            q = scipy.c_[0.5,mean(self.Ta[scan][w1:w2])]
            q1, success = scipy.optimize.leastsq(errfunc, q.copy()[0], args=(x,y))
            corrfit = fitfunc(q1, x)
            #
            x=arange(len(self.Ta[scan]))
            corrfit = fitfunc(q1, x)
            xw1=w1+np.zeros(10,float)
            xw2=w2+np.zeros(10,float)
            xw3=w3+np.zeros(10,float)
            xw4=w4+np.zeros(10,float)
            yw1=-50+40*arange(10)
            yw2=-50+40*arange(10)
            plt.figure(2)
            plt.plot(xw1,yw1)
            plt.plot(xw2,yw2)
            plt.plot(xw3,yw1)
            plt.plot(xw4,yw2)
            plt.plot(self.Ta[scan])
            plt.plot(corrfit, 'b-')
            self.Ta[scan]=self.Ta[scan]-corrfit
        #plt.savefig('./figures/baseline2w.ps')
        plt.show()

    def baseline2w2(self, scanList, w1, w2, w3, w4):
        """Remove a baseline for one scan between channels w1 and w2"""
        #
        for scan in scanList:
            x=append(w1+arange(len(self.Ta[scan][w1:w2])), w3+arange(len(self.Ta[scan][w3:w4])))
            y=append(self.Ta[scan][w1:w2], self.Ta[scan][w3:w4])
            #print x
            #print y
            fitfunc = lambda q, x: q[0]*x**2+q[1]*x+q[2]
            errfunc = lambda q, x, y: fitfunc(q,x)-y
            q = scipy.c_[-0.5,0.5,mean(self.Ta[scan][w1:w2])]
            q1, success = scipy.optimize.leastsq(errfunc, q.copy()[0], args=(x,y))
            corrfit = fitfunc(q1, x)
            #
            x=arange(len(self.Ta[scan]))
            corrfit = fitfunc(q1, x)
            xw1=w1+np.zeros(10,float)
            xw2=w2+np.zeros(10,float)
            xw3=w3+np.zeros(10,float)
            xw4=w4+np.zeros(10,float)
            yw1=-50+40*arange(10)
            yw2=-50+40*arange(10)
            plt.figure(2)
            plt.plot(xw1,yw1)
            plt.plot(xw2,yw2)
            plt.plot(xw3,yw1)
            plt.plot(xw4,yw2)
            plt.plot(self.Ta[scan])
            plt.plot(corrfit, 'b-')
            self.Ta[scan]=self.Ta[scan]-corrfit
        #plt.savefig('./figures/baseline2w.ps')
        #plt.show()


    def baseline2w3(self, scanList, w1, w2, w3, w4):
        """Remove a baseline for one scan between channels w1 and w2"""
        #
        for scan in scanList:
            x=append(w1+arange(len(self.Ta[scan][w1:w2])), w3+arange(len(self.Ta[scan][w3:w4])))
            y=append(self.Ta[scan][w1:w2], self.Ta[scan][w3:w4])
            #print x
            #print y
            fitfunc = lambda q, x: q[0]*x**3+q[1]*x**2+q[2]*x+q[3]
            errfunc = lambda q, x, y: fitfunc(q,x)-y
            q = scipy.c_[0.5,0.5,0.5,mean(self.Ta[scan][w1:w2])]
            q1, success = scipy.optimize.leastsq(errfunc, q.copy()[0], args=(x,y))
            corrfit = fitfunc(q1, x)
            #
            x=arange(len(self.Ta[scan]))
            corrfit = fitfunc(q1, x)
            xw1=w1+np.zeros(10,float)
            xw2=w2+np.zeros(10,float)
            xw3=w3+np.zeros(10,float)
            xw4=w4+np.zeros(10,float)
            yw1=-50+40*arange(10)
            yw2=-50+40*arange(10)
            plt.figure(2)
            plt.plot(xw1,yw1)
            plt.plot(xw2,yw2)
            plt.plot(xw3,yw1)
            plt.plot(xw4,yw2)
            plt.plot(self.Ta[scan])
            plt.plot(corrfit, 'b-')
            self.Ta[scan]=self.Ta[scan]-corrfit
        #plt.savefig('./figures/baseline2w.ps')
        plt.show()

    def baseline2w4(self, scanList, w1, w2, w3, w4):
        """Remove a baseline for one scan between channels w1 and w2"""
        #
        for scan in scanList:
            x=append(w1+arange(len(self.Ta[scan][w1:w2])), w3+arange(len(self.Ta[scan][w3:w4])))
            y=append(self.Ta[scan][w1:w2], self.Ta[scan][w3:w4])
            #print x
            #print y
            fitfunc = lambda q, x: q[0]*x**4+q[1]*x**3+q[2]*x**2+q[3]*x+q[4]
            errfunc = lambda q, x, y: fitfunc(q,x)-y
            q = scipy.c_[0.5,0.5,0.5,0.5,mean(self.Ta[scan][w1:w2])]
            q1, success = scipy.optimize.leastsq(errfunc, q.copy()[0], args=(x,y))
            corrfit = fitfunc(q1, x)
            #
            x=arange(len(self.Ta[scan]))
            corrfit = fitfunc(q1, x)
            xw1=w1+np.zeros(10,float)
            xw2=w2+np.zeros(10,float)
            xw3=w3+np.zeros(10,float)
            xw4=w4+np.zeros(10,float)
            yw1=-50+40*arange(10)
            yw2=-50+40*arange(10)
            plt.figure(2)
            plt.plot(xw1,yw1)
            plt.plot(xw2,yw2)
            plt.plot(xw3,yw1)
            plt.plot(xw4,yw2)
            plt.plot(self.Ta[scan])
            plt.plot(corrfit, 'b-')
            self.Ta[scan]=self.Ta[scan]-corrfit
        plt.savefig('./figures/baseline2w.ps')
        #plt.show()


    def average(self, scanList):
        """ Average several scans. Ensure that the baselines are removed"""
        self.scanMean=np.zeros(self.Nchan, float)
        for ichan in range((self.Nchan)):
            for scan in scanList:
                #self.scanMean[ichan]=mean(self.Ta[scan][ichan])
                self.scanMean[ichan]+=self.Ta[scan][ichan]
            self.scanMean[ichan]/=len(scanList)
        #plt.plot(self.chan[8:self.Nchan-8], self.scanMean[8:self.Nchan-8])
        #plt.savefig('./average.ps')
        ##plt.show()

    def averageTemp(self, scanList):
        """ Average several scans and print average temperature. Ensure that the baselines are removed"""
        self.scanMean=np.zeros(self.Nchan, float)
        for ichan in range((self.Nchan)):
            for scan in scanList:
                #self.scanMean[ichan]=mean(self.Ta[scan][ichan])
                self.scanMean[ichan]+=self.Ta[scan][ichan]
            self.scanMean[ichan]/=len(scanList)
        print("Average temperature : %0.2f K" % mean(self.scanMean[10:self.Nchan-10]))
        plt.plot(self.chan, self.scanMean)
        plt.text(20,40,"Bords enleves pour calculer la temperature moyenne.")
        plt.show()

    def plotAveChan(self):
        plt.figure(3)
        plt.plot(self.chan[8:self.Nchan-8], self.scanMean[8:self.Nchan-8])
        plt.ylabel("Ta [K]")
        plt.xlabel("Channels")
        #plt.title("Long"+str(self.lon), fontsize=13)
        #plt.savefig('./figures/averageC.ps')
        plt.show()

    def plotAveVel(self):
        plt.figure(4)
        plt.plot(self.velo[8:self.Nchan-8], self.scanMean[8:self.Nchan-8])
        plt.ylabel("Ta [K]")
        plt.xlabel("Velocity [km/s]")
        #plt.title("Long"+str(self.lon), fontsize=13)
        #plt.savefig('./figures/averageV.ps')
        plt.show()

    def toFits(self, filename, scan):
        """Convert one scan to Fits data format put into filename.fits"""
        #
        nbits=16
        #hdu = pyfits.PrimaryHDU(int16(self.Ta[scan]))
        #hdu.writeto('fits/'+str(filename)+'.fits',clobber=True)
        hdu = fits.PrimaryHDU(int16(self.Ta[scan]))
        hdu.writeto('fits/'+str(filename)+'.fits',clobber=True)
        
        hdulist = fits.open(str(filename)+'.fits',  mode='update')
        #hdulist = pyfits.open(str(filename)+'.fits',mode='update')
        prihdr = hdulist[0].header
        prihdr.update('BITPIX',nbits)
        prihdr.update('BSCALE',1)
        prihdr.update('BZERO',0)
        prihdr.update('BUNIT','K')
        prihdr.update('CTYPE1','FREQ')
        prihdr.update('CRVAL1',round((self.freq0+self.deltaf*self.Nchan/2.),2)*pow(10,6))
        prihdr.update('CDELT1',self.deltaf*pow(10,6))
        prihdr.update('CRPIX1',self.Nchan/2.)
        prihdr.update('CTYPE2','GLON')
        prihdr.update('CRVAL2',self.lon)
        prihdr.update('CDELT2',0.)
        prihdr.update('CRPIX2',0.)
        prihdr.update('CTYPE3','GLAT')
        prihdr.update('CRVAL3',self.lat)
        prihdr.update('CDELT3',0.)
        prihdr.update('CRPIX3',0.)
        prihdr.update('CTYPE4','STOKES')
        prihdr.update('CRVAL4',1.)
        prihdr.update('CDELT4',0.)
        prihdr.update('CRPIX4',0.)
        prihdr.update('TELESCOP','SRT-PARIS')
        prihdr.update('OBJECT','Milky Way')
        prihdr.update('GLAT',self.lat)
        prihdr.update('GLON',self.lon)
        prihdr.update('EPOCH',2000.0)
        prihdr.update('LINE','HI')
        prihdr.update('RESTFREQ',1420.405752*pow(10,6))
        prihdr.update('VELO-LSR',self.Vlsr)
        prihdr.update('DATE-OBS',self.mydate)
        prihdr.update('UT',self.myt)
        prihdr.update('ELEVATIO',self.El)
        prihdr.update('AZIMUTH',self.Az)
        hdulist.flush()


    def aveToFits(self, filename, mydata):
        """Convert the data to Fits data format put into filename.fits"""
        #
 
        nbits=16
        #hdu = pyfits.PrimaryHDU(int16(self.Ta[scan]))
        #hdu.writeto('fits/'+str(filename)+'.fits',clobber=True)
        hdu = fits.PrimaryHDU(int16(self.Ta[scan]))
        hdu.writeto('fits/'+str(filename)+'.fits',clobber=True)
        
        hdulist = fits.open(str(filename)+'.fits',  mode='update')
        #hdulist = pyfits.open(str(filename)+'.fits',mode='update')
        prihdr = hdulist[0].header
        prihdr.update('BITPIX',nbits)
        prihdr.update('BSCALE',1)
        prihdr.update('BZERO',0)
        prihdr.update('BUNIT','K')
        prihdr.update('CTYPE1','FREQ')
        prihdr.update('CRVAL1',round((self.freq0+self.deltaf*self.Nchan/2.),2)*pow(10,6))
        prihdr.update('CDELT1',self.deltaf*pow(10,6))
        prihdr.update('CRPIX1',self.Nchan/2.)
        prihdr.update('CTYPE2','GLON')
        prihdr.update('CRVAL2',self.lon)
        prihdr.update('CDELT2',0.)
        prihdr.update('CRPIX2',0.)
        prihdr.update('CTYPE3','GLAT')
        prihdr.update('CRVAL3',self.lat)
        prihdr.update('CDELT3',0.)
        prihdr.update('CRPIX3',0.)
        prihdr.update('CTYPE4','STOKES')
        prihdr.update('CRVAL4',1.)
        prihdr.update('CDELT4',0.)
        prihdr.update('CRPIX4',0.)
        prihdr.update('TELESCOP','SRT-PARIS')
        prihdr.update('OBJECT','Milky Way')
        prihdr.update('GLAT',self.lat)
        prihdr.update('GLON',self.lon)
        prihdr.update('EPOCH',2000.0)
        prihdr.update('LINE','HI')
        prihdr.update('RESTFREQ',1420.405752*pow(10,6))
        prihdr.update('VELO-LSR',self.Vlsr)
        prihdr.update('DATE-OBS',self.mydate)
        prihdr.update('UT',self.myt)
        prihdr.update('ELEVATIO',self.El)
        prihdr.update('AZIMUTH',self.Az)
        hdulist.flush()

