import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import signal

class Filter: #same class as Driver

        def __init__(self, fp , w, Lctime,nPulses,Dt): 

            self.fp = fp
            self.w = w
            self.Lctime = Lctime
            self.nPulses = nPulses
            self.upperTime = nPulses/fp
            self.Dt = Dt
            self.time = np.arange(0,self.upperTime,  self.Dt )
            self.nPointsPerPulse = 1/fp/self.Dt   
            self.L = len(self.time) 
            
        def pulses(self, Driver):
            
            t_norm = self.time % (1/self.fp)
            sig = np.zeros(self.time.size)
            risingEdge_indeces = (t_norm < self.w)
            
            for i in range(0,sig.size):
                if ( risingEdge_indeces[i] == 1 ):
                    sig[i] = 1 - math.exp(-t_norm[i] / self.Lctime)          
           
            fallingEdge_indeces = np.invert(risingEdge_indeces)
            sig_max = 1 - math.exp(-self.w / self.Lctime)    
            for i in range(0,sig.size):
                if (fallingEdge_indeces[i] == 1):
                    sig[i] = sig_max * math.exp( -(t_norm[i] - self.w ) / self.Lctime ) 
            return sig   
        
        def wdm_sigfilter(self,Driver,field,filterBandwidth,filterOrder,wdm_shift,ftype):
    
            nOfSamples = 2*len(field)-1
            samplingFrequency   = 1 / self.Dt   
            freq = np.linspace(-1/2, 1/2,nOfSamples)*samplingFrequency
            single_freq = freq[1]-freq[0]
            
            if (ftype == "Butterworth"):
                a , b = signal.butter(filterOrder, filterBandwidth,"low", analog = 'true') 
            elif (ftype == "Bessel"):
                a , b = signal.bessel(filterOrder , filterBandwidth , "low", analog = 'true')
            #filterBandwidth /(smaplingfrequancy/2)
            W,H  =  signal.freqs(a,b,freq)
            if( wdm_shift >= 0 ):
                H_z = np.zeros(round(wdm_shift/single_freq) )
                H_shift =np.concatenate((H_z, H))
                H_shift = H_shift[0:nOfSamples]

            if (wdm_shift < 0):
                H_z = np.zeros(abs(round(wdm_shift/single_freq)))
                H_shift = np.concatenate((H,H_z))
                H_shift = H_shift[0+abs(round(wdm_shift/single_freq)):len(H_shift)]

            norm_pulses =np.zeros(nOfSamples,dtype = complex)
            for i in range(len(field)):
                norm_pulses[i] = field[i]
            

            wdmfield =  (np.fft.ifft( np.fft.ifftshift(H_shift *np.fft.fftshift( np.fft.fft(norm_pulses)  )  ) ))           
            wdmfield = wdmfield[0:len(field)]
           
            
            return wdmfield
        
        def photodetector(self,Driver,field,filterBandwidth,filterOrder,wdm_shift,ftype):
    
            nOfSamples = 2*len(field)-1
            samplingFrequency   = 1 / self.Dt   
            freq = np.linspace(-1/2, 1/2,nOfSamples)*samplingFrequency
            single_freq = freq[1]-freq[0]
            
            if (ftype == "Butterworth"):
                a , b = signal.butter(filterOrder, filterBandwidth,"low", analog = 'true') 
            elif (ftype == "Bessel"):
                a , b = signal.bessel(filterOrder , filterBandwidth , "low", analog = 'true')
            #filterBandwidth /(smaplingfrequancy/2)
            W,H  =  signal.freqs(a,b,freq)
            if( wdm_shift >= 0 ):
                H_z = np.zeros(round(wdm_shift/single_freq) )
                H_shift = np.concatenate((H_z, H))
                H_shift = H_shift[0:nOfSamples]
 
            if (wdm_shift < 0):
                H_z = np.zeros(abs(round(wdm_shift/single_freq)))
                H_shift = np.concatenate((H,H_z))
               
                H_shift = H_shift[0+abs(round(wdm_shift/single_freq)):len(H_shift)]
 
            norm_pulses =np.zeros(nOfSamples)
            for i in range(len(field)):
                norm_pulses[i] = abs(field[i])**2
                            
            Qout =  (np.fft.ifft( np.fft.ifftshift(H_shift *np.fft.fftshift( np.fft.fft( norm_pulses)  )  ) )).real
            Qout = Qout[0:len(field)]
                       
            return Qout
 
        def laser_spectrum(self,Driver,E,Q1,phi1,filterBandwidth,filterOrder,nPointsPerPulse,nPulse1,nPulse2,ftype):
            
            arrsize = int(2*nPointsPerPulse-1) #increased size of the input
            print(arrsize)
            
            samplingFrequency   = 1 / self.Dt   
            freq = np.linspace(-1/2, 1/2,arrsize)*samplingFrequency
            
            #getting field of the first pulse from the array
            field = E[int(nPointsPerPulse * (nPulse1-1) ) : int(nPointsPerPulse * nPulse1)]
                       
            norm_pulses =np.zeros(arrsize, dtype = complex)#,dtype = complex) #adding zero's for higher accurancy
            for i in range(0,int(nPointsPerPulse)):
                norm_pulses[i] = field[i]

            Spectrum = abs(np.fft.fftshift( np.fft.fft( norm_pulses)  ) )**2 
           
            field2 = E[int(nPointsPerPulse * (nPulse2-1) ) : int(nPointsPerPulse * nPulse2)]
            norm_pulses2 =np.zeros(arrsize,dtype = complex)
            for i in range(0,int(nPointsPerPulse)): #adding zero's for higher accurancy
                norm_pulses2[i] = field2[i]
            Spectrum2 = abs(np.fft.fftshift( np.fft.fft( norm_pulses2)  ) )**2
            
            datawrite = []
            datawrite.append(freq)
            datawrite.append(Spectrum)   
            datawrite.append(Spectrum2)  
            text = np.c_[np.stack(datawrite, axis = 1)]  
            np.savetxt('last_spectrum.txt', text)
            output = [Spectrum, Spectrum2,freq]
            return output
              
        def wdm_bessel(self,Driver,pulse,filterBandwidth,filterOrder,centralfreq,wdm_shift ):
    
            nOfSamples = 2*self.time.size-1
            samplingFrequency   = 1 / self.Dt   
            freq = np.linspace(-1/2, 1/2,nOfSamples)*samplingFrequency
            single_freq = freq[1]-freq[0]
            lowfreq = centralfreq - filterBandwidth/2
            highfreq = centralfreq + filterBandwidth/2
            a , b = signal.bessel(filterOrder , filterBandwidth , "low", analog = 'true') # is ok

            W,H  =  signal.freqs(a,b,freq)
            if( wdm_shift >= 0 ):
                H_shift = np.zeros(round(wdm_shift/single_freq) )
                H_shift2 =np.concatenate((H_shift, H))
                H_shift2 = H_shift2[0:2*self.time.size-1]
  
            plt.plot(freq,abs(H_shift2)**2)
            plt.show()
   
            norm_pulses =np.zeros(nOfSamples)
            for i in range(self.time.size):
                norm_pulses[i] = pulse[i]
            

            if (len(H) != len(norm_pulses)):
                print("array error")
                return 0
            I = self.Ib + self.Ip * (np.fft.ifft( np.fft.ifftshift(H *np.fft.fftshift( np.fft.fft( norm_pulses)  )  ) )).real
            I = I[0:self.time.size]
            return I      
        