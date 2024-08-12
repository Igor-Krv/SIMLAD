import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy import signal

class Driver:

        def __init__(self, Ip , Ib, fp , w, Lctime,nPulses,Dt):

            self.Ip = Ip
            self.Ib = Ib
            self.fp = fp
            self.w = w
            self.Lctime = Lctime
            self.nPulses = nPulses
            self.upperTime = nPulses/fp
            self.Dt = Dt
            self.time = np.arange(0,self.upperTime,  self.Dt )
            #print("driver time array length : " +str(len(self.time)) )
            self.nPointsPerPulse = 1/fp/self.Dt   
            self.L = len(self.time)  
            self.arbtime =[]
            
        def pulses(self, Driver, form = None):
            
            t_norm = self.time % (1/self.fp)
            sig = np.zeros(self.time.size)
            risingEdge_indeces = (t_norm < self.w)
            
            if(form == 'Cos'):
                for i in range(0,sig.size):
                    sig[i] = (1-np.cos(self.time[i]*self.fp))
                   
            else:
                for i in range(0,sig.size):
                    if ( risingEdge_indeces[i] == 1 ):
                        sig[i] = 1 - math.exp(-t_norm[i] / self.Lctime)          
               
                fallingEdge_indeces = np.invert(risingEdge_indeces)
                sig_max = 1 - math.exp(-self.w / self.Lctime)    
                for i in range(0,sig.size):
                    if (fallingEdge_indeces[i] == 1):
                        sig[i] = sig_max * math.exp( -(t_norm[i] - self.w ) / self.Lctime ) 
            return sig  
        def cospulses(self, Driver):
            
            t_norm = self.time % (1/self.fp)
            sig = np.zeros(self.time.size)
            for i in range(0,sig.size):
                sig[i] = (1-np.cos(self.time[i]*self.fp))
            sig_max = 1 - np.cos(-self.w / self.Lctime)    

            return sig
        def cossigfilter(self,Driver,pulse,filterBandwidth,filterOrder,ftype):
    
            nOfSamples = 2*self.time.size-1
            samplingFrequency   = 1 / self.Dt   
            #freq = np.arange(-1/2, 1/2  , 1/nOfSamples)*samplingFrequency   #there was 1/2+ 1/nOfSamples
            freq = np.linspace(-1/2, 1/2,nOfSamples)*samplingFrequency

            
            a , b = signal.butter(filterOrder, filterBandwidth,"low", analog = 'true') # is ok
            if (ftype == "Bessel"):
                a , b = signal.bessel(filterOrder , filterBandwidth , "low", analog = 'true')
            W , H  =  signal.freqs(a,b,freq)  

            norm_pulses =np.zeros(nOfSamples)
            for i in range(self.time.size):
                norm_pulses[i] = pulse[i]
            
            if (len(H) != len(norm_pulses)):
                print("array error")
                return 0
            I = self.Ib + self.Ip/2 * (np.fft.ifft( np.fft.ifftshift(H *np.fft.fftshift( np.fft.fft( norm_pulses)  )  ) )).real
            I = I[0:self.time.size]
            return I      
            
            
        def arbpulses(self,Driver,Ib,Ip,w,Id,wd,w1,fp,shift,form,incline = None):
            
            Lctime = self.Lctime
            nPulses = 1
            upperTime = nPulses/fp
            Dt = self.Dt
            time = np.arange(0,upperTime,Dt )
            #print("arb time array length : " +str(len(time)) )
            L = len(time)
            loctime = L
            nPointsPerPulse = 1/fp/Dt   
            
            t_norm = time % (1/fp)
            sig = np.zeros(int(nPointsPerPulse))
            risingEdge_indeces = (t_norm < w)
            
            t_norm = time % (1/fp)
            
            sig = np.zeros(time.size)
            risingEdge_indeces = (t_norm < w)
            if (form == 'cos'):
                for i in range(0,sig.size):
                    sig[i] = (1-np.cos(self.time[i]*self.w))
            else:
                for i in range(0,sig.size):
                    if ( risingEdge_indeces[i] == 1 ):
                        sig[i] = 1 - math.exp(-t_norm[i] / Lctime)         
               
                fallingEdge_indeces = np.invert(risingEdge_indeces) 
                sig_max = 1 - math.exp(-w / Lctime)    
                for i in range(0,sig.size):
                    if (fallingEdge_indeces[i] == 1):
                        sig[i] = sig_max * math.exp( -(t_norm[i] - w ) / Lctime ) 
            sigadd = np.zeros(int(shift/Dt))
            sigadd=np.concatenate([sigadd,sig])
            sigadd = sigadd[0:len(sig)]
            k_rise = np.linspace(1, incline,len(time))

            #print(str(incline) + " pulse incline")
    
            if (incline == None):
                return sigadd*Ip + Ib
            else:
                return sigadd*Ip*k_rise+Ib 
        def arbminor(self,Driver,w,Id,wd,w1,fp,filterBandwidth,filterOrder,ftype):
       
            Lctime = self.Lctime
            nPulses = 1
            upperTime = nPulses/fp
            Dt = self.Dt
            time = np.arange(0,upperTime,Dt )
           
            L = len(time)
            loctime = L
            nPointsPerPulse = 1/fp/Dt   
            
            t_norm = time % (1/fp)
            sig = np.zeros(int(nPointsPerPulse))
            risingEdge_indeces = (t_norm < wd)
            
            t_norm = time % (1/fp)
            sig = np.zeros(time.size)
            risingEdge_indeces = (t_norm < wd)
            
            for i in range(0,sig.size):
                if ( risingEdge_indeces[i] == 1 ):
                    sig[i] = 1 - math.exp(-t_norm[i] / Lctime)          
           
            fallingEdge_indeces = np.invert(risingEdge_indeces)       
            sig_max = 1 - math.exp(-wd / Lctime)    
            for i in range(0,sig.size):
                if (fallingEdge_indeces[i] == 1):
                    sig[i] = sig_max * math.exp( -(t_norm[i] - wd ) / Lctime ) 
            sigadd = np.zeros(int(w1/Dt))
            sigadd=np.concatenate([sigadd,sig])
            pulse = sigadd
            return sigadd*Id
        def arbincline(self,Driver,w,Ib,wd,w1,fp,filterBandwidth,filterOrder,ftype):
       
            Lctime = self.Lctime
            nPulses = 1
            upperTime = nPulses/fp
            Dt = self.Dt
            time = np.arange(0,upperTime,Dt )
            
            L = len(time)
            loctime = L
            nPointsPerPulse = 1/fp/Dt   
            
            t_norm = time % (1/fp)
            sig = np.zeros(int(nPointsPerPulse))
            risingEdge_indeces = (t_norm < wd)
            
            t_norm = time % (1/fp)
            sig = np.zeros(time.size)
            risingEdge_indeces = (t_norm < wd)
            
            for i in range(0,sig.size):
                if ( risingEdge_indeces[i] == 1 ):
                    sig[i] = 1 - math.exp(-t_norm[i] / Lctime)          
            fallingEdge_indeces = np.invert(risingEdge_indeces)
            sig_max = 1 - math.exp(-wd / Lctime)    
            for i in range(0,sig.size):
                if (fallingEdge_indeces[i] == 1):
                    sig[i] = sig_max * math.exp( -(t_norm[i] - wd ) / Lctime ) 
            sigadd = np.zeros(int(w1/Dt))
            sigadd=np.concatenate([sigadd,sig])

            pulse = sigadd
            return sigadd*Ib  
        
        def arbbetterfilter(self,Ib,Ip,Id,wd,w1,pulse,filterBandwidth,filterOrder,ftype):
            
            nOfSamples = 2*len(pulse)-1
            samplingFrequency   = 1 / self.Dt   
            freq = np.linspace(-1/2, 1/2,nOfSamples)*samplingFrequency
            a , b = signal.butter(filterOrder, filterBandwidth,"low", analog = 'true') # is ok
            if (ftype == "Bessel"):
                a , b = signal.bessel(filterOrder , filterBandwidth , "low", analog = 'true')
            W , H  =  signal.freqs(a,b,freq)  
            norm_pulses =np.zeros(nOfSamples)
            for i in range(len(pulse)):
                norm_pulses[i] = pulse[i]
            
            I =  (np.fft.ifft( np.fft.ifftshift(H *np.fft.fftshift( np.fft.fft( norm_pulses)  )  ) )).real
            
            I = I[0:len(pulse)]            
            return I 
        
        
        
        def arbfilter(self, Ib, Ip, pulse, filterBandwidth, filterOrder, ftype):
    
            nOfSamples = 2*len(pulse)-1
            samplingFrequency   = 1 / self.Dt   
            freq = np.linspace(-1/2, 1/2,nOfSamples)*samplingFrequency
            a , b = signal.butter(filterOrder, filterBandwidth,"low", analog = 'true') # is ok
            if (ftype == "Bessel"):
                a , b = signal.bessel(filterOrder , filterBandwidth , "low", analog = 'true')
            W , H  =  signal.freqs(a,b,freq)  
            norm_pulses =np.zeros(nOfSamples)
            for i in range(len(pulse)):
                norm_pulses[i] = pulse[i]
            
            if (len(H) != len(norm_pulses)):
                print("array error")
                return 0
            I = Ib + Ip * (np.fft.ifft( np.fft.ifftshift(H *np.fft.fftshift( np.fft.fft( norm_pulses)  )  ) )).real
            I = I[0:len(pulse)]
            return I 
        
            #I = Ib + Ip/2*(1 + cos())
        def sigfilter(self, Driver, pulse, filterBandwidth, filterOrder, ftype, form):
    
            nOfSamples = 2*self.time.size-1
            samplingFrequency   = 1 / self.Dt   
            freq = np.linspace(-1/2, 1/2,nOfSamples)*samplingFrequency
            
            
            a , b = signal.butter(filterOrder, filterBandwidth,"low", analog = 'true') # is ok
            if (ftype == "Bessel"):
                a , b = signal.bessel(filterOrder , filterBandwidth , "low", analog = 'true')
            W , H  =  signal.freqs(a,b,freq)  
            norm_pulses =np.zeros(nOfSamples)
            for i in range(self.time.size):
                norm_pulses[i] = pulse[i]
            
            if (len(H) != len(norm_pulses)):
                print("array error")
                return 0
            if (form == 'cos'):
                I = self.Ib + self.Ip/2 * (np.fft.ifft( np.fft.ifftshift(H *np.fft.fftshift( np.fft.fft( norm_pulses)  )  ) )).real
            else:
                I = self.Ib + self.Ip * (np.fft.ifft( np.fft.ifftshift(H *np.fft.fftshift( np.fft.fft( norm_pulses)  )  ) )).real
            I = I[0:self.time.size]
            return I  
        def bessel(self,Driver,pulse,filterBandwidth,filterOrder ):
    
            nOfSamples = 2*self.time.size-1
            samplingFrequency   = 1 / self.Dt   
            freq = np.linspace(-1/2, 1/2,nOfSamples)*samplingFrequency
            a , b = signal.bessel(filterOrder , filterBandwidth , "low", analog = 'true')
            W , H  =  signal.freqs(a,b,freq) 
            norm_pulses =np.zeros(nOfSamples)
            for i in range(self.time.size):
                norm_pulses[i] = pulse[i]
            
            if (len(H) != len(norm_pulses)):
                print("array error")
                return 0
            I = self.Ib + self.Ip * (np.fft.ifft( np.fft.ifftshift(H *np.fft.fftshift( np.fft.fft( norm_pulses)  )  ) )).real
            I = I[0:self.time.size]
            return I 
        
        def wdm_butter(self,Driver,pulse,filterBandwidth,filterOrder,centralfreq):
    
            nOfSamples = 2*self.time.size-1

            samplingFrequency   = 1 / self.Dt   
            freq = np.linspace(-1/2, 1/2,nOfSamples)*samplingFrequency

            lowfreq = centralfreq - filterBandwidth/2
            highfreq = centralfreq + filterBandwidth/2
            a , b = signal.butter(filterOrder, [lowfreq/(samplingFrequency/2),highfreq/(samplingFrequency/2)],"bandpass", 'ba',samplingFrequency,analog = 'true') # is ok
            W , H  =  signal.freqs(a,b,freq)   

            norm_pulses =np.zeros(nOfSamples)
            for i in range(self.time.size):
                norm_pulses[i] = pulse[i]

            if (len(H) != len(norm_pulses)):
                print("array error")
                return 0
            I = self.Ib + self.Ip * (np.fft.ifft( np.fft.ifftshift(H *np.fft.fftshift( np.fft.fft( norm_pulses)  )  ) )).real
            I = I[0:self.time.size]
            return I  
        
        # def getfig(self, x, y, root, scatter, figure, xl = None, yl = None, lim = None):
        #     try: 
        #         scatter.get_tk_widget().forget_pack()
        #     except AttributeError: 
        #         pass 
        #     scatter.get_tk_widget().delete("all")
        #     figure.clear()
        #     ax3 = figure.add_subplot(111)
        #     ax3.plot(x,y)
        #     ax3.set_xlabel('T')
        #     ax3.set_title('I ma')
        #     scatter3 = FigureCanvasTkAgg(figure, root) 
        #     scatter3.get_tk_widget().place(x=240,y=60,width=800,height=700)
        #     root.update()
            
        def getfig_n(self,toolbar,n,x,y,name,leg , root1,canvas, figure, shift = None,frame = None,W = None, H = None,L = None,Ib = None,Ip = None):
            try: 
                canvas.get_tk_widget().place_forget()
            except AttributeError: 
                pass
            
            figure.clear()
            root1.update()
            color = ["tab:orange","tab:green","tab:red","tab:blue"]
            axs = figure.subplots(n)
            left  = 0.11            # the left side of the subplots of the figure
            right = 0.783           # the right side of the subplots of the figure
            bottom = 0.297-0.004*n  # the bottom of the subplots of the figure
            top =0.91+0.009*n       # the top of the subplots of the figure
            wspace = 0              # the amount of width reserved for blank space between subplots
            hspace = 1              # the amount of height reserved for white space between subplots
            print("plot data : n size " + str(n) + " x len " + str(len(x)) + " y len " + str(len(y)) )
            figure.subplots_adjust(left =left , bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
            if(n > 1):
                for i in range(0,n): 
                    if( shift != 0 and name[i].find('Interference')!= -1):
                        axs[i].plot(( x[0:len(y[i])-shift]),(y[i][0:len(y[i])-shift]),color[i%len(color)])
                    else:
                        axs[i].plot(x[0:len(y[i])],y[i],color[i%len(color)])
                   # print(str(name[i]) + " " + str(len(y[i])) +" " +  str(len(x[0:len(y[i])])))
                    axs[i].yaxis.major.formatter._useMathText = True
                    axs[i].set_title(name[i])#, fontsize = 8)
                    axs[i].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
                    #axs[i].ticklabel_format(style='plain')
                    #axs[i].ticklabel_format(useMathText = True)
                    #axs[i].ticklabel_format(useOffset=False)
                    #serif #math
                    if ( i == n-1):
                        #axs[i].set(xlabel= r'Time, ns $\frac{T_i}{2}$', ylabel=leg[i])
                        axs[i].set(xlabel= r'Time, $ns$', ylabel=leg[i])
            elif(len(y) == 1):
                for i in range(0,n): 
                    axs.plot(x[0:len(y[i])-shift],y[i][0:len(y[i])-shift],color[i%len(color)])
                    axs.set_title(name[i])#,fontsize=8)
                    axs.yaxis.major.formatter._useMathText = True
                    axs.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
                    axs.set(xlabel=r'Time, $ns$' , ylabel=leg[i])#, fontsize = 8)
            else:
                 
                 axs.plot(x,y,"go",markersize = 3)
                 axs.set_title(name)
                 axs.grid(which = 'both' , linewidth =1)

                 axs.set(xlabel=r'Injection Ratio, $dB$', ylabel=leg)
            if (Ip != None and L != None):
                Iarr = Ip*np.ones(L)
                axs[0].plot(x,Iarr*1e12, 'b--' )
            if (Ib != None and L != None):
                Iarrb = Ib*np.ones(L)
                axs[0].plot(x,Iarrb*1e12, 'g--')          
            figure.set_facecolor('w')        
            
            if (n > 7):
                n = 7
            canvas.get_tk_widget().place(x=10,y=10,width= W ,height= H)
            toolbar.place(x = 10 , y = H + 30)          
            canvas.draw()
            root1.update()

        def getplot(self,x,y,name = None, xl = None, yl = None,lim = None):
            plt.plot(x,y)
            plt.xlabel(xl  )
            plt.ylabel(yl)
            plt.title(name  )

            plt.show()

            
        def getplot_n(self,n,x,y,name,leg,L = None,Ib = None,Ip = None):
            color = ["tab:orange","tab:green","tab:red","tab:blue","tab:red"]
            fig ,axs = plt.subplots(n,1)
            for i in range(0,n):

                axs[i].plot(x,y[i],color[i%n])
                axs[i].set_title(name[i])
                axs[0].set(xlabel='T ns', ylabel=leg)
                
            if (Ip != None and L != None):
                Iarr = Ip*np.ones(L)
                axs[1].plot(x,Iarr*1e12, 'b--' )
            if (Ib != None and L != None):
                Iarrb = Ib*np.ones(L)
                axs[1].plot(x,Iarrb*1e12, 'g--')    
   
            plt.show()