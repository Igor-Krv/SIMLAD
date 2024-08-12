import numpy as np
from scipy import optimize
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
e = 1.6e-19        #The electron charge, C
h = 6.63e-34       #The Planck constant, J*s

nu0 = 1.93548e14
#@jitclass()

class Laser:
            def __init__(self, Ip,Driver,arg = None, master = None):
                self.Ib = Driver.Ib
                self.Ip = Driver.Ip
                self.eta = arg[0]     
                self.time = Driver.time
                self.tau_ph =arg[1]                     #The photon lifetime, ns 
                self.Ntr = arg[2]*1e7                   #The transparency carrier number
                self.Nth = arg[3]*1e7                   #The threshold carrier number
                self.chi = arg[4]                       #The gain compression factor, 1/W
                self.Csp = arg[5]                       #The fraction of spontaneously emitted photons that end up in the active mode
                self.Gamma = arg[6]                     #The confinement factor
                self.tau_e = arg[7]                     #The effective carrier lifetime, ns
                self.alpha = arg[8]                     #The linewidth enhancement factor                
                self.L = Driver.L
                self.Ith = e * self.Nth/self.tau_e      #Threshold current
                self.Itr = e * self.Ntr/self.tau_e      #Transparency current
                #print("Threshold current :" , self.Ith)
                self.nPointsPerPulse = Driver.nPointsPerPulse
                #print("N poits per pulse :" , self.nPointsPerPulse)
                #print("Array size" , self.L)
                self.hbar = 1.054 *1e-25                #J*ns
                self.nu0 = nu0
                self.Dt = Driver.Dt
                self.Master = master
                if (self.Master):
                    self.Qst = self.Gamma * self.tau_ph/e 
                
                self.N = np.zeros(self.L)               #carrier number array
                self.Q = np.zeros(self.L)               #laser intensity
                self.phi = np.zeros(self.L)             #phase 
                self.nu = np.zeros(self.L)              #frequency
                self.DT = np.zeros(self.L)              #temperature difference
                
                self.phiS = np.zeros(self.L)
                self.DT1 = np.zeros(self.L)
                self.DT2 = np.zeros(self.L)
                self.DTJ = np.zeros(self.L)
                self.DTS = np.zeros(self.L)
                
                self.r1 = arg[18]
                self.r2 = arg[19]

                self.rh_add = arg[20] * 1e9 
                self.k12 = arg[21]
                self.tauj = arg[22]
                self.lambdas = arg[23]
                self.Cj = arg[24]
                self.RL = 8 * 1e9                       #resistance
                self.E = np.zeros(self.L,dtype = complex)
              
                self.Nav = np.zeros(self.L)             #array for average solutions without noise
                self.Qav = np.zeros(self.L)
                self.phiav = np.zeros(self.L)
                
                self.S = np.zeros(self.L)               #interference
                self.teta1 = 0                          #params for phase in interference
                self.teta2 = 0
                
                self.N[0] = self.Nth                    #initial values
                self.Q[0] = 0
                #!TODO
                self.DT[0] = 0
                self.S[0] = 0
                
                #temperature effects
                self.rh = arg[9]*1e9                    #thermal resistance , K/W
                self.th = arg[10]                       #thermal rise time , ns
                self.lambd = arg[11]                    #wavelength , micrometers
                self.kap_om = arg[12]#2 * np.pi * 10    #temperature coefficient of frequency , GHZ                
                #injection parameters
                self.t1 = arg[13]                       # Field transmission coefficient of the facet through which the external optical field is injected into the laser cavity
                self.w_inj = 2*np.pi*arg[14]            #2*np.pi*(self.w_M - self.om_th)    use it instead of w_M and om_th
                self.om_th = arg[15] 
                self.Ll = arg[16]                       # length of the resonator , micrometers
                self.kap_inj =  arg[17]                 # coupling coefficient.
               
                self.delay = arg[25]                    #interferometer delay
                self.phase_m1 = arg[26]                 #phase modulation 1
                self.phase_m2 = arg[27]
                #constant values for optimization 
                self.v1 = 2/self.tau_e
                self.v2 = 1e9*self.chi*self.eta*h*nu0/2/self.Gamma/self.tau_ph #chi*eta*h*nu0/2/Gamma/tau_ph = chi_Q is the dimensionless compression factor
                self.gammaQ = self.chi*self.eta*h*nu0/2/self.Gamma/self.tau_ph
                
                self.v3 = self.tau_e*1/self.Dt 
                self.v4 = self.tau_ph*(1/self.Dt)
                self.v5 = (self.Dt*2*np.pi)
                self.c1 = self.Gamma/(1/self.tau_ph)
                self.c2 = self.Nth - self.Ntr
                self.c3 = 1e9*self.chi * self.eta * h * nu0/self.Gamma/self.tau_ph
                self.c4 = self.alpha/2/self.tau_ph  
        
            #unstable, turned off in Main core
            def sci_ODE_slave(self,I,EM,phase = None,freq = None, noise = None, temp = None, interference = None, biexp = None, field = None):
                #create local values
                Csp = self.Csp       
                gammaQ = self.gammaQ
                time = self.time
                tau_e = self.tau_e
                Ntr = self.Ntr
                Nth = self.Nth
                tau_ph = self.tau_ph
                Gamma = self.Gamma
                alpha = self.alpha
                c4 = alpha/2/tau_ph  
                chi = self.chi
                eta = self.eta  
                th = self.th
                rh = self.rh
                Ith = self.Ith 
                RL = 8 * 1e9
                kap_om = self.kap_om
                kap_inj = self.kap_inj
                w_inj = self.w_inj
                         
                def rate(t,arr,dI,dQM,dphiM):
                       
                    #N = arr[0]
                    #Q = arr[1]                   
                    N , Q, DT, phi = arr
                
                    if Q < 0:
                        Q = 1e-12                
                
                    dDT = 0
                    if(temp == 1):
                        dDT =  - DT/ th  + 0.8 * rh/th * (1 - eta) * (dI(t) - Ith)  + RL*rh/th*((dI(t))**2 - Ith**2)
                
                    dN =   dI(t)/e  -N/(tau_e)  - Q/Gamma/(tau_ph)*( N -  Ntr)/( Nth - Ntr) \
                        / np.sqrt(1 + 2*1e9*chi*eta*h*nu0/Gamma/tau_ph*abs(Q))
                         
                    dQ =  (( N - Ntr)/(Nth - Ntr)/ np.sqrt(1 + 2*1e9*chi*eta*h*nu0/Gamma/tau_ph*abs(Q) ) - 1)*Q/tau_ph + \
                          Csp* N/tau_e +\
                          kap_inj*np.sqrt((Q) * dQM(t))*np.cos(w_inj*t + dphiM(t) - phi ) 
                          
                    dphi =  alpha/2/tau_ph*((N - Ntr)/(Nth - Ntr)-1) - kap_om*dDT +\
                        kap_inj* np.sqrt(dQM(t)/(Q))*np.sin(w_inj*t + dphiM(t) - phi) 

                    dydt = [dN, dQ,dDT,dphi]                
                    return dydt
                
                
                N0 = self.Nth
                Q0 = np.float64(1)
                phi0 = np.float64(0)
                DT0 = 0
                
                phiM = np.zeros(self.L)
                QM = np.zeros(self.L)
                for i in range(self.L):
                    QM[i] = abs(EM[i])**2
                for i in range(self.L-1):
                    phiM[i+1] = phiM[i] +  self.Dt/abs(EM[i])**2 *  ( (EM[i+1] - EM[i])/self.Dt * np.conjugate(EM[i])).imag
                                 
                Interpolated = interp1d(self.time,I,fill_value = "extrapolate")
                QM_in = interp1d(self.time,QM,fill_value = "extrapolate")
                phiM_in = interp1d(self.time,phiM,fill_value = "extrapolate")
                
                print(f"t_eval time {len(time)}")
                print(f"Injection {len(QM)}")  
                output = solve_ivp(rate,[0,time[-1]],[N0,Q0,DT0,phi0], method = "RK45" ,t_eval = time,args = ( Interpolated,QM_in,phiM_in,), rtol=1e-4, atol = 1e-6)
            
                E = np.sqrt(output.y[1]) * np.exp(1j*output.y[2])
               
                return [output.y[0],E ,output.y[2],output.y[3],output.t]
            
            def sci_ODE(self,I,phase = None,freq = None, noise = None, temp = None, interference = None, field = None, biexp = None):

                Csp = self.Csp       

                gammaQ = self.gammaQ
                time = self.time
                feed = np.zeros(0)
                count  = 0 
                tau_e = 1
                Ntr = self.Ntr
                Nth = self.Nth
                tau_ph = self.tau_ph
                Gamma = self.Gamma
                alpha = self.alpha
                c4 = alpha/2/tau_ph  
                chi = self.chi
                eta = self.eta
                th = self.th
                rh = self.rh
                Ith = self.Ith 
                RL = 8 * 1e9
                kap_om = self.kap_om
                print("temp is" + str(temp))
                                   
                def rate(t,arr,dI):
                       
                    #N = arr[0]
                    #Q = arr[1]                   
                    N , Q, DT, phi = arr
                    dDT = 0
                    if(temp == 1):
                        dDT =  - DT/ th  + 0.8 * rh/th * (1 - eta) * (dI(t) - Ith)  + RL*rh/th*((dI(t))**2 - Ith**2)
                
                    dN =   dI(t)/e  -N/(tau_e)  - Q/Gamma/(tau_ph)*( N -  Ntr)/( Nth - Ntr) \
                        / np.sqrt(1 + 1e9*chi*eta*h*nu0/Gamma/tau_ph*abs(Q))
                    dQ =  (( N - Ntr)/(Nth - Ntr)/ np.sqrt(1 + 1e9*chi*eta*h*nu0/Gamma/tau_ph*abs(Q) ) - 1)*Q/tau_ph + \
                          Csp* N/tau_e 
                    dphi =  alpha/2/tau_ph*((N - Ntr)/(Nth - Ntr)-1) - kap_om*dDT

                    dydt = [dN, dQ,dDT,dphi]                
                    return dydt
                
                N0 = self.Nth
                Q0 = np.float64(1)
                phi0 = np.float64(0)
                DT0 = 0
                 
                Interpolated = interp1d(self.time,I,fill_value = "extrapolate")
                
                output = solve_ivp(rate,[0,self.time[-1]],[N0,Q0,DT0,phi0], method = "RK45" ,t_eval = time,args = ( Interpolated,), rtol=1e-4, atol = 1e-6)
                E = np.sqrt(output.y[1]) * np.exp(1j*output.y[2]) #electric field
                
                return [output.y[0],E,output.y[2],output.y[3],output.t]    
            
            def runge_kutta_N(self, In, Nn, Qn):
                return   In/e*self.Dt - Nn/self.v3 - \
                        Qn/self.c1*(Nn - self.Ntr)/(self.c2) \
                        / np.sqrt((1+2*self.gammaQ*1e9*Qn))*self.Dt 
                        
            def runge_kutta_Q(self, Nn, Qn):  
                return  ((Nn- self.Ntr)/(self.c2) \
                        / np.sqrt((1 + 2*self.gammaQ*1e9*Qn)) - 1)*Qn/self.v4 + \
                        self.Csp*Nn/self.v3  
                                   
            def runge_kutta_phi(self, phin, Nn,DTn):
                return  ( (self.c4*((Nn - self.Ntr)/(self.c2) - 1) -  self.kap_om *DTn)*self.Dt )
             
            def runge_kutta_QS(self, QM, phiM, Nn, Qn, phin, i, injfix):                         
                return  ((Nn- self.Ntr)/(self.c2) \
                        / np.sqrt((1 + 2*self.gammaQ*1e9*Qn)) - 1)*Qn/self.v4 + self.Csp*Nn/self.v3   + \
                          + injfix*self.kap_inj*np.sqrt(abs(Qn * QM))*np.cos(self.w_inj*(i+1)*self.Dt + phiM - phin ) *self.Dt
                          
            def runge_kutta_phiS(self,phin,Nn,Qn,DTn,QM,phiM,i ):
                return    (self.c4*((Nn - self.Ntr)/(self.c2) - 1)-  self.kap_om *self.DT[i] )*self.Dt    +\
                            self.kap_inj* np.sqrt(QM/Qn)\
                            *np.sin(self.w_inj*(i+1)*self.Dt + phiM - phin)*self.Dt   
                
            def runge_kutta(self,I,phase = None, freq =None, noise = None, temp = None, interference = None, biexp = None):
                self.N = np.zeros(self.L)
                self.Q = np.zeros(self.L)
                self.phi = np.zeros(self.L)
                self.nu = np.zeros(self.L)
                self.DT = np.zeros(self.L)
                
                self.Nav[0] = self.Nth
                self.Qav[0] = 1
                self.N[0] = self.Nth            #initial values
                self.Q[0] = 1
                self.DT[0] = 0
                self.DT1[0] = 0
                self.DT2[0] = 0
                self.DTJ[0] = 0
                self.DTS[0] = 0
                
                self.S[0] = 0
                self.E[0] = 1
        
                if (temp == True):
                    print("1.24/self.lambd :" + str(1.24/self.lambd))
                    for i in range(0, self.L-1):
                        if (biexp == False):
                            self.DT[i+1] = self.DT[i] - self.DT[i]/ self.th * self.Dt + 1.24/self.lambd * self.rh/self.th * (1 -self.eta) * (I[i] - self.Ith) * self.Dt + self.RL*self.rh/self.th*(I[i]**2 - self.Ith**2)*self.Dt
                        else:
                            self.DT1[i+1] = self.DT1[i] +( -(self.DT1[i] - self.rh*1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.r1 - self.k12*(self.DT1[i] -self.DT2[i]) )*self.Dt
                            self.DT2[i+1] = self.DT2[i] +( -(self.DT2[i] - self.rh*1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.r2 - self.k12*(self.DT2[i] -self.DT1[i]) )*self.Dt
                            self.DTJ[i+1] = self.DTJ[i] +( -self.DTJ[i]/self.tauj + self.rh*self.Cj/self.tauj * 1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.Dt
                for i in range(0, self.L - 1):
                    k1N = self.runge_kutta_N(I[i], self.N[i], self.Q[i]) 
                    k2N = self.runge_kutta_N(I[i] , self.N[i] + k1N/2, self.Q[i] )
                    k3N = self.runge_kutta_N(I[i] , self.N[i] + k2N/2, self.Q[i]  )
                    k4N = self.runge_kutta_N(I[i]  , self.N[i] + k3N, self.Q[i] )
                    self.N[i+1] = self.N[i] + 1/6*(k1N + 2*k2N + 2*k3N + k4N)
                    
                    k1Q = self.runge_kutta_Q(self.N[i], self.Q[i])
                    k2Q = self.runge_kutta_Q(self.N[i]  , self.Q[i] + k1Q/2)
                    k3Q = self.runge_kutta_Q(self.N[i] , self.Q[i] + k2Q/2)
                    k4Q = self.runge_kutta_Q(self.N[i] , self.Q[i] + k3Q)
                    
                    self.Q[i+1] = self.Q[i] + 1/6*(k1Q + 2*k2Q + 2*k3Q + k4Q)
                    
                    if(self.N[i+1] < 0): self.N[i+1] = 1e-12
                    if(self.Q[i+1] < 0): self.Q[i+1] = 1e-3
                    
                    k1phi = self.runge_kutta_phi(self.phi[i], self.N[i], self.DT[i])
                    k2phi = self.runge_kutta_phi(self.phi[i] + k1phi/2, self.N[i], self.DT[i])
                    k3phi = self.runge_kutta_phi(self.phi[i] + k2phi/2, self.N[i], self.DT[i])
                    k4phi = self.runge_kutta_phi(self.phi[i] + k3phi/2, self.N[i], self.DT[i])
                    
                    self.phi[i+1] = self.phi[i] + 1/6*(k1phi + 2*k2phi + 2*k3phi + k4phi)
                
                #if(freq):
                for i in range(0,self.L-1):
                    self.nu[i] = (self.phi[i+1] - self.phi[i]) /self.v5
                self.nu[self.L-1] = self.nu[self.L-2]

                self.E = np.sqrt(self.Q) * np.exp(1j*self.phi)
                return [self.N, self.E,self.DT,self.phi,self.nu]      

            def runge_kutta_slave(self,EM,I,injfix, phase = None, freq =None, noise = None, temp = None, interference = None, biexp = None):
                self.N = np.zeros(self.L)
                self.Q = np.zeros(self.L)
                self.phi = np.zeros(self.L)
                self.nu = np.zeros(self.L)
                self.DT = np.zeros(self.L)
                
                self.Nav[0] = self.Nth
                self.Qav[0] = 1
                self.N[0] = self.Nth            #initial values
                self.Q[0] = 1
                self.DT[0] = 0
                self.DT1[0] = 0
                self.DT2[0] = 0
                self.DTJ[0] = 0
                self.DTS[0] = 0
                
                self.S[0] = 0
                self.E[0] = 1
                if (injfix): inj_m = 2
                else: inj_m = 1 # multiply by 2 in Q equation
                phiM = np.zeros(self.L)
                QM = np.zeros(self.L)
                for i in range(self.L):
                    QM[i] = abs(EM[i])**2
                for i in range(self.L-1):
                    phiM[i+1] = phiM[i] +  self.Dt/abs(EM[i])**2 *  ( (EM[i+1] - EM[i])/self.Dt * np.conjugate(EM[i])).imag
                if (temp):
                    for i in range(0, self.L - 1):               
                        self.DT[i+1] = self.DT[i] - self.DT[i]/ self.th * self.Dt + 0.8 * self.rh/self.th * (1 -self.eta) * (I[i] - self.Ith) * self.Dt + self.RL*self.rh/self.th*(I[i]**2 - self.Ith**2)*self.Dt

                if (temp == True):
                    for i in range(0, self.L-1):
                        if (biexp == False):
                            self.DT[i+1] = self.DT[i] - self.DT[i]/ self.th * self.Dt + 0.8 * self.rh/self.th * (1 -self.eta) * (I[i] - self.Ith) * self.Dt + self.RL*self.rh/self.th*(I[i]**2 - self.Ith**2)*self.Dt
                        else:
                            self.DT1[i+1] = self.DT1[i] +( -(self.DT1[i] - self.rh*1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.r1 - self.k12*(self.DT1[i] -self.DT2[i]) )*self.Dt
                            self.DT2[i+1] = self.DT2[i] +( -(self.DT2[i] - self.rh*1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.r2 - self.k12*(self.DT2[i] -self.DT1[i]) )*self.Dt
                            self.DTJ[i+1] = self.DTJ[i] +( -self.DTJ[i]/self.tauj + self.rh*self.Cj/self.tauj * 1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.Dt
                            
                for i in range(0, self.L - 1):
                    k1N = self.runge_kutta_N(I[i], self.N[i], self.Q[i]) 
                    k2N = self.runge_kutta_N(I[i] , self.N[i] + k1N/2, self.Q[i] )
                    k3N = self.runge_kutta_N(I[i] , self.N[i] + k2N/2, self.Q[i]  )
                    k4N = self.runge_kutta_N(I[i]  , self.N[i] + k3N, self.Q[i] )
                    self.N[i+1] = self.N[i] + 1/6*(k1N + 2*k2N + 2*k3N + k4N)
                    
                    k1Q = self.runge_kutta_QS(QM[i], phiM[i], self.N[i], self.Q[i], self.phi[i], i, inj_m)
                    k2Q = self.runge_kutta_QS(QM[i], phiM[i], self.N[i]  , self.Q[i] + k1Q/2, self.phi[i], i, inj_m)
                    k3Q = self.runge_kutta_QS(QM[i], phiM[i], self.N[i] , self.Q[i] + k2Q/2, self.phi[i], i, inj_m)
                    k4Q = self.runge_kutta_QS(QM[i], phiM[i], self.N[i] , self.Q[i] + k3Q, self.phi[i], i, inj_m)
                    
                    self.Q[i+1] = self.Q[i] + 1/6*(k1Q + 2*k2Q + 2*k3Q + k4Q)
                    
                    # if(self.N[i+1] < 0): self.N[i+1] = 1e-12
                    # if(self.Q[i+1] < 0): self.Q[i+1] = 1e-3
                    
                    if(self.N[i+1] < 0): self.N[i+1] = 1e-12
                    if(self.Q[i+1] < 0): self.Q[i+1] = self.Qav[i+1]
                        
                    k1phi = self.runge_kutta_phiS(self.phi[i], self.N[i],self.Q[i], self.DT[i],QM[i],phiM[i],i)
                    k2phi = self.runge_kutta_phiS(self.phi[i] + k1phi/2, self.N[i],self.Q[i], self.DT[i],QM[i],phiM[i],i)
                    k3phi = self.runge_kutta_phiS(self.phi[i] + k2phi/2, self.N[i],self.Q[i], self.DT[i],QM[i],phiM[i],i)
                    k4phi = self.runge_kutta_phiS(self.phi[i] + k3phi/2, self.N[i],self.Q[i], self.DT[i],QM[i],phiM[i],i)
                    self.phi[i+1] = self.phi[i] + 1/6*(k1phi + 2*k2phi + 2*k3phi + k4phi)
                                           
                if(freq):
                    for i in range(0,self.L-1):
                        self.nu[i] = (self.phi[i+1] - self.phi[i]) /self.v5
                    self.nu[self.L-1] = self.nu[self.L-2]

                self.E = np.sqrt(self.Q) * np.exp(1j*self.phi)
                return [self.N, self.E,self.DT,self.phi,self.nu] 
                    
            def euler_m(self,I, phase = None,freq = None, noise = None, temp = None, interference = None, field = None, biexp = None):
                # if(interference):
                #     print ('interference enabled') 
                # if(temp):
                #     print('temperature effects enabled')
                # if(noise):
                #     print('stochastics enabled')
                # if(field):
                #     print('Field calculation enabled')   
                # if(biexp):
                #     print('biexponential tempertature') 
                
                self.N = np.zeros(self.L)
                self.Q = np.zeros(self.L)
                self.phi = np.zeros(self.L)
                self.nu = np.zeros(self.L)
                
                self.Nav[0] = self.Nth
                self.Qav[0] = 1
                self.N[0] = self.Nth            #initial values
                self.Q[0] = 1
                self.DT[0] = 0
                self.DT1[0] = 0
                self.DT2[0] = 0
                self.DTJ[0] = 0
                self.DTS[0] = 0
                
                self.S[0] = 0
                self.E[0] = 1
                if ( field == False):
                    self.N[0] = self.Nav[0]
                    self.Q[0] = self.Qav[0]
                    self.phi[0] = self.phiav[0]
                    for i in range(0, self.L - 1):
                        self.Nav[i+1] = self.Nav[i] + I[i]/e*self.Dt - self.Nav[i]/self.v3 - self.Qav[i]/self.Gamma/self.tau_ph*(self.Nav[i]-self.Ntr) \
                        /(self.c2) / np.sqrt(1 + self.c3*self.Qav[i])*self.Dt
    
                        self.Qav[i+1] = self.Qav[i] + ((self.Nav[i] - self.Ntr)/(self.c2)/ np.sqrt(1+ self.c3 \
                        * self.Qav[i]) - 1) * self.Qav[i] / self.tau_ph * self.Dt + self.Csp*self.Nav[i] / self.v3
    
                        self.phiav[i+1] = self.phiav[i] + self.c4 * ((self.Nav[i] - self.Ntr)/(self.c2)-1) * self.Dt
            
                self.xi_N = np.random.randn(self.L)
                self.xi_Q = np.random.randn(self.L)
                self.xi_phi = np.random.randn(self.L)
                self.xi_1 = np.random.randn(self.L)
                self.xi_2 = np.random.randn(self.L)
                self.xi_3 = np.random.randn(self.L)
                self.Q_E = np.zeros(self.L)
                #datawrite = []
               
                if (temp == True):
                    for i in range(0, self.L-1):
                        if (biexp == False):
                            self.DT[i+1] = self.DT[i] - self.DT[i]/ self.th * self.Dt + 0.8 * self.rh/self.th * (1 -self.eta) * (I[i] - self.Ith) * self.Dt + self.RL*self.rh/self.th*(I[i]**2 - self.Ith**2)*self.Dt
                        else:
                            self.DT1[i+1] = self.DT1[i] +( -(self.DT1[i] - self.rh*1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.r1 - self.k12*(self.DT1[i] -self.DT2[i]) )*self.Dt
                            self.DT2[i+1] = self.DT2[i] +( -(self.DT2[i] - self.rh*1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.r2 - self.k12*(self.DT2[i] -self.DT1[i]) )*self.Dt
                            self.DTJ[i+1] = self.DTJ[i] +( -self.DTJ[i]/self.tauj + self.rh*self.Cj/self.tauj * 1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.Dt
                    if (biexp == True):
                      for i in range(0, self.L):
                          self.DT[i] = (self.DT1[i] + self.DT2[i])/2 + self.DTJ[i]
                for i in range(0, self.L-1):
                    if(noise):
                        if (self.Nav[i] <= 0):
                                self.Nav[i] = 1e-12
                        if (self.Qav[i] <= 0):
                                self.Qav[i] = 1e-12
                        noise_N = -2*np.sqrt(self.Csp*self.Nav[i]*self.Qav[i]/2/self.tau_e) * (np.cos(self.phiav[i]) * self.xi_Q[i] + np.sin(self.phiav[i])*self.xi_phi[i]) \
                        * np.sqrt(self.Dt) + np.sqrt(2*self.Nav[i]/self.tau_e)*self.xi_N[i]*np.sqrt(self.Dt)
                        
                        noise_Q = 2*np.sqrt( self.Csp*self.Nav[i] *self.Qav[i]/2/self.tau_e)*(self.xi_Q[i]*np.cos(self.phiav[i]) + self.xi_phi[i]*np.sin(self.phiav[i]) )*np.sqrt(self.Dt) 
                        
                        noise_phi = np.sqrt(self.Csp * self.Nav[i]/2/self.tau_e/self.Qav[i])*(self.xi_phi[i]*np.cos(self.phiav[i]) - self.xi_Q[i]*np.sin(self.phiav[i]))*np.sqrt(self.Dt)
                        
                    else:
                        noise_N = 0
                        noise_Q = 0
                        noise_phi = 0                        
                    
                    if (field == False):
                                               
                         # counter = counter + 1
                         # N[i+1] = N[i] + I[i]/e* Dt -  N[i]/ v3 - \
                         # Q[i]/ c1*( N[i] -  Ntr)/( c2) \
                         # / np.sqrt((1+2* gammaQ*1e9* Q[i]))* Dt + noise_N
                        
                         # Q[i+1] =  Q[i] + (( N[i]-  Ntr)/( c2) \
                         # / np.sqrt((1 + 2* gammaQ*1e9* Q[i])) - 1)* Q[i]/ v4 + \
                         # Csp* N[i]/ v3 + noise_Q  
                        
                         # if( N[i+1] < 0):  N[i+1] = 1e-3
                         # if( Q[i+1] < 0):  Q[i+1] = Qav[i+1]
                         
                        
                         # if(  temp == False):
                         #     phi[i+1] =  phi[i] +  c4*(( N[i] -  Ntr)/( c2) - 1)* Dt + noise_phi
                         # elif ( temp == True):
                         #    phi[i+1] =  phi[i] + ( c4*(( N[i] -  Ntr)/( c2) - 1) -   kap_om * DT[i])* Dt + noise_phi
                               
                        self.N[i+1] = self.N[i] + I[i]/e*self.Dt - self.N[i]/self.v3 - \
                        self.Q[i]/self.c1*(self.N[i] - self.Ntr)/(self.c2) \
                        / np.sqrt((1+2*self.gammaQ*1e9*self.Q[i]))*self.Dt + noise_N
                        
                        self.Q[i+1] = self.Q[i] + ((self.N[i]- self.Ntr)/(self.c2) \
                        / np.sqrt((1 + 2*self.gammaQ*1e9*self.Q[i])) - 1)*self.Q[i]/self.v4 + \
                        self.Csp*self.N[i]/self.v3 + noise_Q  
                        
                        if(self.N[i+1] < 0): self.N[i+1] = 1e-3
                        if(self.Q[i+1] < 0): self.Q[i+1] = self.Qav[i+1]
                        
                        #self.phi[i+1] = self.phi[i] + (self.c4*((self.N[i] - self.Ntr)/(self.c2) - 1) -  self.kap_om *self.DT[i])*self.Dt + noise_phi                   
                        if(  temp == False):
                            self.phi[i+1] = self.phi[i] + self.c4*((self.N[i] - self.Ntr)/(self.c2) - 1)*self.Dt + noise_phi
                        elif ( temp == True):
                            self.phi[i+1] = self.phi[i] + (self.c4*((self.N[i] - self.Ntr)/(self.c2) - 1) -  self.kap_om *self.DT[i])*self.Dt + noise_phi
                    # if (field):
                    #     #if(i != 0): 
                    #     self.Q[i] = abs(self.E[i])**2
                    #     #Q_E = abs(self.E[i])**2
                    #     self.N[i+1] = self.N[i] + I[i]/e*self.Dt - self.N[i]/self.v3 - \
                    #     self.Q[i]/self.c1*(self.N[i] - self.Ntr)/(self.c2)* \
                    #     self.Dt/np.sqrt((1 + 2*self.v2*abs(self.E[i])**2)) -(np.sqrt(2*self.Csp*self.N[i]*abs(self.E[i])**2/self.tau_e)*self.xi_1[i]-np.sqrt(2*self.N[i]/self.tau_e)*self.xi_3[i])*np.sqrt(self.Dt)
                    #     if (temp == True):
                    #         self.E[i+1] = self.E[i] + self.E[i] * self.Dt / self.tau_ph/2 *( ((self.N[i] - self.Ntr)/self.c2)/np.sqrt((1 + 2*self.v2*abs(self.E[i])**2 )) - 1 +  1j*self.alpha*((self.N[i] - self.Ntr)/self.c2-1) \
                    #         -2*1j*self.tau_ph*self.kap_om*self.DT[i]) +\
                    #         np.sqrt( self.Csp * self.N[i] /2/self.tau_e )* (self.xi_1[i] *np.sqrt(self.Dt) + self.xi_2[i]*np.sqrt(self.Dt)*1j )
                    #     else:
                    #         self.E[i+1] = self.E[i] + self.E[i]*self.Dt/self.tau_ph/2 * ( (self.N[i] - self.Ntr)/self.c2/np.sqrt((1 + 2*self.v2*abs(self.E[i])**2)) - 1 +  1j*self.alpha*( (self.N[i] - self.Ntr)/self.c2-1) ) \
                    #         + np.sqrt( self.Csp * self.N[i] /2/self.tau_e )* (self.xi_1[i]**2 *np.sqrt(self.Dt) + self.xi_2[i]**2*np.sqrt(self.Dt)*1j )
                    #     if( phase == True):
                    #         self.phi[i+1] = self.phi[i] + self.Dt/abs(self.E[i])**2 *  ( (self.E[i+1] - self.E[i])/self.Dt * np.conjugate(self.E[i]) ).imag
                        
               # if (temp == True and phase == True):

                  #  for i in range(0,self.L-1):
                      #  self.phi[i+1] = self.phi[i] + (self.c4*((self.N[i] - self.Ntr)/(self.c2) - 1) -  self.kap_om *self.DT[i])*self.Dt + noise_phi
                        
                       # self.phiS[i+1] = self.phiS[i] + (self.c4*((self.N[i] - self.Ntr)/(self.c2) - 1) -  self.kap_om *self.DTS[i])*self.Dt + noise_phi
                        #datawrite.append(noise_Q)
                 
                #np.savetxt('addM.txt', datawrite) 
                if(freq):
                    for i in range(0,self.L-1):
                        self.nu[i] = (self.phi[i+1] - self.phi[i]) /self.Dt/2/np.pi
                    self.nu[self.L-1] = self.nu[self.L-2]
     
                if (interference == True):
                    #0.8 sec - len(8000)
                    # delay / Dt instead of nppp
                    shift = int(self.delay / self.Dt)
                    #print('master interferometer shift '+ str(shift))
                    #nppp = int(self.nPointsPerPulse)
                    Q1 = self.Q[0:self.L-shift]
                    Q2 = self.Q[shift:self.L]
                    phi1 = self.phi[0:self.L-shift]
                    phi2 = self.phi[shift:self.L]
                    for i in range(0,self.L-shift):
                        self.S[i] = abs( np.sqrt(Q1[i]) * np.exp(1j * phi1[i]+ 1j*self.phase_m1) + np.sqrt(Q2[i]) * np.exp(1j * phi2[i] + 1j*self.phase_m2) )**2
                # self.P = np.zeros(self.L)
                # self.P = self.eta*h*(2*np.pi*3*10**8/self.lambd )/(2*self.Gamma*self.tau_ph)*self.Q
                # if (power_out == 0):
                self.E = np.sqrt(self.Q) * np.exp(1j*self.phi)
                # else:
                    # self.E = np.sqrt(self.P) * np.exp(1j*self.phi)
                Data =  [self.N, self.E,self.DT,self.phi,self.nu]
                return Data
            
            def euler_m_slave(self,EM,driver,I,injfix,noise = None,temp = None,interference = None, field = None ,injection = None,phase = None,freq = None,biexp =None):   
                print("Optical injection enabled")
                if(interference): print ('interference enabled ')
                if(temp): print('temperature effects enabled ')
                if(noise): print('stochastics enabled')
                if(field): print('Field calculation enabled') 
                
                
                #self.kap_inj =  arg[17]#33.095             # coupling coefficient.
                #self.w_inj = 2*np.pi*(self.w_M - self.om_th)
                #print ("INJECTION PARAMETERS: kappa_inj " + str(self.kap_inj) + " w_inj " + str(self.w_inj))
                #print ("INJECTION PARAMETERS: w_M " + str(self.w_M) + " om_th " + str(self.om_th))
                #print("driver time is" + str(len(driver.time)))
                #print("QM " + str(len(laserM.Q)))
                self.N = np.zeros(self.L)
                self.Q = np.zeros(self.L)
                self.phi = np.zeros(self.L)
                self.nu = np.zeros(self.L)
                self.S = np.zeros(self.L)
                self.E = np.zeros(self.L)

                self.N[0] = self.Nth            #initial values
                self.Q[0] = 1
                self.DT[0] = 0
                self.S[0] = 0
                self.E[0] = 1
                
                self.xi_A = np.random.randn(self.L)
                self.xi_B = np.random.randn(self.L)
                self.xi_C = np.random.randn(self.L)
                self.xi_1 = np.random.randn(self.L)
                self.xi_2 = np.random.randn(self.L)
                noise_N = np.zeros(self.L)
                noise_Q = np.zeros(self.L)
                noise_phi = np.zeros(self.L) 
                
                self.Nav[0] = self.Nth
                self.Qav[0] = 1
                self.N[0] = self.Nth            #initial values
                self.Q[0] = 1
                self.DT[0] = 0
                self.DT1[0] = 0
                self.DT2[0] = 0
                self.DTJ[0] = 0
                self.DTS[0] = 0
                
                self.S[0] = 0
                self.E[0] = 1
                if ( field == False):
                    self.N[0] = self.Nav[0]
                    self.Q[0] = self.Qav[0]
                    self.phi[0] = self.phiav[0]
                    for i in range(0, self.L - 1): 
                        if(noise):
                            self.Nav[i+1] = self.Nav[i] + I[i]/e*self.Dt - self.Nav[i]/self.v3 - \
                            self.Qav[i]/self.Gamma/self.tau_ph*(self.Nav[i]-self.Ntr) /(self.c2) / np.sqrt(1 + self.c3*self.Qav[i])*self.Dt
        
                            self.Qav[i+1] = self.Qav[i] + ((self.Nav[i] - self.Ntr)/(self.c2)/ np.sqrt(1+ self.c3 \
                            * self.Qav[i]) - 1) * self.Qav[i] / self.tau_ph * self.Dt + self.Csp*self.Nav[i] / self.v3
    
                            self.phiav[i+1] = self.phiav[i] + self.c4 * ((self.Nav[i] - self.Ntr)/(self.c2)-1) * self.Dt 
                        if (self.Nav[i] <= 0):
                                self.Nav[i] = 1e-12
                        if (self.Qav[i] <= 0):
                                self.Qav[i] = 1e-3
                        noise_N[i] = -2*np.sqrt(self.Csp*self.Nav[i]*self.Qav[i]/2/self.tau_e) * (np.cos(self.phiav[i]) * self.xi_A[i] + np.sin(self.phiav[i])*self.xi_B[i]) \
                        * np.sqrt(self.Dt) + np.sqrt(2*self.Nav[i]/self.tau_e)*self.xi_C[i]*np.sqrt(self.Dt)                        
                        noise_Q[i] = 2*np.sqrt( self.Csp*self.Nav[i] *self.Qav[i]/2/self.tau_e)*(self.xi_B[i]*np.cos(self.phiav[i]) + self.xi_B[i]*np.sin(self.phiav[i]) )* np.sqrt(self.Dt)                    
                        noise_phi[i] = np.sqrt(self.Csp * self.Nav[i]/(2*self.tau_e*self.Qav[i]))*(self.xi_B[i]*np.cos(self.phiav[i]) - self.xi_A[i]*np.sin(self.phiav[i]))*np.sqrt(self.Dt)
                
                #extract Q nad Phi from Master Field
                phiM = np.zeros(self.L)
                nuM = np.zeros(self.L)
                QM = np.zeros(self.L)
                add = np.zeros(self.L)
                for i in range(self.L):
                    QM[i] = abs(EM[i])**2
                for i in range(self.L-1):
                    phiM[i+1] = phiM[i] +  self.Dt/abs(EM[i])**2 *  ( (EM[i+1] - EM[i])/self.Dt * np.conjugate(EM[i])).imag
                for i in range(0,self.L-1):
                    nuM[i] = (phiM[i+1] - phiM[i]) /self.Dt/2/np.pi
                nuM[self.L-1] = nuM[self.L-2]    
                    
                if (temp == True):
                    for i in range(0, self.L-1):
                        if (biexp == False):
                            self.DT[i+1] = self.DT[i] - self.DT[i]/ self.th * self.Dt + 0.8 * self.rh/self.th * (1 -self.eta) * (I[i] - self.Ith) * self.Dt + self.RL*self.rh/self.th*(I[i]**2 - self.Ith**2)*self.Dt
                        else:
                            self.DT1[i+1] = self.DT1[i] +( -(self.DT1[i] - self.rh*1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.r1 - self.k12*(self.DT1[i] -self.DT2[i]) )*self.Dt
                            self.DT2[i+1] = self.DT2[i] +( -(self.DT2[i] - self.rh*1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.r2 - self.k12*(self.DT2[i] -self.DT1[i]) )*self.Dt
                            self.DTJ[i+1] = self.DTJ[i] +( -self.DTJ[i]/self.tauj + self.rh*self.Cj/self.tauj * 1.24/self.lambd*(1 -self.eta)*(I[i]-self.Ib) )*self.Dt
                    if (biexp == True):
                      for i in range(0, self.L):
                          self.DT[i] = (self.DT1[i] + self.DT2[i])/2 + self.DTJ[i]          
                      
                #datawrite = []
                for i in range(0, self.L-1):        

                    if (field == False):
                        
                        if( injfix ):
                            c = 2
                        else:
                            c = 1
                            
                        self.N[i+1] = self.N[i] + I[i]/e*self.Dt - self.N[i]/self.v3 - \
                        self.Q[i]/self.c1*(self.N[i] - self.Ntr)/(self.c2) \
                        / np.sqrt((1+2*self.gammaQ*1e9*self.Q[i]))*self.Dt + noise_N[i]

                        self.Q[i+1] = self.Q[i] + ((self.N[i]- self.Ntr)/(self.c2) \
                        / np.sqrt((1 + 2*self.gammaQ*1e9*self.Q[i])) - 1)*self.Q[i]/self.v4 + \
                        self.Csp*self.N[i]/self.tau_e*self.Dt  + c*self.kap_inj*np.sqrt(self.Q[i] * QM[i])*np.cos(self.w_inj*driver.time[i] + phiM[i] - self.phi[i] ) *self.Dt + noise_Q[i]  
                        
                        if(self.N[i+1] < 0): self.N[i+1] = 1e-12
                        if(self.Q[i+1] < 0): self.Q[i+1] = self.Qav[i+1]
                        
                        if( temp == False):
                            self.phi[i+1] = self.phi[i] +self.alpha/2/self.tau_ph*((self.N[i] - self.Ntr)/(self.c2) - 1)*self.Dt   +\
                            self.kap_inj* np.sqrt(QM[i]/self.Q[i])*np.sin(self.w_inj*driver.time[i+1] + phiM[i] - self.phi[i])*self.Dt + noise_phi[i]  
                            
                            add[i] = np.sqrt(QM[i]/self.Q[i])*self.kap_inj*np.sin(self.w_inj*driver.time[i+1]+ phiM[i] - self.phi[i])*self.Dt
                        elif ( temp == True):
                            self.phi[i+1] = self.phi[i] +(self.c4*((self.N[i] - self.Ntr)/(self.c2) - 1) -  self.kap_om *self.DT[i])*self.Dt +noise_phi[i] +\
                            self.kap_inj* np.sqrt(QM[i]/self.Q[i])*np.sin(self.w_inj*driver.time[i] + phiM[i] - self.phi[i])*self.Dt  
                            ''' self.kap_inj* np.sqrt(laserM.Q[i]/self.Q[i])\
                            *np.sin(self.w_inj*driver.time[i] + laserM.phi[i] - self.phi[i])*self.Dt + '''  
                    #TODO Field calculatio
                    # if (field):
                    #     #if(i != 0):
                    #     self.Q[i] = abs(self.E[i])**2
                    #     #Q_E = abs(self.E[i])**2
                    #     self.N[i+1] = self.N[i] + I[i]/e*self.Dt - self.N[i]/self.v3 - \
                    #     self.Q[i]/self.c1*(self.N[i] - self.Ntr)/(self.c2)* \
                    #     self.Dt/np.sqrt((1 + 2*self.v2*abs(self.E[i])**2)) -(np.sqrt(2*self.Csp*self.N[i]*abs(self.E[i])**2/self.tau_e)*self.xi_1[i]-np.sqrt(2*self.N[i]/self.tau_e)*self.xi_3[i])*np.sqrt(self.Dt)
                    #     if (temp == True):
                    #         self.E[i+1] = self.E[i] + self.E[i] * self.Dt / self.tau_ph/2 *( ((self.N[i] - self.Ntr)/self.c2)/np.sqrt((1 + 2*self.v2*abs(self.E[i])**2 )) - 1 +  1j*self.alpha*((self.N[i] - self.Ntr)/self.c2-1) \
                    #         -2*1j*self.tau_ph*self.kap_om*self.DT[i]) +\
                    #         np.sqrt( self.Csp * self.N[i] /2/self.tau_e )* (self.xi_1[i] *np.sqrt(self.Dt) + self.xi_2[i]*np.sqrt(self.Dt)*1j )
                    #     else:
                    #         self.E[i+1] = self.E[i] + self.E[i]*self.Dt/self.tau_ph/2 * ( (self.N[i] - self.Ntr)/self.c2/np.sqrt((1 + 2*self.v2*abs(self.E[i])**2)) - 1 +  1j*self.alpha*( (self.N[i] - self.Ntr)/self.c2-1) ) \
                    #         + np.sqrt( self.Csp * self.N[i] /2/self.tau_e )* (self.xi_1[i]**2 *np.sqrt(self.Dt) + self.xi_2[i]**2*np.sqrt(self.Dt)*1j )
                    #     if( phase == True):
                    #         self.phi[i+1] = self.phi[i] + self.Dt/abs(self.E[i])**2 *  ( (self.E[i+1] - self.E[i])/self.Dt * np.conjugate(self.E[i]) ).imag
                
                #datawrite.append()
                #text = np.c_[np.stack(datawrite, axis = 1)]  
                #np.savetxt('add.txt', datawrite)         
     
                for i in range(0,self.L-1):
                    self.nu[i] = (self.phi[i+1] - self.phi[i]) /self.v5
                self.nu[self.L-1] = self.nu[self.L-2]

                if (interference == True):
                    #0.8 sec - len(8000)

                    shift = int(self.delay / self.Dt)
                    print( str(shift) + " slave interferometer delay")
                    Q1 = self.Q[0:self.L-shift]
                    Q2 = self.Q[shift:self.L]
                    phi1 = self.phi[0:self.L-shift]
                    phi2 = self.phi[shift:self.L]
                    for i in range(0,self.L-shift):
                        self.S[i] = abs( np.sqrt(Q1[i]) * np.exp(1j * phi1[i] + 1j*self.phase_m1) + np.sqrt(Q2[i]) * np.exp(1j * phi2[i] + 1j*self.phase_m2) )**2
                        
                self.E = np.sqrt(self.Q) * np.exp(1j*self.phi)
                Data =  [self.N,self.E,self.DT,self.phi,self.nu]                   
                return Data
            def phasepack(self,field):

                phi = np.zeros(len(field))
                for i in range(0, len(field) - 1):
                    phi[i+1] = phi[i] + self.Dt/abs(field[i])**2 *  ( (field[i+1] - field[i])/self.Dt * np.conjugate(field[i]) ).imag
                return phi
            def frequency(self,phi):
                for i in range(0,self.L-1):
                    self.nu[i] = (phi[i+1] - phi[i]) /self.v5
                self.nu[self.L-1] = self.nu[self.L-2]
                return self.nu  
            def interference(self, field):
                
                shift = int(self.delay / self.Dt)

                field1 = field[0:self.L-shift]
                field2 = field[shift:self.L]

                S = np.zeros(len(field))
                for i in range(0,self.L-shift):
                    S[i] = abs(field1[i]*np.exp(1j*(self.phase_m1-self.phase_m2)) + field2[i] )**2
                    
                    
                return S
            def interference_Qphi(self, Q, phi):
                S = np.zeros(len(Q))
                shift = int(self.delay / self.Dt)
                print( str(shift) + " slave interferometer delay")
                Q1 = Q[0:self.L-shift]
                Q2 = Q[shift:self.L]
                phi1 = phi[0:self.L-shift]
                phi2 = phi[shift:self.L]
                for i in range(0,self.L-shift):
                    S[i] = abs( np.sqrt(Q1[i]) * np.exp(1j * phi1[i] + 1j*self.phase_m1) + np.sqrt(Q2[i]) * np.exp(1j * phi2[i] + 1j*self.phase_m2) )**2
                                
                return S
            def transfer(self,x,Is):

                print( "Is  " + str(Is) +  " Ith " + str(self.Ith) +  " Itr " + str(self.Itr))
                print( "gq " + str(self.v2) + " tph " + str(self.tau_ph) + " Г " + str(self.Gamma) + " Ntr/Nth " + str(self.Ntr/self.Nth))
                q  = (self.Gamma * self.tau_ph)/e * (self.Ith - self.Itr)
                B = np.sqrt(1 + (self.v2**2)*(q**2)+ 2*self.Gamma*self.tau_ph/e*(Is-self.Itr)*self.v2)                
                A = ( B - 1 - q *self.v2)

                Qs = (self.Gamma *self.tau_ph)/e *(Is - self.Ith) - q*A
                Ns = self.Nth + self. tau_e/(self.Gamma * self.tau_ph)*q*A
                
                dGQ = -(Ns - self.Ntr)/(self.Nth - self.Ntr) *self.v2/((1+2*self.v2*Qs)**(3/2))
                            
                dGN = 1/(self.Nth-self.Ntr)*1/np.sqrt(1+2*self.v2*Qs)
                            
                Gs = (Ns - self.Ntr)/(self.Nth - self.Ntr) * 1/(np.sqrt(1+2*self.v2*Qs))
                             
                b1 = 1/e
                a11 = -1/self.tau_e - Qs/(self.tau_ph * self.Gamma)*dGN
                            
                a12= -1/(self.Gamma * self.tau_ph)*( Gs + Qs*dGQ)
                            
                a21 = Qs/self.tau_ph * dGN
                            
                a22 = 1/self.tau_ph * (Gs - 1 + Qs*dGQ)
                wr = np.sqrt(Qs/(self.Gamma*self.tau_ph**2)*dGN - 1/(self.tau_e*self.tau_ph)*(Gs-1+Qs*dGQ) )
                wd = (wr**2)/(-a11-a22)
                
                
                return  abs(1 /( (1j*x/wr)**2 + (1j*x/wd) + 1))**2      
            def stability_diagram(self,Qb,Qe,phib,phie,guess,step):
                Nx = step
                print(Nx)
                eta = self.eta
                tau_ph = self.tau_ph 
                Ntr = self.Ntr
                Nth = self.Nth
                chi = self.chi
                
                Csp = self.Csp
                Gamma = self.Gamma
                tau_e = self.tau_e
                alpha = self.alpha
                
                gammaQ = self.gammaQ
                i = 20.*1e-12
                k_inj = self.kap_inj
                if (int(Qb) < 1):
                    Qb = 1
                phi_l = np.linspace(phib, phie,int((abs(phib)+abs(phie)+0.005)/0.005))

                def fun(Qs,phis):
                    return (Csp*(e*Ntr*Qs+ i*(Nth-Ntr)*Gamma*np.sqrt(1+2*Qs*gammaQ)*tau_ph))/(e*(Qs*tau_e +(Nth-Ntr)*Gamma*np.sqrt(1+2*Qs*gammaQ)*tau_ph )) +\
                        Qs/tau_ph*(-1 + (-Ntr + (tau_e*(e*Ntr*Qs+i*(Nth-Ntr)*Gamma*np.sqrt(1+2*Qs*gammaQ)*tau_ph ))/(e*(Qs*tau_e+(Nth-Ntr)*Gamma*np.sqrt(1+2*Qs*gammaQ)*tau_ph)) ) \
                        /((Nth-Ntr)*np.sqrt(1+2*abs(Qs)*gammaQ)) ) + 2*np.sqrt(abs(Qs) * Qinj)*k_inj*np.cos(phis)
                #TODO !  is it ok to put QS in abs?
                maxp = int((2*3.14/0.02*10)*(10**4/50))
                points=[ np.empty(maxp),np.empty(maxp),np.empty(maxp),np.empty(maxp),np.empty(maxp)]
                pointsnum = 0
                Qinj = int(Qb)
               
                while (Qinj < int(Qe)):
                    inStablePrevious = 0
                    inStable = 0
                    for j in range(0,len(phi_l)):
                
                        Qs = optimize.root(fun, [guess],phi_l[j], method='hybr').x[0]        
                        Ns = (tau_e*(e*Ntr*Qs + i*(Nth-Ntr)*Gamma*np.sqrt(1+2*Qs*gammaQ)*tau_ph))/(e*(Qs*tau_e+(Nth-Ntr)*Gamma*np.sqrt(1+2*gammaQ*Qs)*tau_ph))        
                        dOmega = np.sqrt(Qinj/Qs)*k_inj*(np.sin(phi_l[j]) - alpha*np.sqrt(1+2*gammaQ*Qs)*np.cos(phi_l[j]))
                   
                        c0 = 1        
                        c1 = ( (Qs*tau_e)/((Nth - Ntr)*np.sqrt(1+ 2*Qs*gammaQ)) - Gamma*(-(((Ns - Ntr)*Qs*gammaQ)/ ((Nth - Ntr)*(1 + 2*Qs*gammaQ)**(3/2))) + \
                            (Ns - Ntr)/((Nth - Ntr)*np.sqrt(1+2*Qs*gammaQ)))*tau_e + Gamma*(tau_e + tau_ph))/(Gamma*tau_e*tau_ph)
                        c2 = ((Qs**2 * tau_e) / ((Nth - Ntr) * np.sqrt(1+2*Qs*gammaQ)) + Qs*Gamma*tau_ph + np.sqrt(Qinj*Qs)*Gamma*k_inj*tau_e*tau_ph*np.cos(phi_l[j]) + \
                            ((Ns - Ntr)*tau_ph*(Qs *(Csp - Gamma) - np.sqrt(Qinj*Qs)*Gamma*k_inj*tau_e*np.cos(phi_l[j]))) / ((Nth - Ntr)*np.sqrt(1+2*Qs*gammaQ)) - \
                            ((Ns - Ntr)* gammaQ * tau_ph * (Qs**2*(Csp - Gamma) - np.sqrt(Qinj*Qs**3)*Gamma*k_inj*tau_e*np.cos(phi_l[j]))) / ((Nth-Ntr)*(1+2*Qs*gammaQ)**(3/2)) -\
                            Qinj*Gamma*k_inj**2*tau_e*tau_ph**2*np.cos(2*phi_l[j]))/(Qs * Gamma* tau_e*tau_ph**2)
                        c3 = ( k_inj*((np.sqrt(Qinj*Qs**3)*tau_e*np.cos(phi_l[j])) / ((Nth - Ntr)*np.sqrt(1+2*Qs*gammaQ)) + np.sqrt(Qinj*Qs)*Gamma*tau_ph*np.cos(phi_l[j]) - \
                            (Qinj*Qs*k_inj*tau_e*tau_ph*np.cos(2*phi_l[j]))/((Nth - Ntr)*np.sqrt(1+2*Qs*gammaQ)) - Qinj*Gamma*k_inj*tau_ph**2*np.cos(2*phi_l[j]) + \
                            ((Ns - Ntr)*np.sqrt(Qinj*Qs) * ((Csp - Gamma)*tau_ph*np.cos(phi_l[j]) + (Qs* alpha* tau_e *np.sin(phi_l[j]))/ (Nth - Ntr))) / ((Nth - Ntr)*np.sqrt(1 + 2*Qs*gammaQ)) - \
                            ((Ns - Ntr)* gammaQ * (np.sqrt(Qinj*Qs**3) *(Csp - Gamma)*tau_ph*np.cos(phi_l[j]) + (np.sqrt(Qinj*Qs**5)*alpha*tau_e*np.sin(phi_l[j]))/ (Nth - Ntr))) / \
                            ((Nth - Ntr)*(1 + 2*Qs*gammaQ)**(3/2)))) / (Qs * Gamma * tau_e * tau_ph**2)
                        if (c1 > 0 and c2 > 0 and c3 > 0 and (c1*c2 - c0*c3) > 0):
                            inStable = 1
                        else:
                            inStable = 0                        
                        if (inStable != inStablePrevious):
                            points[0][pointsnum] = Qinj
                            points[1][pointsnum] = dOmega/2/np.pi
                            points[2][pointsnum] = phi_l[j]
                            points[3][pointsnum] = Qs
                            points[4][pointsnum] = Ns         
                            pointsnum = pointsnum +1
                            A = 10*np.log10(Qinj/Qs) + 1/Nx*(10*np.log10(int(Qe)/Qs) - 10*np.log10(int(Qb)/Qs) )
                            Qinj = Qs*10**(A/10)
                        inStablePrevious = inStable 
                arr = points[0][0:pointsnum]/points[3][0:pointsnum]
                print(pointsnum)
                return [10*np.log10(arr),points[1][0:pointsnum]]
                                                         