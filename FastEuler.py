from numba import jit,types
import numba
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy import optimize

e = 1.6e-19        #The electron charge, C
h = 6.63e-34       #The Planck constant, J*s
nu0 = 1.93548e14
@numba.jit(nopython=True)

def fast_euler(Ib,Ip,I,Dt,L,nPointsPerPulse,arg = np.array(38,dtype = np.float64),phase = None,freq = None, noise = None, temp = None, interference = None,field = None,biexp = None):
                print(arg)
                Ib = Ib
                Ip = Ip
               
                eta = arg[0] #0.3    

                tau_ph =arg[1] #0.001              #The photon lifetime, ns 
                Ntr = arg[2]*1e7 #4e7#6e7          #The transparency carrier number
                Nth = arg[3]*1e7#5.5e7#6.5e7       #The threshold carrier number
                chi = arg[4]#20                    #The gain compression factor, 1/W; 
                                               # chi*eta*h*nu0/2/Gamma/tau_ph = chi_Q is the dimensionless compression factor
                Csp = arg[5]*1e-5                  #The fraction of spontaneously emitted photons that end up in the active mode
                Gamma = arg[6]                     #The confinement factor
                tau_e = arg[7]                     #The effective carrier lifetime, ns
                alpha = arg[8]#5.0#6.0             #The linewidth enhancement factor                
                
                Ith = e *  Nth/ tau_e      #Threshold current
                Itr = e *  Ntr/ tau_e      #Transparency current
                print("Threshold current :" ,  Ith)
                #nPointsPerPulse = Driver.nPointsPerPulse
                print("N poits per pulse :" ,  nPointsPerPulse)
                print("Array size" ,  L)
                
                #Dt = Driver.Dt
               

                #I = np.zeros(L,dtype = types.float64)
                N = np.zeros( L,dtype = types.float64)               #carrier number array
                Q = np.zeros( L,dtype = types.float64)               #laser intensity
                phi = np.zeros( L,dtype = types.float64)             #phase 
                nu = np.zeros( L,dtype = types.float64)              #frequency
                DT = np.zeros( L,dtype = types.float64)              #temperature difference
                E = np.zeros(L, dtype=types.complex128)
                phiS = np.zeros( L,dtype = types.float64)
                DT1 = np.zeros( L,dtype = types.float64)
                DT2 = np.zeros( L,dtype = types.float64)
                DTJ = np.zeros( L,dtype = types.float64)
                DTS = np.zeros( L,dtype = types.float64)
                
                
                
    
                Nav = np.zeros( L,dtype = types.float64)             #array for average solutions without noise
                Qav = np.zeros( L,dtype = types.float64)
                phiav = np.zeros( L,dtype = types.float64)
                
                S = np.zeros( L,dtype = types.float64)               #interference
                teta1 = 0                          #params for phase in interference
                teta2 = 0
                
                N[0] =  Nth                    #initial values
                Q[0] = 1
                DT[0] = 0
                S[0] = 0
                
                #temperature effects
                rh = arg[9]*1e9#10 * 1e9           #thermal resistance , K/W
                th = arg[10]#10                    #thermal rise time , ns
                lambd = arg[11]#1.55               #wavelength , micrometers
                kap_om = arg[12]#2 * np.pi * 10    #temperature coefficient of frequency , GHZ                
                #injection parameters
                t1 = arg[13]#0.7                   # Field transmission coefficient of the facet through which the external optical field is injected into the laser cavity
                w_M = arg[14]#0
                om_th = arg[15]#0
                Ll = arg[16]#0.001                 # length of the resonator , micrometers
                kap_inj =  arg[17]#33.095          # coupling coefficient.
                r1 = arg[18]
                r2 = arg[19]

                rhs = arg[20] * 1e9 
                k12 = arg[21]
                tauj = arg[22]
                lambdas = arg[23]
                Cj = arg[24]
                RL = 8 * 1e9 #resistance
                w_inj = 2*np.pi*arg[14]            #2*np.pi*( w_M -  om_th)    use it instead of om_M and Om_th
                 #interferometer delay
                delay = arg[15]
                #constant values for optimization 
                v1 = 2/ tau_e
                v2 = 1e9* chi* eta*h*nu0/2/ Gamma/ tau_ph #chi*eta*h*nu0/2/Gamma/tau_ph = chi_Q is the dimensionless compression factor
                gammaQ =  chi* eta*h*nu0/2/ Gamma/ tau_ph
                
                v3 =  tau_e*1/ Dt 
                v4 =  tau_ph*(1/ Dt)
                v5 = ( Dt*2*np.pi)
                c1 =  Gamma/(1/ tau_ph)
                c2 =  Nth -  Ntr
                c3 = 1e9* chi *  eta * h * nu0/ Gamma/ tau_ph
                c4 =  alpha/2/ tau_ph 

               

                Nav[0] =  Nth
                Qav[0] = 1
                #I = np.zeros(L)
                N[0] =  Nth            #initial values
                Q[0] = 1
                DT[0] = 0
                DT1[0] = 0
                DT2[0] = 0
                DTJ[0] = 0
                DTS[0] = 0
                
                S[ 0] = 0
                if (temp == True):
                    for i in range(0,  L-1):
                        if (biexp == False):
                             DT[i+1] =  DT[i] -  DT[i]/  th *  Dt + 0.8 *  rh/ th * (1 - eta) * (I[i] -  Ith) *  Dt +  RL* rh/ th*(I[i]**2 -  Ith**2)* Dt
                        else:
                             DT1[i+1] =  DT1[i] +( -( DT1[i] -  rhs*1.24/ lambd*(1 - eta)*(I[i]- Ib) )* r1 -  k12*( DT1[i] - DT2[i]) )* Dt
                             DT2[i+1] =  DT2[i] +( -( DT2[i] -  rhs*1.24/ lambd*(1 - eta)*(I[i]- Ib) )* r2 -  k12*( DT2[i] - DT1[i]) )* Dt
                             DTJ[i+1] =  DTJ[i] +( - DTJ[i]/ tauj +  rhs* Cj/ tauj * 1.24/ lambd*(1 - eta)*(I[i]- Ib) )* Dt
                    if (biexp == True):
                      for i in range(0,  L):
                           DT[i] = ( DT1[i] +  DT2[i])/2 +  DTJ[i]  
                          
                if ( field == False):
                     N[0] =  Nav[0]
                     Q[0] =  Qav[0]
                     phi[0] =  phiav[0]
                     for i in range(0,  L - 1):
                         Nav[i+1] =  Nav[i] + I[i]/e* Dt -  Nav[i]/ v3 -  Qav[i]/ Gamma/ tau_ph*( Nav[i]- Ntr) \
                        /( c2) / np.sqrt(1 +  c3* Qav[i])* Dt
    
                         Qav[i+1] =  Qav[i] + (( Nav[i] -  Ntr)/( c2)/ np.sqrt(1+  c3 \
                        *  Qav[i]) - 1) *  Qav[i] /  tau_ph *  Dt +  Csp* Nav[i] /  v3
    
                         phiav[i+1] =  phiav[i] +  c4 * (( Nav[i] -  Ntr)/( c2)-1) *  Dt
            
          
                counter = 0
               
                xi_N = np.random.randn( L)
                xi_Q = np.random.randn( L)
                xi_phi = np.random.randn( L)
                xi_1 = np.random.randn( L)
                xi_2 = np.random.randn( L)
                xi_3 = np.random.randn( L)
                if(noise):
                    print("NOISE ON")
                else:
                    print("NOISE OFF")
                for i in range(0,  L-1):
                    if(noise):
                        if ( Nav[i] <= 0):
                                 Nav[i] = 1e-12
                        if ( Qav[i] <= 0):
                                 Qav[i] = 1e-12
                                 
                                 
                        
                        noise_N = -2*np.sqrt( Csp* Nav[i]* Qav[i]/2/ tau_e) * (np.cos( phiav[i]) *  xi_Q[i] + np.sin( phiav[i])* xi_phi[i]) \
                        * np.sqrt( Dt) + np.sqrt(2* Nav[i]/ tau_e)* xi_N[i]*np.sqrt( Dt)
                        
                        noise_Q = 2*np.sqrt(  Csp* Nav[i] * Qav[i]/2/ tau_e)*( xi_Q[i]*np.cos( phiav[i]) +  xi_phi[i]*np.sin( phiav[i]) )*np.sqrt( Dt) 
                        
                        noise_phi = np.sqrt( Csp *  Nav[i]/2/ tau_e/ Qav[i])*( xi_phi[i]*np.cos( phiav[i]) -  xi_Q[i]*np.sin( phiav[i]))*np.sqrt( Dt)
                    else:
                        noise_N = 0
                        noise_Q = 0
                        noise_phi = 0                        
                    
                    if (1 > 0):
                         counter = counter + 1
                         N[i+1] = N[i] + I[i]/e* Dt -  N[i]/ v3 - \
                         Q[i]/ c1*( N[i] -  Ntr)/( c2) \
                         / np.sqrt((1+2* gammaQ*1e9* Q[i]))* Dt + noise_N
                        
                         Q[i+1] =  Q[i] + (( N[i]-  Ntr)/( c2) \
                         / np.sqrt((1 + 2* gammaQ*1e9* Q[i])) - 1)* Q[i]/ v4 + \
                         Csp* N[i]/ v3 + noise_Q  
                        
                         if( N[i+1] < 0):  N[i+1] = 1e-3
                         if( Q[i+1] < 0):  Q[i+1] = Qav[i+1]
                        
                               
                        
                        
                        
                         if(  temp == False):
                             phi[i+1] =  phi[i] +  c4*(( N[i] -  Ntr)/( c2) - 1)* Dt + noise_phi
                         elif ( temp == True):
                            phi[i+1] =  phi[i] + ( c4*(( N[i] -  Ntr)/( c2) - 1) -   kap_om * DT[i])* Dt + noise_phi
                    
                if(freq):
                    for i in range(0, L-1):
                         nu[i] = ( phi[i+1] -  phi[i]) / Dt/2/np.pi
                    nu[ L-1] =  nu[ L-2]
                #EO = np.multiply(np.sqrt(Q),np.exp( np.multiply(phi, 1j )))
               # E = np.exp(1j*phi)
                #print(counter)
                Data =  [ N,  Q, DT, phi, nu]
                       
                return Data

@numba.jit(nopython=True) 
def euler_m_slave(Ib,Ip,I,Dt,L,nPointsPerPulse,time,EM,injfix,arg = np.array(38,dtype = np.float64),noise = None,temp = None,interference = None, field = None ,injection = None,phase = None,freq = None,biexp =None):   

                print(arg)
                Ib = Ib
                Ip = Ip
               
                eta = arg[0] #0.3    

                tau_ph =arg[1] #0.001              #The photon lifetime, ns 
                Ntr = arg[2]*1e7 #4e7#6e7          #The transparency carrier number
                Nth = arg[3]*1e7#5.5e7#6.5e7       #The threshold carrier number
                chi = arg[4]#20                    #The gain compression factor, 1/W; 
                                               # chi*eta*h*nu0/2/Gamma/tau_ph = chi_Q is the dimensionless compression factor
                Csp = arg[5]*1e-5                  #The fraction of spontaneously emitted photons that end up in the active mode
                Gamma = arg[6]                     #The confinement factor
                tau_e = arg[7]                     #The effective carrier lifetime, ns
                alpha = arg[8]#5.0#6.0             #The linewidth enhancement factor                
                
                Ith = e *  Nth/ tau_e      #Threshold current
                Itr = e *  Ntr/ tau_e      #Transparency current
                print("Threshold current :" ,  Ith)
                #nPointsPerPulse = Driver.nPointsPerPulse
                print("N poits per pulse :" ,  nPointsPerPulse)
                print("Array size" ,  L)
                
                #Dt = Driver.Dt
               
                #I = np.zeros(L,dtype = types.float64)
                N = np.zeros( L,dtype = types.float64)               #carrier number array
                Q = np.zeros( L,dtype = types.float64)               #laser intensity
                phi = np.zeros( L,dtype = types.float64)             #phase 
                nu = np.zeros( L,dtype = types.float64)              #frequency
                DT = np.zeros( L,dtype = types.float64)              #temperature difference
                E = np.zeros(L, dtype=types.complex128)
                phiS = np.zeros( L,dtype = types.float64)
                DT1 = np.zeros( L,dtype = types.float64)
                DT2 = np.zeros( L,dtype = types.float64)
                DTJ = np.zeros( L,dtype = types.float64)
                DTS = np.zeros( L,dtype = types.float64)            
            
            
                Nav = np.zeros( L,dtype = types.float64)             #array for average solutions without noise
                Qav = np.zeros( L,dtype = types.float64)
                phiav = np.zeros( L,dtype = types.float64)
                
                S = np.zeros( L,dtype = types.float64)               #interference
                teta1 = 0                          #params for phase in interference
                teta2 = 0
                
                N[0] =  Nth                    #initial values
                Q[0] = 1
                DT[0] = 0
                S[0] = 0
                
                #temperature effects
                rh = arg[9]*1e9#10 * 1e9           #thermal resistance , K/W
                th = arg[10]#10                    #thermal rise time , ns
                lambd = arg[11]#1.55               #wavelength , micrometers
                kap_om = arg[12]#2 * np.pi * 10    #temperature coefficient of frequency , GHZ                
                #injection parameters
                t1 = arg[13]#0.7                   # Field transmission coefficient of the facet through which the external optical field is injected into the laser cavity
                w_M = arg[14]#0
                om_th = arg[15]#0
                Ll = arg[16]#0.001                 # length of the resonator , micrometers
                kap_inj =  arg[17]#33.095          # coupling coefficient.
                r1 = arg[18]
                r2 = arg[19]

                rhs = arg[20] * 1e9 
                k12 = arg[21]
                tauj = arg[22]
                lambdas = arg[23]
                Cj = arg[24]
                RL = 8 * 1e9 #resistance
                w_inj = arg[14]#2*np.pi*arg[14]            #2*np.pi*( w_M -  om_th)    use it instead of om_M and Om_th
                 #interferometer delay
                delay = arg[15]
                #constant values for optimization 
                v1 = 2/ tau_e
                v2 = 1e9* chi* eta*h*nu0/2/ Gamma/ tau_ph #chi*eta*h*nu0/2/Gamma/tau_ph = chi_Q is the dimensionless compression factor
                gammaQ =  chi* eta*h*nu0/2/ Gamma/ tau_ph
                
                v3 =  tau_e*1/ Dt 
                v4 =  tau_ph*(1/ Dt)
                v5 = ( Dt*2*np.pi)
                c1 =  Gamma/(1/ tau_ph)
                c2 =  Nth -  Ntr
                c3 = 1e9* chi *  eta * h * nu0/ Gamma/ tau_ph
                c4 =  alpha/2/ tau_ph 

               

                Nav[0] =  Nth
                Qav[0] = 1
                #I = np.zeros(L)
                N[0] =  Nth            #initial values
                Q[0] = 1
                DT[0] = 0
                DT1[0] = 0
                DT2[0] = 0
                DTJ[0] = 0
                DTS[0] = 0
                
                S[ 0] = 0
                
                
                phiM = np.zeros( L,dtype = types.float64)
                nuM = np.zeros( L,dtype = types.float64)
                QM = np.zeros( L,dtype = types.float64)
                add = np.zeros( L,dtype = types.float64)
                
                
                for i in range(L):
                    QM[i] = abs(EM[i])**2
                for i in range(L-1):
                    phiM[i+1] = phiM[i] +  Dt/abs(EM[i])**2 *  ( (EM[i+1] - EM[i])/Dt * np.conjugate(EM[i])).imag
                for i in range(0,L-1):
                    nuM[i] = (phiM[i+1] - phiM[i]) /Dt/2/np.pi
                nuM[L-1] = nuM[L-2] 
            
            
                if (temp == True):
                    for i in range(0,  L-1):
                        if (biexp == False):
                             DT[i+1] =  DT[i] -  DT[i]/  th *  Dt + 0.8 *  rh/ th * (1 - eta) * (I[i] -  Ith) *  Dt +  RL* rh/ th*(I[i]**2 -  Ith**2)* Dt
                        else:
                             DT1[i+1] =  DT1[i] +( -( DT1[i] -  rhs*1.24/ lambd*(1 - eta)*(I[i]- Ib) )* r1 -  k12*( DT1[i] - DT2[i]) )* Dt
                             DT2[i+1] =  DT2[i] +( -( DT2[i] -  rhs*1.24/ lambd*(1 - eta)*(I[i]- Ib) )* r2 -  k12*( DT2[i] - DT1[i]) )* Dt
                             DTJ[i+1] =  DTJ[i] +( - DTJ[i]/ tauj +  rhs* Cj/ tauj * 1.24/ lambd*(1 - eta)*(I[i]- Ib) )* Dt
                    if (biexp == True):
                      for i in range(0,  L):
                           DT[i] = ( DT1[i] +  DT2[i])/2 +  DTJ[i]    
                           
                for i in range(0, L-1):        
                    # if(noise):
                    #     if (Nav[i] <= 0):
                    #             Nav[i] = 1e-12
                    #     if (Qav[i] <= 0):
                    #             Qav[i] = 1e-12
                    #     noise_N = -2*np.sqrt(Csp*Nav[i]*Qav[i]/2/tau_e) * (np.cos(phiav[i]) * xi_A[i] + np.sin(phiav[i])*xi_B[i]) \
                    #     * np.sqrt(Dt) + np.sqrt(2*Nav[i]/tau_e)*xi_C[i]*np.sqrt(Dt)                        
                    #     noise_Q = 2*np.sqrt( Csp*Nav[i] *Qav[i]/2/tau_e)*(xi_B[i]*np.cos(phiav[i]) + xi_B[i]*np.sin(phiav[i]) )*np.sqrt( Dt)                         
                    #     noise_phi = np.sqrt(Csp * Nav[i]/(2*tau_e*Qav[i]))*(xi_B[i]*np.cos(phiav[i]) - xi_A[i]*np.sin(phiav[i]))*np.sqrt(Dt)                        
                                  

                    # if (field == False):
                        
                    #     if( injfix ):
                    #         c = 1
                    #     else:
                    #         c = 2
                            
                    #     N[i+1] = N[i] + I[i]/e*Dt - N[i]/v3 - \
                    #     Q[i]/Gamma/tau_ph*(N[i] - Ntr)/(c2) \
                    #     / np.sqrt((1+2*gammaQ*1e9*Q[i]))*Dt # + noise_N[i]

                       
                    #     Q[i+1] = Q[i] + ((N[i]- Ntr)/(c2) \
                    #     / np.sqrt((1 + 2*gammaQ*1e9*Q[i])) - 1)*Q[i]/tau_ph*Dt + \
                    #     Csp*N[i]/tau_e*Dt  + c*kap_inj*np.sqrt(Q[i] * QM[i])*np.cos(w_inj*time[i] + phiM[i] - phi[i] ) *Dt #+ noise_Q[i]  
                        
                    #     if(N[i+1] < 0): N[i+1] = 1e-12
                    #     if(Q[i+1] < 0): Q[i+1] = Qav[i+1]
                        
                    #     #phi[i+1] = phi[i] + (c4*((N[i] - Ntr)/(c2) - 1) -  kap_om *DT[i])*Dt + noise_phi                   
                    #     if( temp == False):
                    #         phi[i+1] = phi[i] +alpha/2/tau_ph*((N[i] - Ntr)/(c2) - 1)*Dt   +\
                    #         kap_inj* np.sqrt(QM[i]/Q[i])*np.sin(w_inj*time[i+1] + phiM[i] - phi[i])*Dt #+ noise_phi[i]  
                            
                    #         add[i] = np.sqrt(QM[i]/Q[i])*kap_inj*np.sin(w_inj*time[i+1]+ phiM[i] - phi[i])*Dt
                    #     elif ( temp == True):
                    #         phi[i+1] = phi[i] +(c4*((N[i] - Ntr)/(c2) - 1) -  kap_om *DT[i])*Dt  +\
                    #         kap_inj* np.sqrt(QM[i]/Q[i])*np.sin(w_inj*time[i] + phiM[i] - phi[i])*Dt     # +noise_phi[i]      
              #E = np.sqrt(Q) * np.exp(1j*phi)
              
                  if (field == False):
                        
                        if( injfix ):
                            c = 1
                        else:
                            c = 2
                            
                        N[i+1] = N[i] + I[i]/e*Dt - N[i]/v3 - \
                        Q[i]/c1*(N[i] - Ntr)/(c2) \
                        / np.sqrt((1+2*gammaQ*1e9*Q[i]))*Dt #+ noise_N[i]

                        Q[i+1] = Q[i] + ((N[i]- Ntr)/(c2) \
                        / np.sqrt((1 + 2*gammaQ*1e9*Q[i])) - 1)*Q[i]/v4 + \
                        Csp*N[i]/tau_e*Dt  + c*kap_inj*np.sqrt(Q[i] * QM[i])*np.cos(w_inj*time[i] + phiM[i] - phi[i] ) *Dt #+ noise_Q[i]  
                        
                        if(N[i+1] < 0): N[i+1] = 1e-12
                        if(Q[i+1] < 0): Q[i+1] = Qav[i+1]
                        
                        #phi[i+1] = phi[i] + (c4*((N[i] - Ntr)/(c2) - 1) -  kap_om *DT[i])*Dt + noise_phi                   
                        if( temp == False):
                            phi[i+1] = phi[i] +alpha/2/tau_ph*((N[i] - Ntr)/(c2) - 1)*Dt   +\
                            kap_inj* np.sqrt(QM[i]/Q[i])*np.sin(w_inj*time[i+1] + phiM[i] - phi[i])*Dt #+ noise_phi[i]  
                            
                            add[i] = np.sqrt(QM[i]/Q[i])*kap_inj*np.sin(w_inj*time[i+1]+ phiM[i] - phi[i])*Dt
                        elif ( temp == True):
                            phi[i+1] = phi[i] +(c4*((N[i] - Ntr)/(c2) - 1) -  kap_om *DT[i])*Dt  +\
                            kap_inj* np.sqrt(QM[i]/Q[i])*np.sin(w_inj*time[i] + phiM[i] - phi[i])*Dt              
              
               
                for i in range(0,L-1):
                    nu[i] = (phi[i+1] - phi[i])/v5
                nu[L-1] = nu[L-2]
                Data =  [N,Q,DT,phi,nu]                   
                return Data            
                           

@numba.jit(nopython = True)           
def fun(Qs,phis,Qinj,args):
    eta =  args[0]
    tau_ph =  args[1] 
    Ntr =  args[2]
    Nth =  args[3]
    chi =  args[4]
                
    Csp =  args[5]
    Gamma =  args[6]
    tau_e =  args[7]
    alpha =  args[8]
                
    gammaQ =  args[9]
    i = 20.*1e-12
    k_inj =  args[10]
              
    return (Csp*(e*Ntr*Qs+ i*(Nth-Ntr)*Gamma*np.sqrt(1+2*Qs*gammaQ)*tau_ph))/(e*(Qs*tau_e +(Nth-Ntr)*Gamma*np.sqrt(1+2*Qs*gammaQ)*tau_ph )) +\
        Qs/tau_ph*(-1 + (-Ntr + (tau_e*(e*Ntr*Qs+i*(Nth-Ntr)*Gamma*np.sqrt(1+2*Qs*gammaQ)*tau_ph ))/(e*(Qs*tau_e+(Nth-Ntr)*Gamma*np.sqrt(1+2*Qs*gammaQ)*tau_ph)) ) \
        /((Nth-Ntr)*np.sqrt(1+2*Qs*gammaQ)) ) + 2*np.sqrt(Qs * Qinj)*k_inj*np.cos(phis)            
            
