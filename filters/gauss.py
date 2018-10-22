from . import base_filter
import numpy as np
from scipy import signal
from scipy import special

class Gauss(base_filter):

    #Filter initialization with initial parameters received
    def __init__(self, *args):
        if (args[0]=='Group Delay'):
            self.name =args[0]
            self.gain=args[1]
            self.tao0=args[2]/1000
            self.wrg=args[3]
            self.palm=args[4]/100
            self.n=args[5]
            self.qmax=args[6]
            self.nmax=1000 #VALOR NO DEFINITIVO, PROBAR CUAL ES EL VALOR MAXIMO PARA EL QUE EMPIEZA A MORIR LA APROXIMACION
            self.error=self.check_input()
            if self.error is False:
                self.poles=[]
                self.zeroes=[]
                self.den=np.poly1d([1])
                self.num=np.poly1d([1])
                self.normalize()#Normalizes current template
                self.do_approximation()#Does normalized approximation and realizes
                self.denormalize()#Denormalizes the approximation to match desired template
    
    def do_approximation(self):
        gamma=1
        if self.n == None:
            self.n=0    #N is calculated through iterations, starting from 0
            Tn_Wrgn=0   #Group Delay at normalized frequency(wrgn)(starting value is fixed for algorithmic purposes)
            while  (((1-self.palm))<= Tn_Wrgn)==False: #Check if tolerance limits are met
                self.n += 1
                g_coef=np.poly1d(self.gauss_coef(self.n,gamma))#Gauss n-th order polynomial is calculated
                g_roots=np.roots(g_coef)# This polynomial roots are |H(s)|^2 poles
                den=np.poly1d([1])
                for i in range(0,len(g_roots)):
                    if np.real(g_roots[i])<=0:#Only left plane poles are considered
                        den=den*np.poly1d([-1/g_roots[i],1])#H(s) transfer function is created
                w=np.array([0,0.0000001,self.wrgn*0.9999999,self.wrgn,self.wrgn*1.0000001])#Vector containing wrgn and 0 and proximate values for calculating group delay numerically
                w, h = signal.freqs([1],den,w)#Frecuency response for these given frequencies is calculated
                dphase=np.ediff1d(np.angle(h,deg=False))
                dw=np.ediff1d(w)
                gd=-dphase/dw#Group delay is calculated
                Tn_0=gd[0]#For normalization purposes, group delay at wn=0 should be 1s
                den=np.poly1d([1])
                for i in range(0,len(g_roots)):#Using group delay at 0 calculated previously, filter is normalized to meet previous requirement(group delay at w=0 equal to 1s)
                    if np.real(g_roots[i])<=0:#This is realized as in contrast to bessel filters, group delay at w=0 is not assured to be 1s
                        den=den*np.poly1d([-1/(g_roots[i]*Tn_0),1])#Normalized transfer function is calculated
                w, h = signal.freqs([1],den,w)
                dphase=np.ediff1d(np.angle(h,deg=False))
                dw=np.ediff1d(w)
                gd=-dphase/dw#Normalized Group_delay is calculated
                Tn_0=gd[0]#Debug value, used to check if filter was correctly denormalized
                Tn_Wrgn=gd[3]#Group delay at w=wrgn is calculated
        else:
            g_coef=np.poly1d(self.gauss_coef(self.n,gamma))#Gauss n-th order polynomial is calculated
            g_roots=np.roots(g_coef)# This polynomial roots are |H(s)|^2 poles
            den=np.poly1d([1])
            for i in range(0,len(g_roots)):
                if np.real(g_roots[i])<=0:#Only left plane poles are considered
                    den=den*np.poly1d([-1/g_roots[i],1])#H(s) transfer function is created
            w=np.array([0,0.001])#Vector containing 0 and proximate values for calculating group delay numerically
            w, h = signal.freqs([1],den,w)#Frecuency response for these given frequencies is calculated
            dphase=np.ediff1d(np.angle(h,deg=False))
            dw=np.ediff1d(w)
            gd=-dphase/dw#Group delay is calculated
            Tn_0=gd[0]#For normalization purposes, group delay at wn=0 should be 1s
            den=np.poly1d([1])
            for i in range(0,len(g_roots)):#Using group delay at 0 calculated previously, filter is normalized to meet previous requirement(group delay at w=0 equal to 1s)
                if np.real(g_roots[i])<=0:#This is realized as in contrast to bessel filters, group delay at w=0 is not assured to be 1s
                    den=den*np.poly1d([-1/(g_roots[i]*Tn_0),1])#Normalized transfer function is calculated
        self.aprox_gain=1
        self.norm_sys = signal.TransferFunction([1],den) #Filter system is obtained
        self.poles=np.roots(den)

    def gauss_coef (self,n,gamma):#Function calculates k-th bessel coefficients for a n-order bessel polynomial. Receives a vector of decreasing integers starting from n
        g_k=[]
        for i in range(0,n+1):
            g_k.append(np.power(gamma,i)*np.power(1/1j,i*2)/special.factorial(i))
            g_k.append(0)
        #b_k is a vector with all of the polynomial coefficients, with decreasing order
        return np.flip(g_k)


    #def factorial(self,n):
    #    if n <= 0:
    #        return 1
    #    else:
    #        return n*self.factorial(n-1)

    #def do_approximation(self):
    #    #Me pasan Tao(0) WRG y palmerita
    #    wrgn = self.tao0*self.wrg
    #    j = np.complex(0,1)
    #    NoCumplo=True
    #    self.n=1
        
    #    while NoCumplo or self.n<=20:
    #        Terminos= np.zeros(self.n*2)
    #        Terminos[len(Terminos)]=1
           
    #        for i in range (0,self.n+1):
    #            Terminos[len(Terminos)-2*i] = (-1*(j**i))/self.factorial(i)

    #        POL=np.poly1d(Terminos)
    #        den=Terminos
    #        w,TaoWRGN = signal.group_delay((1,den),[wrgn])

    #        if TaoWRGN >= (1-self.palm):
    #            NoCumplo=False
    #        else:
    #            NoCumplo=True

    #        self.poles=np.roots(POL)
    #        self.zeroes=1
    #        self.n=self.n+1
                