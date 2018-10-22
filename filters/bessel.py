from . import base_filter
import numpy as np
from scipy import signal
from scipy import special

class Bessel(base_filter):

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
        if self.n == None:
            self.n=0    #N is calculated through iterations, starting from 0
            Tn_Wrgn=0   #Group Delay at normalized frequency(wrgn)(starting value is fixed for algorithmic purposes)
            while  ((1-self.palm)<= Tn_Wrgn)==False: #Check if tolerance limits are met
                self.n += 1
                k = np.linspace(0, self.n, num=(self.n + 1))  #vector of increasing integers from 0 to n is created
                k=np.flip(k)
                bess_coef=np.poly1d(self.bessel_coef(k))# K-th Bessel of N-th order is calculated
                w, h = signal.freqs(np.poly1d(self.bessel_coef(0)),bess_coef,self.wrgn)#Transfer function is calculated for wrgn
                Tn_Wrgn = 1-np.power(np.abs(h),2)*np.power(self.wrgn,2*self.n)/np.power(self.bessel_coef(0),2)#Group delay is calculated for wrgn
        else:
            k = np.linspace(0, self.n, num=(self.n + 1))  #vector of increasing integers from 0 to n is created
            k=np.flip(k)
            bess_coef=np.poly1d(self.bessel_coef(k))# K-th Bessel of N-th order is calculated
        self.aprox_gain=self.bessel_coef(0)#Transfer function has a gain constant
        self.norm_sys = signal.TransferFunction(self.aprox_gain,bess_coef) #Filter system is obtained
        self.poles=np.roots(bess_coef)#Poles are obtained

    def bessel_coef (self,k):#Function calculates k-th bessel coefficients for a n-order bessel polynomial. Receives a vector of decreasing integers starting from n
        b_k = (special.factorial(2*self.n-k))/((special.factorial(k))*(special.factorial(self.n-k))*(np.power(2,(self.n-k))))
        #b_k is a vector with all of the polynomial coefficients, with decreasing order from left to right
        return b_k