from . import base_filter
import numpy as np
from scipy import signal
from scipy import special

class Bessel(base_filter):

    #Filter initialization with initial parameters received
    def __init__(self, name,Ap,Ao,wpl,wph,wal,wah,gain,n,tao0=None,wrg=None,palm=None):
        if name:
            self.name = name
            self.Ap=Ap
            self.Ao=Ao
            self.wpl=wpl
            self.wph=wph
            self.wal=wal
            self.wah=wah
            self.n=n
            self.gain=gain
            self.tao0=tao0
            self.wrg=wrg
            self.palm=palm
            self.wan=1
            self.poles=[]
            self.zeroes=[]
            self.den=np.poly1d([1])
            self.num=np.poly1d([1])
            self.normalize()#Normalizes current template
            self.do_approximation()#Does normalized approximation and realizes
            self.denormalize()#Denormalizes the approximation to match desired template

    def do_approximation(self):
        self.n=0    #N is calculated through iterations, starting from 0
        Tn_Wrgn=0   #Group Delay at normalized frequency(wrgn)(starting value is fixed for algorithmic purposes)
        while  ((1-self.palm)<= Tn_Wrgn)==False: #Check if tolerance limits are met
            self.n += 1
            k = np.linspace(0, self.n, num=(self.n + 1))  #vector of increasing integers from 0 to n is created
            k=np.flip(k)
            bess_coef=np.poly1d(self.bessel_coef(k))# K-th Bessel of N-th order is calculated
            w, h = signal.freqs(np.poly1d(self.bessel_coef(0)),bess_coef,self.wrgn)#Transfer function is calculated for wrgn
            Tn_Wrgn = 1-np.power(np.abs(h),2)*np.power(self.wrgn,2*self.n)/np.power(self.bessel_coef(0),2)#Group delay is calculated for wrgn
        self.aprox_gain=self.bessel_coef(0)#Transfer function has a gain constant
        self.norm_sys = signal.TransferFunction([1],bess_coef) #Filter system is obtained
        self.poles=np.roots(bess_coef)#Poles are obtained

    def bessel_coef (self,k):#Function calculates k-th bessel coefficients for a n-order bessel polynomial. Receives a vector of decreasing integers starting from n
        b_k = (special.factorial(2*self.n-k))/((special.factorial(k))*(special.factorial(self.n-k))*(np.power(2,(self.n-k))))
        #b_k is a vector with all of the polynomial coefficients, with decreasing order from left to right
        return b_k