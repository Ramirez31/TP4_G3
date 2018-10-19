from . import base_filter
import numpy as np
from scipy import signal
from scipy import special
import mpmath as mp

class Cauer(base_filter):

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
            self.poles=[]
            self.zeroes=[]
            self.den=np.poly1d([1])
            self.num=np.poly1d([1])
            self.normalize()#Normalizes current template
            self.do_approximation()#Does normalized approximation and realizes
            self.denormalize()#Denormalizes the approximation to match desired template

    def do_approximation(self):
       epsilonp=np.sqrt(np.power(10,self.Ap/10)-1)
       gp=np.power(10,-self.Ap/20)
       k1=np.sqrt((np.power(10,self.Ap/10)-1)/(np.power(10,self.Ao/10)-1))
       k=1/self.wan
       a=mp.asin(1j/epsilonp)
       vo=mp.ellipf(1j/epsilonp,k1)/(1j*self.n)
       self.n=np.ceil(special.ellipk(np.sqrt(1-np.power(k1,2)))*special.ellipk(k)/(special.ellipk(k1)*special.ellipk(np.sqrt(1-np.power(k,2)))))
       for i in range(1,int(np.floor(self.n/2))+1):
            cd=mp.ellipfun('cd',(2*i-1)/self.n,k)
            zero = (1j/(float(k*cd.real)+1j*float(k*cd.imag)))
            self.num=self.num*np.poly1d([1/zero,1])
            self.num=self.num*np.poly1d([1/np.conj(zero),1])
            cd= mp.ellipfun('cd',(2*i-1)/self.n-1j*vo,k)
            pole=1j*(float(cd.real)+1j*float(cd.imag))
            if np.real(pole)<=0:    
                self.den=self.den*np.poly1d([-1/pole,1])
                self.den=self.den*np.poly1d([-1/np.conj(pole),1])
       if np.mod(self.n,2) == 1:
            sn= 1j*mp.ellipfun('sn',1j*vo,k)
            pole=1j*(float(sn.real)+1j*float(sn.imag))
            if np.real(pole)<=0:    
                self.den=self.den*np.poly1d([-1/pole,1])
                self.den=self.den*np.poly1d([-1/np.conj(pole),1])
       self.zeroes=np.roots(self.num)
       self.poles=np.roots(self.den)
       self.aprox_gain=np.power(gp,1-(self.n-2*np.floor(self.n/2)))
       self.num=self.num*self.aprox_gain
       self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained
                