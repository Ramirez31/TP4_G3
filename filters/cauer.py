from . import base_filter
import numpy as np
from scipy import signal
from scipy import special

class Cauer(base_filter):

    #Filter initialization with initial parameters received
    def __init__(self, name,Ap,Ao,wpl,wph,wal,wah,gain,n):
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
            self.poles=[]
            self.zeroes=[]
            self.den=np.poly1d([1])
            self.num=np.poly1d([1])
            self.normalize()#Normalizes current template
            self.do_approximation()#Does normalized approximation and realizes
            self.denormalize()#Denormalizes the approximation to match desired template

    def do_approximation(self):
       epsilonp=1/np.sqrt(np.power(10,self.Ap/10)-1)
       k1=np.sqrt((np.power(10,self.Ap/10)-1)/(np.power(10,self.Ao/10)-1))
       k=1/self.wan
       a=np.arcsin(1j/epsilonp)
       vo=special.ellipkinc(a,k1)
       self.n=np.ceil(special.ellipk(np.sqrt(1-np.power(k1,2)))*special.ellipk(k)/(special.ellipk(k1)*special.ellipk(np.sqrt(1-np.power(k,2)))))
       for i in range(1,int(np.floor(self.n/2))+1):
            sn,cn,dn=special.ellipj((2*i-1)/self.n,k)
            cd=cn/dn
            zero = 1j/(k*cd)
            self.num=self.num*np.poly1d([1/zero,1])
            vo=special.ellipkinc(np.arcsin(1j/epsilonp),k1)
            sn1,cn1,dn1 = special.ellipj(((2*i-1)/self.n)-1j*vo,k)
            wi=cn1/dn1
            pole=1j*wi
            if np.real(pole)<=0:    
                self.den=self.den*np.poly1d([-1/pole,1])

       self.zeroes=np.roots(self.num)
       self.poles=np.roots(self.den)
       self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained
                