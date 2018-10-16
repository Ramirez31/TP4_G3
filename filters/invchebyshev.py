from . import base_filter
import numpy as np
from scipy import signal

class Invchebyshev(base_filter):

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
       self.epsilon=1/np.sqrt(np.power(10,self.Ao/10)-1)
       self.n = int(np.ceil(np.arccosh(1/(self.epsilon*np.sqrt(np.power(10,self.Ap/10)-1)))/np.arccosh(self.wan)))
       for i in range(1,2*self.n+1):
           alfa=(2*i-1)*np.pi/(2*self.n)
           beta=np.absolute(np.arcsinh(1/self.epsilon)/self.n)
           pole=self.wan/(np.sin(alfa)*np.sinh(beta)+1j*np.cos(alfa)*np.cosh(beta))
           if np.real(pole)<=0:
                pol= np.poly1d([-1/pole, 1])
                self.den= self.den*pol
           if i<=self.n:
               if self.is_odd(self.n) and (i == (np.floor(self.n/2)+1)):
                   pass
               else:
                zero=self.wan*1j/np.cos(alfa)
                pol= np.poly1d([-1/zero, 1])
                self.num= self.num*pol
       self.zeroes=np.roots(self.num)
       self.poles=np.roots(self.den) 
       self.aprox_gain=1
       self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained 

    def is_odd(self,num):
        return num & 0x1