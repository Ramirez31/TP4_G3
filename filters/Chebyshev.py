from . import base_filter
import numpy as np
from scipy import signal

class Chebyshev(base_filter):

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
       self.epsilon=np.sqrt(np.power(10,self.Ap/10)-1)
       self.n = int(np.ceil(np.arccosh(np.sqrt(np.power(10,self.Ao/10)-1)/self.epsilon)/np.arccosh(self.wan)))
       for i in range(1,self.n+1):
           alfa=(2*i-1)*np.pi/(2*self.n)
           beta1= np.absolute(np.arcsinh(1/self.epsilon)/self.n)
           beta2=-np.absolute(np.arcsinh(1/self.epsilon)/self.n) 
           pole1= np.sin(alfa)*np.sinh(beta1)+1j*np.cos(alfa)*np.cosh(beta1)
           pole2= np.sin(alfa)*np.sinh(beta2)+1j*np.cos(alfa)*np.cosh(beta2)
           if np.real(pole1)<=0:
                pol= np.poly1d([-1/pole1, 1])
                self.den= self.den*pol
           if np.real(pole2)<=0:
                pol= np.poly1d([-1/pole2, 1])
                self.den= self.den*pol
       self.zeroes=np.roots(self.num)
       self.poles=np.roots(self.den) 
       self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained