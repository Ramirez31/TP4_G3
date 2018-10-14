from . import base_filter
import numpy as np
from scipy import signal

class Invchebyshev(base_filter):

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
       self.epsilon=1/np.sqrt(np.power(10,self.Ao/10)-1)
       self.n = int(np.ceil(np.arccosh(1/(self.epsilon*np.sqrt(np.power(10,self.Ap/10)-1)))/np.arccosh(self.wan)))
       for i in range(1,self.n+1):
            alfa=(2*i-1)*np.pi/(2*self.n)
            beta1=np.arcsinh(1/self.epsilon)/self.n
            beta2=-np.arcsinh(1/self.epsilon)/self.n
            pole1=self.wan*(1/(np.sin(alfa)*np.sinh(beta1)+1j*np.cos(alfa)*np.cosh(beta1)))
            pole2=self.wan*(1/(np.sin(alfa)*np.sinh(beta2)+1j*np.cos(alfa)*np.cosh(beta2)))
            if(np.real(pole1)<=0):
                self.poles.append(pole1)
                pol= np.poly1d([-1/pole1, 1])
                self.den= self.den*pol
            if(np.real(pole2)<=0):
                self.poles.append(pole2)
                pol= np.poly1d([-1/pole2, 1])
                self.den= self.den*pol
            zero1=self.wan*1j/np.cos(alfa)
            zero2=-self.wan*1j/np.cos(alfa)
            if zero1 not in self.zeroes:
                if(self.is_odd(self.n) and (np.floor(self.n/2)+1 ==i)):
                    self.den=self.den*np.poly1d([1, 0])
                    self.poles.append(0)
                    bool=True
                else:
                    self.zeroes.append(zero1)
                    pol= np.poly1d([-1/zero1, 1])
                    self.num= self.num*pol
            if zero2 not in self.zeroes:
                    self.zeroes.append(zero2)
                    pol= np.poly1d([-1/zero2, 1])
                    self.num= self.num*pol
       self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained

    def is_odd(self,num):
        return num & 0x1











