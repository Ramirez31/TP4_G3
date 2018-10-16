from . import base_filter
import numpy as np
from scipy import signal

class Butterworth(base_filter):

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

    #
    def do_approximation(self):
        self.epsilon=np.sqrt(np.power(10,(self.Ap)/10)-1)
        self.n = int(np.ceil(np.log10((np.power(10,(self.Ao)/10)-1)/np.power(self.epsilon,2))/(2*np.log10(self.wan))))
        ro=np.power(self.epsilon,-1/self.n)
        for i in range(1,self.n+1):
            root= ro*(-np.sin(np.pi*(2*i-1)/(2*self.n))+1j*np.cos(np.pi*(2*i-1)/(2*self.n)))
            if np.real(root)<=0: #Only left-plane poles are utilized
                self.poles.append(root)
                pol= np.poly1d([-1/root, 1])
                self.den= self.den*pol
        self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained
        self.aprox_gain=1