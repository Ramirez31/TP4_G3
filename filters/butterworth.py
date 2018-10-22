from . import base_filter
import numpy as np
from scipy import signal

class Butterworth(base_filter):

    #Filter initialization with initial parameters received
    def __init__(self, *args):
        if (args[0]=='LowPass')|(args[0]=='HighPass'):
            self.name = args[0]
            self.gain=args[1]
            self.wpl=args[2]
            self.wal=args[3]
            self.Ap=args[4]
            self.Ao=args[5]
            self.n=args[6]
            self.input_qmax=args[7]
        elif (args[0]=='BandPass')|(args[0]=='StopBand'):
            self.name = args[0]
            self.gain=args[1]
            self.wpl=args[2]
            self.wph=args[3]
            self.wal=args[4]
            self.wah=args[5]
            self.Ap=args[6]
            self.Ao=args[7]
            self.n=args[8]
            self.input_qmax=args[9]
        if self.name:
            self.nmax=1000 #VALOR NO DEFINITIVO, PROBAR CUAL ES EL VALOR MAXIMO PARA EL QUE EMPIEZA A MORIR LA APROXIMACION
            self.errormsg=self.check_input()
            if self.errormsg == '':
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
        self.maxq=0
        if self.n==None:
            self.n = int(np.ceil(np.log10((np.power(10,(self.Ao)/10)-1)/np.power(self.epsilon,2))/(2*np.log10(self.wan))))
        if (self.n>self.nmax) is False:
            ro=np.power(self.epsilon,-1/self.n)
            for i in range(1,self.n+1):
                root= ro*(-np.sin(np.pi*(2*i-1)/(2*self.n))+1j*np.cos(np.pi*(2*i-1)/(2*self.n)))
                if np.real(root)<=0: #Only left-plane poles are utilized
                    self.poles.append(root)
                    pol= np.poly1d([-1/root, 1])
                    self.den= self.den*pol
            self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained
            self.aprox_gain=1
        else:
            self.errormsg='Required order surpasses maximum order limit for this approximation type.'