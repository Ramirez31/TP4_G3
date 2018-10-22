from . import base_filter
import numpy as np
from scipy import signal

class Chebyshev(base_filter):

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

    def do_approximation(self):
       self.epsilon=np.sqrt(np.power(10,(self.Ap)/10)-1)
       if self.n ==None:
           self.n = int(np.ceil(np.arccosh(np.sqrt(np.power(10,(self.Ao)/10)-1)/self.epsilon)/np.arccosh(self.wan)))
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
       if np.mod(self.n,2) == 0:
           self.aprox_gain=1/np.power(10,self.Ap/20)
           self.num=self.num*self.aprox_gain
       else:
           self.aprox_gain=1
       
       self.zeroes=np.roots(self.num)
       self.poles=np.roots(self.den)
       self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained