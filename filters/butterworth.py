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
            self.denorm_percent=args[8]/100
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
            self.denorm_percent=args[10]/100
        if self.name:
            if self.name == 'BandPass':#Order limit for fix order (Max denormalized order)
                self.nmax=18
            elif self.name == 'StopBand':
                self.nmax=10
            else:
                self.nmax=20
            self.errormsg=self.check_input()
            if self.n !=None:
                self.set_fix_order()
            if self.name == 'BandPass':#Order limit for normalized approximation (taking into account that denormalized limits are not met)
                self.nmax=9
            elif self.name == 'StopBand':
                self.nmax=10
            else:
                self.nmax=20
            if self.errormsg == '':
                while True:
                    self.q=0
                    self.poles=[]
                    self.zeroes=[]
                    self.den=np.poly1d([1])
                    self.num=np.poly1d([1])
                    self.normalize()#Normalizes current template
                    self.do_approximation()#Does normalized approximation and realizes
                    if self.errormsg =='':
                        self.denormalize()#Denormalizes the approximation to match desired template

                    if(self.input_qmax==None) or (self.input_qmax>=self.q):
                        break
                    else:
                        if self.n>1:
                            self.n=self.n-1
                        else:
                            self.errormsg=self.errormsg+'Required Q factor can not be achieved with current template\n'
                            break

    #
    def do_approximation(self):
        self.epsilon=np.sqrt(np.power(10,(self.Ap)/10)-1)
        self.maxq=0
        if self.n==None:
            self.n = int(np.ceil(np.log10((np.power(10,(self.Ao)/10)-1)/np.power(self.epsilon,2))/(2*np.log10(self.wan))))
        if (self.n>self.nmax) is False:
            ro=np.power(self.epsilon,-1/self.n)
            for i in range(1,int(self.n)+1):
                root= ro*(-np.sin(np.pi*(2*i-1)/(2*self.n))+1j*np.cos(np.pi*(2*i-1)/(2*self.n)))
                if np.real(root)<=0: #Only left-plane poles are utilized
                    self.poles.append(root)
                    pol= np.poly1d([-1/root, 1])
                    self.den= self.den*pol
            self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained
            self.denormalize_range()
            self.aprox_gain=1
        else:
            self.errormsg='Required order surpasses maximum order limit for this approximation type.\n'

    def set_fix_order(self):
        if (self.name =='LowPass') or (self.name =='HighPass'):
            self.n=self.n
        elif (self.name =='BandPass') or (self.name =='StopBand'):
            if np.mod(self.n,2)==0:
                self.n=int(self.n/2)
            else:
                self.errormsg=self.errormsg+'Required filter cannot be realized with the specified order\n'
