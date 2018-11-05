from . import base_filter
import numpy as np
from scipy import signal

class Chebyshev(base_filter):
    #Chebyshev approximation class.In order to instantiate receives variadic arguments depending on filter type required.
    #If required filter is Low Pass or High Pass object receives in order:
    #string filter_type : Type of filter required to implement (LowPass,HighPass,BandPass,StopBand).
    #float gain : Desired DC gain.
    #float wp : Desired passband frequency.
    #float wa : Desired attenuation frequency.
    #float Ap : Desired passband maximum attenuation.
    #float Aa : Desired stopband minimum attenuation.
    #int n : Fixed approximation order. If no limitations regarding this value are met function should receive n=None, class will find smaller order that meets template requirements.
    #float max Q : Desired maximum pole Q factor. If no limitations regarding this value are met function should receive max_q=None. When Q is a limiting factor, it is not guaranteed that template requirements will be met.
    #float denorm Percent : Porcentual value, if equals 0 attenuation at w=wp will be Ap. If value increases, approximation will shift rigthwards until, at w=wa attenuation is Aa

    #If required filter is BandPass or StopBand function receives wp+ and wa+ after the previously defined wp and wa (this values will be used as wp- and wa-)


    #Filter initialization with initial parameters received
    def __init__(self, *args):
        if (args[0]=='LowPass')|(args[0]=='HighPass'):
            self.name = args[0]
            self.gain=args[1]
            self.wpl=args[2]*2*np.pi
            self.wal=args[3]*2*np.pi
            self.Ap=args[4]
            self.Ao=args[5]
            self.n=args[6]
            self.input_qmax=args[7]
            self.denorm_percent=args[8]/100
        elif (args[0]=='BandPass')|(args[0]=='StopBand'):
            self.name = args[0]
            self.gain=args[1]
            self.wpl=args[2]*2*np.pi
            self.wph=args[3]*2*np.pi
            self.wal=args[4]*2*np.pi
            self.wah=args[5]*2*np.pi
            self.Ap=args[6]
            self.Ao=args[7]
            self.n=args[8]
            self.input_qmax=args[9]
            self.denorm_percent=args[10]/100
        if self.name:
            if self.check_4_infs_and_nans(args) is False:
                if (self.name =='LowPass') or (self.name =='HighPass'):
                    self.nmax=25
                else:
                    self.nmax=36
                self.errormsg=self.check_input()
                if self.n !=None:
                    self.set_fix_order()
                    if (self.name =='LowPass') or (self.name =='HighPass'):
                        self.nmax=25
                    else:
                        self.nmax=18
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
            else:
                self.errormsg='Input cannot be inf or NaN\n'

    def do_approximation(self):
       self.epsilon=np.sqrt(np.power(10,(self.Ap)/10)-1)
       if self.n ==None:
           self.n = int(np.ceil(np.arccosh(np.sqrt(np.power(10,(self.Ao)/10)-1)/self.epsilon)/np.arccosh(self.wan)))
       if (self.n>self.nmax) is False:
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
               self.num=self.num*1/np.power(10,self.Ap/20)
           else:
               self.aprox_gain=1
           
           self.zeroes=np.roots(self.num)
           self.poles=np.roots(self.den)
           self.den=np.real(self.den)
           self.num=np.real(self.num)
           self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained
           self.denormalize_range()
       else:
           self.errormsg='Required order surpasses maximum order limit for this approximation type.\n'

    def set_fix_order(self):
        if (self.name =='LowPass') or (self.name =='HighPass'):
            self.n=int(self.n)
        elif (self.name =='BandPass') or (self.name =='StopBand'):
            if np.mod(self.n,2)==0:
                self.n=int(self.n/2)
            else:
                self.errormsg=self.errormsg+'Required filter cannot be realized with the specified order\n'