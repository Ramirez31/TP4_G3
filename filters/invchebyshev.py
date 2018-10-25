from . import base_filter
import numpy as np
from scipy import signal

class Invchebyshev(base_filter):
    #Inverse Chebyshev approximation class.In order to instantiate receives variadic arguments depending on filter type required.
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
            self.nmax=1000
            if self.check_4_infs_and_nans(args) is False:
                self.errormsg=self.check_input()
                if self.n !=None:
                    self.set_fix_order()
                if self.errormsg == '':
                    while True:
                        self.q=0
                        self.poles=[]
                        self.zeroes=[]
                        self.den=np.poly1d([1])
                        self.num=np.poly1d([1])
                        self.normalize()#Normalizes current template
                        self.do_approximation()#Does normalized approximation and realizes
                        self.denormalize()#Denormalizes the approximation to match desired template

                        if(self.input_qmax==None) or (self.input_qmax>=self.q):
                            break
                        else:
                            if self.n>2:
                                self.n=self.n-1
                            else:
                                self.errormsg=self.errormsg+'Required Q factor can not be achieved with current template\n'
                                break
            else:
                self.errormsg='Input cannot be inf or NaN\n'

    def do_approximation(self):
       self.epsilon=1/np.sqrt(np.power(10,self.Ao/10)-1)
       if self.n == None:
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
       self.denormalize_range()

    def is_odd(self,num):
        return num & 0x1

    def set_fix_order(self):#If an order is fixed for denormalized filter, n for normalized filter is obtained
        if (self.name =='LowPass') :
            self.n=self.n
        elif (self.name =='HighPass'):
            if np.mod(self.n,4)==0:
                self.n=int(self.n/2)
            elif np.mod(self.n-1,4)==0:
                self.n=int((self.n+1)/2)
            else:
                self.errormsg=self.errormsg+'Required filter cannot be realized with the specified order\n'
        elif (self.name=='BandPass'):
            if np.mod(self.n,6)==0:
                self.n=int(self.n/3)
            elif np.mod(self.n-2,6)==0:
                self.n=int((self.n+1)/3)
            else:
                self.errormsg=self.errormsg+'Required filter cannot be realized with the specified order\n'
        elif (self.name=='StopBand'):
            if np.mod(self.n,8)==0:
                self.n=int(self.n/4)
            elif np.mod(self.n-2,8)==0:
                self.n=int((self.n+2)/4)
            else:
                self.errormsg=self.errormsg+'Required filter cannot be realized with the specified order\n'


    def denormalize_range(self):

        if self.denorm_percent != 0:
            w=np.linspace(1,self.wan,100000)
            w,mag,phase=signal.bode(self.norm_sys,w)
            for i in range(0,len(mag)):
                a=np.absolute(-mag[i]-self.Ap)
                if np.absolute(-mag[i]-self.Ap)<0.001:
                    break
            max_denorm_w=w[i]
            diff=max_denorm_w-1
            denorm_w=max_denorm_w-diff*(1-self.denorm_percent)
            den=np.poly1d([1])
            num=np.poly1d([1])
            for pole in self.poles:
                den=den*np.poly1d([-denorm_w/(pole),1])
            for zero in self.zeroes:
                num=num*np.poly1d([denorm_w/(zero),1])
            self.den=den
            self.poles=np.roots(self.den)
            self.num=num
            self.zeroes=np.roots(self.num)
            self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained