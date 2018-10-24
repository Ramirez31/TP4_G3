from . import base_filter
import numpy as np
from scipy import signal

class Butterworth(base_filter):
    #Butterworth approximation class.In order to instantiate receives variadic arguments depending on filter type required.
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
    #
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
            if self.check_4_infs_and_nans(args) is False:
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
                        self.do_approximation()#Does normalized approximation
                        if self.errormsg =='':
                            self.denormalize()#Denormalizes the approximation to match desired template

                        if(self.input_qmax==None) or (self.input_qmax>=self.q):#If no Q limitations are received or if there was any, it was met, approximation is complete
                            break
                        else:
                            if self.n>1:# Filter order cannot be smaller than 1
                                self.n=self.n-1
                            else:
                                self.errormsg=self.errormsg+'Required Q factor can not be achieved with current template\n'
                                break
            else:
                self.errormsg='Input cannot be inf or NaN\n'

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
            self.denormalize_range()# Denormalize normalized template meeting required denormalization range.
            self.aprox_gain=1
        else:
            self.errormsg='Required order surpasses maximum order limit for this approximation type.\n'

    def set_fix_order(self):
        if (self.name =='LowPass') or (self.name =='HighPass'):
            self.n=self.n
        elif (self.name =='BandPass') or (self.name =='StopBand'):# If a K order SB or BP denormalized filter is required. LP normalized filter will have to have half its order.
            if np.mod(self.n,2)==0:
                self.n=int(self.n/2)
            else:
                self.errormsg=self.errormsg+'Required filter cannot be realized with the specified order\n'
