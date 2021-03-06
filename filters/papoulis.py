from . import base_filter
import numpy as np
from scipy import signal
from scipy import special

class Papoulis(base_filter):
    #Papoulis approximation class.In order to instantiate receives variadic arguments depending on filter type required.
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
                if (self.name=='LowPass') or (self.name=='HighPass'):
                    self.nmax=21
                else:
                    self.nmax=40
                self.errormsg=self.check_input()
                if self.n !=None:
                    self.set_fix_order()
                if (self.name=='LowPass') or (self.name=='HighPass'):
                    self.nmax=21
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
                        if self.errormsg == '':
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
       self.epsilon=np.sqrt(np.power(10,self.Ap/10)-1)
       if self.n == None:
            self.n=0
            Ln2_wan=0
            while (Ln2_wan >= ((np.power(10,self.Ao/10)-1)/np.power(self.epsilon,2)))==False: #While template specifications are not met, order increases
                self.n+=1
                if self.n>self.nmax:
                    self.errormsg='Required order surpasses maximum order limit for this approximation type.\n'
                    break
                Ln=np.poly1d([1])
                if np.mod(self.n,2)==1:#If n is odd 
                     k=(self.n-1)/2
                     ao=1/(np.sqrt(2)*(k+1))
                     v=np.poly1d([ao])
                     for r in range(1,int(k)+1):
                         ar=ao*(2*r+1)
                         vr=ar*special.legendre(r)
                         v=np.polyadd(v,vr)
                     v2=np.polymul(v,v)#Function to be integrated is calculated
                     Ln=np.polyint(v2)#Primitive is calculated, Barrow will be applied to this function
                else:#If n is even
                     k=(self.n-2)/2
                     p_legendre=special.legendre(k+1)
                     temp_pol=np.polymul(np.polyder(p_legendre),np.polyder(p_legendre))
                     phi=np.polymul(temp_pol,np.poly1d([1,1]))#Function to be integrated is calculated
                     Ln=np.polyint(phi)#Primitive is calculated, Barrow will be applied to this function with integration bounds(-S^2-1),-1
                
                Ln2=np.poly1d([0])
                poly2mul=np.poly1d([-2,0,-1])#First (-S^2-1) is replaced in the obtained primitive (upper bound Barrow)
                for i in range(0,len(Ln)+1):
                    polybeingmul=poly2mul
                    for j in range(0,i-1):#Loop used to get the n-th power of polymul
                         polybeingmul=np.polymul(polybeingmul,poly2mul)
                    Ln2=np.polyadd(Ln2,Ln[i]*polybeingmul)
                Ln2=np.polyadd(Ln2,-np.poly1d([np.polyval(Ln,-1)]))#Finally polynomial is evualated in -1 and this is substracted to previos value (lower bound Barrow)
                if np.mod(self.n,2)==0:#If n is even Ln2(w=1)=1, so knowing that S=wj->w=-jS
                    Ln2=Ln2/np.polyval(Ln2,-1j)#So as to make Ln2(w=1)=1
                Ln2_wan=np.polyval(Ln2,-1j*self.wan)#Ln2 is evaulated in Wan, so as to see if template condition is met
       else:
            Ln=np.poly1d([1])
            if np.mod(self.n,2)==1:#If n is odd 
                 k=(self.n-1)/2
                 ao=1/(np.sqrt(2)*(k+1))
                 v=np.poly1d([ao])
                 for r in range(1,int(k)+1):
                     ar=ao*(2*r+1)
                     vr=ar*special.legendre(r)
                     v=np.polyadd(v,vr)
                 v2=np.polymul(v,v)#Function to be integrated is calculated
                 Ln=np.polyint(v2)#Primitive is calculated, Barrow will be applied to this function
            else:#If n is even
                 k=(self.n-2)/2
                 p_legendre=special.legendre(k+1)
                 temp_pol=np.polymul(np.polyder(p_legendre),np.polyder(p_legendre))
                 phi=np.polymul(temp_pol,np.poly1d([1,1]))#Function to be integrated is calculated
                 Ln=np.polyint(phi)#Primitive is calculated, Barrow will be applied to this function with integration bounds(-S^2-1),-1
            
            Ln2=np.poly1d([0])
            poly2mul=np.poly1d([-2,0,-1])#First (-S^2-1) is replaced in the obtained primitive (upper bound Barrow)
            for i in range(0,len(Ln)+1):
                polybeingmul=poly2mul
                for j in range(0,i-1):#Loop used to get the n-th power of polymul
                     polybeingmul=np.polymul(polybeingmul,poly2mul)
                Ln2=np.polyadd(Ln2,Ln[i]*polybeingmul)
            Ln2=np.polyadd(Ln2,-np.poly1d([np.polyval(Ln,-1)]))#Finally polynomial is evualated in -1 and this is substracted to previos value (lower bound Barrow)
            if np.mod(self.n,2)==0:#If n is even Ln2(w=1)=1, so knowing that S=wj->w=-jS
                Ln2=Ln2/np.polyval(Ln2,-1j)#So as to make Ln2(w=1)=1
       roots=np.roots(np.polyadd(Ln2*np.power(self.epsilon,2),np.poly1d([1])))#Roots of the denominator are found
       for i in range(0,len(roots)):#If they are in the left plane, they are kept
           if np.real(roots[i])<=0:
                if roots[i]!=0:
                    pol=np.poly1d([-1/roots[i],1])#H(s) polynomial is constructed
                else:
                    pol=np.poly1d([1,0])
                self.den=self.den*pol
       self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained
       self.poles=np.roots(self.den)
       self.aprox_gain=1
       self.denormalize_range()

    def set_fix_order(self):
        if (self.name =='LowPass') or (self.name =='HighPass'):
            self.n=self.n
        elif (self.name =='BandPass') or (self.name =='StopBand'):
            if np.mod(self.n,2)==0:
                self.n=int(self.n/2)
            else:
                self.errormsg=self.errormsg+'Required filter cannot be realized with the specified order\n'