from abc import ABCMeta, abstractmethod
import numpy as np
from scipy import signal

class base_filter(metaclass=ABCMeta):
    #Base abstract filter class, contains attributes and methods common to all approximation types. Class receives variadic arguments depending on required approximation.

    #Function initializes class
    def __init__(self,*args,):
        pass

    # Function does normalized approximation       
    @abstractmethod
    def do_approximation(self):
        pass

    # Function does normalized approximation, and denormalizes it with initial parameters provided (used for instanciating object)       
    def normalize(self):
        if self.name != 'Group Delay':
            self.Ap=self.Ap+self.gain
            self.Ao=self.Ao+self.gain
        if self.name=='LowPass':
            self.wan=self.wal/self.wpl
        elif self.name=='HighPass':
            self.wan=self.wpl/self.wal
        elif self.name=='BandPass':
            self.wan=(self.wah-self.wal)/(self.wph-self.wpl)
        elif self.name=='StopBand':
            self.wan=(self.wph-self.wpl)/(self.wah-self.wal)
        elif self.name=='Group Delay':
            self.wrgn=self.wrg*self.tao0

    # Function denormalizes the approximation realized previously, it can denormalize to: LP, HP, BP, BS.
    def denormalize(self):
        self.den=np.poly1d([1])
        self.num=np.poly1d([1])

        if self.name=='LowPass': #If required filter is LP
            #Denormalize poles
            for i in range(0,len(self.poles)):#For each pole denormalization is realized
                if self.poles[i] != 0:#If pole is not zero
                    pol=np.poly1d([-1/(self.wpl*self.poles[i]),1])
                    if np.real(np.roots(pol))<=0:
                        self.den= self.den*pol #Filter is denormalized by frequency scaling it S=Sn/wc
                else:
                    self.den=self.den*np.poly1d([1/(self.wpl),0])
            #Denormalize zeroes
            for i in range(0,len(self.zeroes)):#For each zero denormalization is realized
                if self.zeroes[i] != 0:#If zero is not located in origin
                    self.num= self.num*np.poly1d([1/(self.wpl*self.zeroes[i]),1]) #Filter is denormalized by frequency scaling it S=Sn/wc
                else:
                    self.num=self.num*np.poly1d([1/(self.wpl),0])
            self.zeroes=np.roots(self.num)
            self.poles=np.roots(self.den)

        elif self.name=='HighPass': #If required filter is HP
            #Denormalize poles
            self.num=np.poly1d(np.zeros(len(self.poles)),r=True) #Zeroes created after doing HP denormalization to poles
            for i in range(0,len(self.poles)):
                if self.poles[i] != 0:#If pole is not zero
                    pol=np.poly1d([1,-self.wpl/(self.poles[i])])
                    if np.real(np.roots(pol))<=0:
                        self.den= self.den*pol #Filter is denormalized by frequency scaling it S=wc/Sn
                else:
                    self.den= self.den*np.poly1d([-self.wpl]) #Filter is denormalized by frequency scaling it S=wc/Sn
            #Denormalize zeroes
            self.den=self.den*np.poly1d(np.zeros(len(self.zeroes)),r=True) #Poles created after doing HP denormalization to zeroes
            for i in range(0,len(self.zeroes)):
                if self.zeroes[i] != 0:#If zero is not located in origin
                    self.num= self.num*np.poly1d([1,self.wpl/(self.zeroes[i])]) #Filter is denormalized by frequency scaling it S=wc/Sn
                else:
                    self.num= self.num*np.poly1d([self.wpl]) #Filter is denormalized by frequency scaling it S=wc/Sn
            self.zeroes=np.roots(self.num)
            self.poles=np.roots(self.den)

        elif self.name=='BandPass': #If required filter is BP
            wo=np.sqrt(self.wpl*self.wph)
            B=(self.wph-self.wpl)/wo
            #Denormalize poles
            self.num=np.poly1d(np.zeros(len(self.poles)),r=True) #Zeroes created after doing HP denormalization to poles
            for i in range (0,len(self.poles)):
                if self.poles[i] != 0:#If pole is not zero
                    if np.absolute(np.real(self.poles[i]))<0.01:
                        self.poles[i]=1j*np.imag(self.poles[i])
                    self.den=self.den*np.poly1d([-1/(wo*B*self.poles[i]), 1, -wo/(B*self.poles[i])])
                else:
                    self.den=self.den*np.poly1d([-1/(wo*B), 0, -wo/B])
            #Denormalize zeroes
            self.den=self.den*np.poly1d(np.zeros(len(self.zeroes)),r=True) #Zeroes created after doing HP denormalization to poles
            for i in range (0,len(self.zeroes)):
                if self.zeroes[i] != 0:#If zero is not located in origin
                    if np.absolute(np.real(self.zeroes[i]))<0.01:
                        self.zeroes[i]=1j*np.imag(self.zeroes[i])
                    self.num=self.num*np.poly1d([1/(wo*B*self.zeroes[i]), 1, wo/(B*self.zeroes[i])])
                else:
                    self.num=self.num*np.poly1d([1/(wo*B), 0, wo/B])
            self.poles=np.roots(self.den)
            self.zeroes=np.roots(self.num)

        elif self.name=='StopBand': #If required filter is SB
            wo=np.sqrt(self.wpl*self.wph)
            B=(self.wph-self.wpl)/wo
            #Denormalize poles
            for i in range (0,len(self.poles)):
                if self.poles[i] != 0:#If pole is not zero
                    if np.absolute(np.real(self.poles[i]))<0.01:
                        self.poles[i]=1j*np.imag(self.poles[i])
                    self.den=self.den*np.poly1d([1/wo, B/(-self.poles[i]), wo])
                else:
                    self.den=self.den*np.poly1d([B,0])
                self.num=self.num*np.poly1d([1/wo,0,wo])
            #Denormalize zeroes
            for i in range (0,len(self.zeroes)):
                if self.zeroes[i] != 0:#If pole is not zero
                    if np.absolute(np.real(self.zeroes[i]))<0.01:
                        self.zeroes[i]=1j*np.imag(self.zeroes[i])
                    self.num=self.num*np.poly1d([1/wo, B/(self.zeroes[i]), wo])
                else:
                    self.num=self.num*np.poly1d([B,0])
                self.den=self.den*np.poly1d([1/wo,0,wo])
            self.zeroes=np.roots(self.num)
            self.poles=np.roots(self.den)


        elif self.name=='Group Delay':
            #Denormalize poles
            for i in range(0,len(self.poles)):#For each pole denormalization is realized
                if self.poles[i] != 0:#If pole is not zero
                    self.den= self.den*np.poly1d([-self.tao0/(self.poles[i]),1]) #Filter is denormalized by frequency scaling it S=Sn/wc
                else:
                    self.den=self.den*np.poly1d([self.tao0,0])
            #Denormalize zeroes
            for i in range(0,len(self.zeroes)):#For each zero denormalization is realized
                if self.zeroes[i] != 0:#If zero is not located in origin
                    self.num= self.num*np.poly1d([self.tao0/self.zeroes[i],1]) #Filter is denormalized by frequency scaling it S=Sn/wc
                else:
                    self.num=self.num*np.poly1d([self.tao0,0])
            self.zeroes=np.roots(self.num)
            self.poles=np.roots(self.den)
        else:
            pass
        for sing_pole in self.poles:
            if np.real(sing_pole) !=0:
                if (-np.absolute(sing_pole)/(2*np.real(sing_pole)))>self.q:
                    self.q=(-np.absolute(sing_pole)/(2*np.real(sing_pole)))
            else:
                self.q=float('inf')
        K=np.power(10,self.gain/20)
        self.num=self.num*K*self.aprox_gain
        self.denorm_n=len(self.den)
        self.denorm_sys = signal.TransferFunction(self.num,self.den) #Denormalized system is obtained

    # Function returns current denormalized filter step response
    def get_step(self):
        return signal.step(self.denorm_sys,N=1000)

    # Function returns current denormalized filter impulse response
    def get_impulse(self):
        return signal.impulse(self.denorm_sys,N=1000)

    # Function returns current filter frequency response (frec,magnitude,phase)
    def get_bode(self):
        if self.name=='Group Delay':
            self.w,self.mag,self.phase = signal.bode(self.denorm_sys)
            return self.w, self.mag, self.phase
        if self.name=='LowPass': #If required filter is LP
            w=np.hstack((np.logspace(np.log10(self.wpl/1000),np.log10(self.wpl),num=500),np.logspace(np.log10(self.wpl*1.01),np.log10(self.wal*0.99),num=500)))
            w=np.hstack((w,np.logspace(np.log10(self.wal),np.log10(self.wal*1000),num=500)))
        elif self.name=='HighPass': #If required filter is HP
            w=np.hstack((np.logspace(np.log10(self.wal/1000),np.log10(self.wal),num=500),np.logspace(np.log10(self.wal*1.01),np.log10(self.wpl*0.99),num=500)))
            w=np.hstack((w,np.logspace(np.log10(self.wpl),np.log10(self.wpl*1000),num=500)))
        elif self.name=='BandPass': #If required filter is BP
            w=np.hstack((np.logspace(np.log10(self.wal/1000),np.log10(self.wal),num=500),np.logspace(np.log10(self.wal*1.01),np.log10(self.wpl*0.99),num=500)))
            w=np.hstack((w,np.logspace(np.log10(self.wpl),np.log10(self.wph),num=500)))
            w=np.hstack((w,np.logspace(np.log10(self.wph*1.01),np.log10(self.wah),num=500)))
            w=np.hstack((w,np.logspace(np.log10(self.wah*1.01),np.log10(self.wah*1000),num=500)))
        elif self.name=='StopBand': #If required filter is SB
            w=np.hstack((np.logspace(np.log10(self.wpl/1000),np.log10(self.wpl),num=500),np.logspace(np.log10(self.wpl*1.01),np.log10(self.wal*0.99),num=500)))
            w=np.hstack((w,np.logspace(np.log10(self.wal),np.log10(self.wah),num=500)))
            w=np.hstack((w,np.logspace(np.log10(self.wah*1.01),np.log10(self.wph),num=500)))
            w=np.hstack((w,np.logspace(np.log10(self.wph*1.01),np.log10(self.wph*1000),num=500)))
        self.w,self.mag,self.phase = signal.bode(self.denorm_sys,w)
        return self.w, self.mag, self.phase

    # Function returns current filter normalized frequency response (frec,magnitude,phase)
    def get_norm_bode(self):
        if self.name!='Group Delay':
            w=np.hstack((np.logspace(-3,0,num=500),np.logspace(0.01,np.log10(self.wan*0.99),num=500)))
            w=np.hstack((w,np.logspace(np.log10(self.wan),np.log10(self.wan*1000),num=500)))
            self.nw,self.nmag,self.nphase = signal.bode(self.norm_sys,w)
        else:
            self.nw,self.nmag,self.nphase = signal.bode(self.norm_sys)
        return self.nw,self.nmag,self.nphase

    # Function returns current filter group delay
    def get_group_delay(self):
        dphase=np.ediff1d(self.phase)#Phase is in deegres
        dw=np.ediff1d(self.w)
        gd=-dphase/dw
        gd=np.append(gd,gd[len(gd)-1])
        return gd/(180/np.pi)#Group delay is -dphi/dw, phi being in radians. Deegres are transformed to rads.
        #return -np.gradient(self.phase,self.w)

    # Function returns current filter zeroes and poles (zeroes, poles)
    def get_zeroes_poles(self):
        return self.zeroes, self.poles

    # Function returns current filter zeroes and poles (zeroes, poles)
    def get_poles(self):
        return self.poles

    # Function returns current filter zeroes and poles (zeroes, poles)
    def get_zeroes(self):
        return self.zeroes


    # Function returns current filter template limitations
    def get_template(self):
        if (self.name=='LowPass')|(self.name=='HighPass'):
            return [self.Ap,self.Ao,self.wan,self.wpl,self.wal,self.gain]
        if (self.name=='BandPass')|(self.name=='StopBand'):
            return [self.Ap,self.Ao,self.wan,self.wpl,self.wph,self.wal,self.wah,self.gain]
        if (self.name=='Group Delay'):
            return [self.tao0,self.wrg,self.palm,self.gain]
        

    # Function returns current filter type:LP,HP,BP or SB
    def filter_is(self):
        return self.name

    def check_input(self):
        errormsg=''
        if (self.gain >=0) and ((self.n==None) or (self.n<=self.nmax)) and((self.input_qmax==None) or (self.input_qmax>0)):
            if self.name == 'Group Delay':
                if (self.palm<=0) or (self.palm>=1):
                    errormsg=errormsg+'Error: Tolerance must be a real number higher than 0 and smaller than 1\n'
                if (errormsg =='') and (self.tao0<=0):
                    errormsg =errormsg+'Error: Group delay at 0 must be a positive value\n'
                if (errormsg =='') and (self.wrg<=0):
                    errormsg = errormsg +'Error:Wrg must be a positive real number\n'

            elif (self.Ao>0) and (self.Ap>0) and (self.Ao>=self.Ap) and (self.denorm_percent<=100) and (self.denorm_percent>=0):
                if self.name == 'LowPass':
                    if (self.wpl>0) and (self.wal>0):
                        if(self.wpl>=self.wal):
                            errormsg=errormsg+'Error: Template requiermentes not met, Wp must be smaller than Wa\n'
                    else:
                        errormsg=errormsg+'Error: Wp and Wa must be positive real values\n'
                elif self.name == 'HighPass':
                    if (self.wpl>0) and (self.wal>0):
                        if(self.wal>=self.wpl):
                            errormsg=errormsg+'Error: Template requirements not met, Wa must be smaller than Wp\n'
                    else:
                        errormsg=errormsg+'Error: Wp and Wa must be positive real values\n'
                elif self.name == 'BandPass':
                    if (self.wpl>0) and (self.wal>0) and (self.wah>0) and (self.wph>0):
                        if ((self.wah > self.wph) and (self.wph > self.wpl) and (self.wpl > self.wal) ) is False:
                            errormsg=errormsg+'Error: Template requirements not met.(Remember, Wa+>Wp+>Wp->Wa-)\n'
                    else:
                        errormsg=errormsg+'Error: Wp-, Wp+, Wa- and Wa+ must be positive real values\n'
                elif self.name == 'StopBand':
                    if (self.wpl>0) and (self.wal>0) and (self.wah>0) and (self.wph>0):
                        if ((self.wph > self.wah) and (self.wah > self.wal) and (self.wal > self.wpl)) is False:
                            errormsg=errormsg+'Error: Template requirements not met.(Remember, Wp+>Wa+>Wa->Wp-)\n'
                    else:
                        errormsg=errormsg+'Error: Wp-, Wp+, Wa- and Wa+ must be positive real values\n'
            else:
                if(self.Ao<0) or (self.Ap<0):
                    errormsg=errormsg+'Error: Ao and Ap must be positive valued real numbers\n'
                elif (self.Ao<self.Ap):
                    errormsg=errormsg+'Error: Ap must be smaller than Ao\n'
                elif (self.denorm_percent<0) or (self.denorm_percent>100):
                    errormsg=errormsg+'Error: Denormalization percent cannot be higher than 100% or smaller than 0%.\n'
                else:
                    errormsg=errormsg+'Error: Filter type was not found in our selection\n'
        else:
            if self.gain<0:
                errormsg=errormsg+'Error: Gain must be a positive real value.\n'
            if ((self.n!=None) and (self.n>self.nmax)):
                errormsg=errormsg+'Error: Fixed order surpasses maximum order limit for this approximation type.\n'
            if ((self.input_qmax!=None) and (self.input_qmax<=0)):
                errormsg=errormsg+'Error: Desired maximum Q factor must be a positive real number.\n'
        return errormsg

    def error_was(self):
        return self.errormsg

    def n_and_q_is(self):
        return self.denorm_n, self.q

    def denormalize_range(self):

        if self.denorm_percent != 0:
            w=np.linspace(1,self.wan,100000)
            w,mag,phase=signal.bode(self.norm_sys,w)
            for i in range(0,len(mag)):#Frequency where H(s)=Aa is found
                a=np.absolute(-mag[i]-self.Ao)
                if np.absolute(-mag[i]-self.Ao)<0.001:
                    break
            max_denorm_w=(self.wan/w[i])#In logarithmic scale wan/frec previously found indicates the maximum frequency that can be scaled
            diff=max_denorm_w-1#Difference between max denorm frequency and wp=1
            denorm_w=max_denorm_w-diff*(1-self.denorm_percent)#Depending on denormalization percent denorm frequency will vary betwen 1 and max_denorm
            den=np.poly1d([1])
            num=np.poly1d([1])
            for pole in self.poles:
                den=den*np.poly1d([-1/(denorm_w*pole),1])
            for zero in self.zeroes:
                num=num*np.poly1d([1/(denorm_w*zero),1])
            self.den=den
            self.poles=np.roots(self.den)
            self.num=num
            self.zeroes=np.roots(self.num)
            self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained
    
    def get_gain(self):
        return self.aprox_gain

    def check_4_infs_and_nans(*args):
        for input in args[1]:
            if (input != input) or (input==float('inf')):
                return True
        return False