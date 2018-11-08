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
            if(self.wpl*self.wph==self.wal*self.wah):
                self.wan=(self.wah-self.wal)/(self.wph-self.wpl)
            else:
                tempwah=self.wpl*self.wph/self.wal
                tempwal=self.wpl*self.wph/self.wah
                if (tempwah<self.wah) and (tempwal>self.wal):
                    self.wan=(tempwah-tempwal)/(self.wph-self.wpl)
                elif(tempwah>self.wah) and (tempwal>self.wal):
                    self.wan=(self.wah-tempwal)/(self.wph-self.wpl)
                elif(tempwah<self.wah) and (tempwal<self.wal):
                    self.wan=(tempwah-self.wal)/(self.wph-self.wpl)

        elif self.name=='StopBand':
            if(self.wpl*self.wph==self.wal*self.wah):
                self.wan=(self.wph-self.wpl)/(self.wah-self.wal)
            else:
                tempwph=self.wal*self.wah/self.wpl
                tempwpl=self.wal*self.wah/self.wph
                if (tempwph<self.wph) and (tempwpl>self.wpl):
                    self.wan=(tempwph-tempwpl)/(self.wah-self.wal)
                elif (tempwph>self.wph) and (tempwpl>self.wpl):
                    self.wan=(self.wph-tempwpl)/(self.wah-self.wal)
                elif (tempwph<self.wph) and (tempwpl<self.wpl):
                    self.wan=(tempwph-self.wpl)/(self.wah-self.wal)
        elif self.name=='Group Delay':
            self.wrgn=self.wrg*self.tao0

    # Function denormalizes the approximation realized previously, it can denormalize to: LP, HP, BP, BS.
    def denormalize(self):
        self.den=np.poly1d([1])
        self.num=np.poly1d([1])
        tempPoles=[]
        tempZeros=[]

        if self.name=='LowPass': #If required filter is LP
            #Denormalize poles
            for i in range(0,len(self.poles)):#For each pole denormalization is realized
                if self.poles[i] != 0:#If pole is not zero
                    pol=np.poly1d([-1/(self.wpl*self.poles[i]),1])
                    if np.real(np.roots(pol))<=0:
                        self.den= self.den*pol #Filter is denormalized by frequency scaling it S=Sn/wc
                        tempPoles.append(np.roots(pol)[0])
                else:
                    pol=np.poly1d([1/(self.wpl),0])
                    self.den=self.den*pol
                    tempPoles.append(np.roots(pol)[0])
            #Denormalize zeroes
            for i in range(0,len(self.zeroes)):#For each zero denormalization is realized
                if self.zeroes[i] != 0:#If zero is not located in origin
                    pol=np.poly1d([1/(self.wpl*self.zeroes[i]),1])
                    self.num= self.num*pol #Filter is denormalized by frequency scaling it S=Sn/wc
                    tempZeros.append(np.roots(pol)[0])
                else:
                    pol=np.poly1d([1/(self.wpl),0])
                    self.num=self.num*pol
                    tempZeros.append(np.roots(pol)[0])

        elif self.name=='HighPass': #If required filter is HP
            #Denormalize poles
            self.num=np.poly1d(np.zeros(len(self.poles)),r=True) #Zeroes created after doing HP denormalization to poles
            for i in range(0,len(np.roots(self.num))):
                tempZeros.append(np.roots(self.num)[i])
            for i in range(0,len(self.poles)):
                if self.poles[i] != 0:#If pole is not zero
                    pol=np.poly1d([1,-self.wpl/(self.poles[i])])
                    if np.real(np.roots(pol))<=0:
                        self.den= self.den*pol #Filter is denormalized by frequency scaling it S=wc/Sn
                        tempPoles.append(np.roots(pol)[0])
                else:
                    pol=np.poly1d([-self.wpl])
                    self.den= self.den*pol #Filter is denormalized by frequency scaling it S=wc/Sn
                    #tempPoles.append(np.roots(pol)[0])
            #Denormalize zeroes
            pol=np.poly1d(np.zeros(len(self.zeroes)),r=True)
            self.den=self.den*pol #Poles created after doing HP denormalization to zeroes
            for i in range(0,len(np.roots(pol))):
                tempPoles.append(np.roots(pol)[i])
            for i in range(0,len(self.zeroes)):
                if self.zeroes[i] != 0:#If zero is not located in origin
                    pol=np.poly1d([1,self.wpl/(self.zeroes[i])])
                    self.num= self.num*pol #Filter is denormalized by frequency scaling it S=wc/Sn
                    tempZeros.append(np.roots(pol)[0])
                else:
                    pol=np.poly1d([self.wpl])
                    self.num= self.num*pol #Filter is denormalized by frequency scaling it S=wc/Sn
                    tempZeros.append(np.roots(pol)[0])

        elif self.name=='BandPass': #If required filter is BP
            wo=np.sqrt(self.wpl*self.wph)
            B=(self.wph-self.wpl)/wo
            #Denormalize poles
            self.num=np.poly1d(np.zeros(len(self.poles)),r=True) #Zeroes created after doing HP denormalization to poles
            for i in range(0,len(np.roots(self.num))):
                tempZeros.append(np.roots(self.num)[i])
            for i in range (0,len(self.poles)):
                if self.poles[i] != 0:#If pole is not zero
                    if np.absolute(np.real(self.poles[i]))<0.0000001:
                        self.poles[i]=1j*np.imag(self.poles[i])
                    pol=np.poly1d([-1/(wo*B*self.poles[i]), 1, -wo/(B*self.poles[i])])
                    self.den=self.den*pol
                    tempPoles.append(np.roots(pol)[0])
                    tempPoles.append(np.roots(pol)[1])
                else:
                    pol=np.poly1d([-1/(wo*B), 0, -wo/B])
                    self.den=self.den*pol
                    tempPoles.append(np.roots(pol)[0])
                    tempPoles.append(np.roots(pol)[1])
            #Denormalize zeroes
            pol=np.poly1d(np.zeros(len(self.zeroes)),r=True)
            self.den=self.den*pol #Zeroes created after doing HP denormalization to poles
            for i in range(0,len(np.roots(pol))):
                tempPoles.append(np.roots(pol)[i])
            for i in range (0,len(self.zeroes)):
                if self.zeroes[i] != 0:#If zero is not located in origin
                    if np.absolute(np.real(self.zeroes[i]))<0.0000001:
                        self.zeroes[i]=1j*np.imag(self.zeroes[i])
                    pol=np.poly1d([1/(wo*B*self.zeroes[i]), 1, wo/(B*self.zeroes[i])])
                    self.num=self.num*pol
                    tempZeros.append(np.roots(pol)[0])
                    tempZeros.append(np.roots(pol)[1])
                else:
                    pol=np.poly1d([1/(wo*B), 0, wo/B])
                    self.num=self.num*pol
                    tempZeros.append(np.roots(pol)[0])
                    tempZeros.append(np.roots(pol)[1])

        elif self.name=='StopBand': #If required filter is SB
            wo=np.sqrt(self.wal*self.wah)
            tempwph=self.wal*self.wah/self.wpl
            tempwpl=self.wal*self.wah/self.wph
            if (tempwph<self.wph) and (tempwpl>self.wpl):
                B=(tempwph-tempwpl)/wo
            elif (tempwph>self.wph) and (tempwpl>self.wpl):
                B=(self.wph-tempwpl)/wo
            elif (tempwph<self.wph) and (tempwpl<self.wpl):
                B=(tempwph-self.wpl)/wo            
            else:
                B=(self.wph-self.wpl)/wo
            bool=True
            #Denormalize poles
            for i in range (0,len(self.poles)):
                if self.poles[i] != 0:#If pole is not zero
                    if np.absolute(np.real(self.poles[i]))<0.0001:
                        self.poles[i]=1j*np.imag(self.poles[i])
                    pol=np.poly1d([1/wo, B/(-self.poles[i]), wo])
                    self.den=self.den*pol
                    tempPoles.append(np.roots(pol)[0])
                    tempPoles.append(np.roots(pol)[1])
                else:
                    pol=np.poly1d([B,0])
                    self.den=self.den*pol
                    tempPoles.append(np.roots(pol)[0])
                if (len(self.zeroes)!=0) and np.mod(self.n,2)==0:
                    bool=False
                if bool==True: 
                    pol=np.poly1d([1/wo,0,wo])
                    self.num=self.num*pol
                    tempZeros.append(np.roots(pol)[0])
                    tempZeros.append(np.roots(pol)[1])
                if len(self.zeroes)!=0:
                    bool=False

            #Denormalize zeroes
            for i in range (0,len(self.zeroes)):
                if self.zeroes[i] != 0:#If pole is not zero
                    pol=np.poly1d([1/wo, B/(self.zeroes[i]), wo])
                    self.num=self.num*pol
                    tempZeros.append(np.roots(pol)[0])
                    tempZeros.append(np.roots(pol)[1])
                else:
                    pol=np.poly1d([B,0])
                    self.num=self.num*pol
                    tempZeros.append(np.roots(pol)[0])

        elif self.name=='Group Delay':
            #Denormalize poles
            for i in range(0,len(self.poles)):#For each pole denormalization is realized
                if self.poles[i] != 0:#If pole is not zero
                    pol=np.poly1d([-self.tao0/(self.poles[i]),1])
                    self.den= self.den*pol #Filter is denormalized by frequency scaling it S=Sn/wc
                    tempPoles.append(np.roots(pol)[0])
                else:
                    pol=np.poly1d([self.tao0,0])
                    self.den=self.den*pol
                    tempPoles.append(np.roots(pol)[0])
            #Denormalize zeroes
            for i in range(0,len(self.zeroes)):#For each zero denormalization is realized
                if self.zeroes[i] != 0:#If zero is not located in origin
                    pol=np.poly1d([self.tao0/self.zeroes[i],1])
                    self.num= self.num*pol #Filter is denormalized by frequency scaling it S=Sn/wc
                    tempZeros.append(np.roots(pol)[0])
                else:
                    pol=np.poly1d([self.tao0,0])
                    self.num=self.num*pol
                    tempZeros.append(np.roots(pol)[0])
        else:
            pass
        self.zeroes=tempZeros
        self.poles=tempPoles
        self.den=np.real(self.den)
        self.num=np.real(self.num)

        for i in range(0,len(self.poles)):
            if np.real(self.poles[i])>0:
                self.poles[i]=-np.real(self.poles[i])+1j*np.imag(self.poles[i])
        for sing_pole in self.poles:
            if np.real(sing_pole) !=0:
                if (-np.absolute(sing_pole)/(2*np.real(sing_pole)))>self.q:
                    self.q=(-np.absolute(sing_pole)/(2*np.real(sing_pole)))
            else:
                self.q=float('inf')
        self.K=np.power(10,self.gain/20)*self.aprox_gain
        self.num=self.num*self.K
        self.denorm_n=len(self.den)-1
        #if self.name== 'StopBand':
        #    self.denorm_sys=signal.ZerosPolesGain.to_tf(signal.ZerosPolesGain(self.zeroes,self.poles,self.K))
        #else:
        self.denorm_sys = signal.TransferFunction(self.num,self.den) #Denormalized system is obtained

    # Function returns current denormalized filter step response
    def get_step(self):
        a=signal.step(self.denorm_sys)
        return signal.step2(self.denorm_sys)

    # Function returns current denormalized filter impulse response
    def get_impulse(self):
        return signal.impulse2(self.denorm_sys)

    # Function returns current filter frequency response (frec,magnitude,phase)
    def get_bode(self):
        if self.name=='Group Delay':
            self.w,self.mag,self.phase = signal.bode(self.denorm_sys)
            return self.w, self.mag, self.phase
        if self.name=='LowPass': #If required filter is LP
            w=np.hstack((np.logspace(np.log10(self.wpl/1000000),np.log10(self.wpl),num=1000),np.logspace(np.log10(self.wpl*1.01),np.log10(self.wal*0.99),num=1000)))
            w=np.hstack((w,np.logspace(np.log10(self.wal),np.log10(self.wal*1000000),num=1000)))
        elif self.name=='HighPass': #If required filter is HP
            w=np.hstack((np.logspace(np.log10(self.wal/1000000),np.log10(self.wal),num=1000),np.logspace(np.log10(self.wal*1.01),np.log10(self.wpl*0.99),num=1000)))
            w=np.hstack((w,np.logspace(np.log10(self.wpl),np.log10(self.wpl*1000000),num=1000)))
        elif self.name=='BandPass': #If required filter is BP
            w=np.hstack((np.logspace(np.log10(self.wal/1000000),np.log10(self.wal),num=1000),np.logspace(np.log10(self.wal*1.01),np.log10(self.wpl*0.99),num=1000)))
            w=np.hstack((w,np.logspace(np.log10(self.wpl),np.log10(self.wph),num=1000)))
            w=np.hstack((w,np.logspace(np.log10(self.wph*1.01),np.log10(self.wah),num=1000)))
            w=np.hstack((w,np.logspace(np.log10(self.wah*1.01),np.log10(self.wah*1000000),num=1000)))
        elif self.name=='StopBand': #If required filter is SB
            w=np.hstack((np.logspace(np.log10(self.wpl/1000000),np.log10(self.wpl),num=1000),np.logspace(np.log10(self.wpl*1.01),np.log10(self.wal*0.99),num=1000)))
            w=np.hstack((w,np.logspace(np.log10(self.wal),np.log10(self.wah),num=1000)))
            w=np.hstack((w,np.logspace(np.log10(self.wah*1.01),np.log10(self.wph),num=1000)))
            w=np.hstack((w,np.logspace(np.log10(self.wph*1.01),np.log10(self.wph*1000000),num=1000)))
        self.w,self.mag,self.phase = signal.bode(self.denorm_sys,w)
        return self.w, self.mag, self.phase

    # Function returns current filter normalized frequency response (frec,magnitude,phase)
    def get_norm_bode(self):
        if self.name!='Group Delay':
            w=np.hstack((np.logspace(-6,0,num=1000),np.logspace(0.01,np.log10(self.wan*0.99),num=1000)))
            w=np.hstack((w,np.logspace(np.log10(self.wan),np.log10(self.wan*100000),num=500)))
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
            return [self.Ap,self.Ao,self.wan/(2*np.pi),self.wpl/(2*np.pi),self.wal/(2*np.pi),self.gain]
        if (self.name=='BandPass')|(self.name=='StopBand'):
            return [self.Ap,self.Ao,self.wan/(2*np.pi),self.wpl/(2*np.pi),self.wph/(2*np.pi),self.wal/(2*np.pi),self.wah/(2*np.pi),self.gain]
        if (self.name=='Group Delay'):
            return [self.tao0,self.wrg/(2*np.pi),self.palm,self.gain]
        

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
            max_denorm_w=(self.wan/w[i])*0.9999#In logarithmic scale wan/frec previously found indicates the maximum frequency that can be scaled
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
        return self.K

    def check_4_infs_and_nans(*args):
        for input in args[1]:
            if (input != input) or (input==float('inf')):
                return True
        return False