from abc import ABCMeta, abstractmethod
import numpy as np
from scipy import signal

class base_filter(metaclass=ABCMeta):
    """Base abstract filter class, contains attributes and methods common to all approximation types. In
    Ap=
    Ao=
    wpl=
    wph=
    wal=
    wah=
    wo=
    n= Filter order
    """
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
                    self.den= self.den*np.poly1d([-1/(self.wpl*self.poles[i]),1]) #Filter is denormalized by frequency scaling it S=Sn/wc
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
                    self.den= self.den*np.poly1d([1,-self.wpl/(self.poles[i])]) #Filter is denormalized by frequency scaling it S=wc/Sn
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
                    self.den=self.den*np.poly1d([-1/(wo*B*self.poles[i]), 1, -wo/(B*self.poles[i])])
                else:
                    self.den=self.den*np.poly1d([-1/(wo*B), 0, -wo/B])
            #Denormalize zeroes
            self.den=self.den*np.poly1d(np.zeros(len(self.zeroes)),r=True) #Zeroes created after doing HP denormalization to poles
            for i in range (0,len(self.zeroes)):
                if self.zeroes[i] != 0:#If zero is not located in origin
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
                    self.den=self.den*np.poly1d([1/wo, B/(-self.poles[i]), wo])
                else:
                    self.den=self.den*np.poly1d([B,0])
                self.num=self.num*np.poly1d([1/wo,0,wo])
            #Denormalize zeroes
            for i in range (0,len(self.zeroes)):
                if self.zeroes[i] != 0:#If pole is not zero
                    self.num=self.num*np.poly1d([1/wo, B/(self.zeroes[i]), wo])
                else:
                    self.num=self.num*np.poly1d([B,0])
                self.den=self.den*np.poly1d([1/wo,0,wo])
            self.poles=np.roots(self.den)
            self.zeroes=np.roots(self.num)

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
        K=np.power(10,self.gain/20)
        self.num=self.num*K*self.aprox_gain
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
        error=False
        if (self.gain>=0) or (self.n>self.nmax):
            if self.name == 'Group delay':
                if (self.palm<=0) or (self.palm>=1):
                    error = True
                if (error is False) and (self.tao0<=0):
                    error = True
                if (error is False) and (self.wrg<=0):
                    error = True

            elif (self.Ao>0) and (self.Ap>0) and (self.Ao>=self.Ap):
                if self.name == 'LowPass':
                    if(self.wpl>=self.wal):
                        error=True
                elif self.name == 'HighPass':
                    if(self.wal>=self.wpl):
                        error=True
                elif self.name == 'BandPass':
                    if ((self.wah > self.wph) and (self.wph > self.wpl) and (self.wpl > self.wal) ) is False:
                        error=True
                elif self.name == 'StopBand':
                    if ((self.wph > self.wah) and (self.wah > self.wal) and (self.wal > self.wpl)) is False:
                        error=True
            else:
                error=True
        else:
            error = True
        return error

    def error_was(self):
        return self.error


