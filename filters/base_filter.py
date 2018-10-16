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
    def __init__(self, name=None,Ap=None,Ao=None,wpl=None,wph=None,wal=None,wah=None,gain=None,n=None,tao0=None,wrg=None,palm=None):
        if name:
            self.name = name
            self.Ap=Ap
            self.Ao=Ao
            self.wpl=wpl
            self.wph=wph
            self.wal=wal
            self._wah=wah
            self.n=n
            self.gain=gain
            self.tao0=tao0
            self.wrg=wrg
            self.palm=palm
            self.zeros=None
            self.poles=None
            self.system=None

    # Function does normalized approximation       
    @abstractmethod
    def do_approximation(self):
        pass

    # Function does normalized approximation, and denormalizes it with initial parameters provided (used for instanciating object)       
    def normalize(self):
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
        self.w,self.mag,self.phase = signal.bode(self.denorm_sys,n=1000)
        return self.w, self.mag, self.phase

    # Function returns current filter normalized frequency response (frec,magnitude,phase)
    def get_norm_bode(self):
        self.w,self.mag,self.nphase = signal.bode(self.norm_sys,n=1000)
        return self.w, self.mag, self.phase

    # Function returns current filter group delay
    def get_group_delay(self):
        dphase=np.ediff1d(self.phase)
        dw=np.ediff1d(self.w)
        gd=np.append(-dphase/dw,0)
        return gd
        #return -np.gradient(self.phase,self.w)

    # Function returns current filter zeroes and poles (zeroes, poles)
    def get_zeroes_poles(self):
        return self.zeroes, self.poles

    # Function returns current filter template limitations
    def get_template(self):
        return self.Ap,self.Ao,self.wpl,self.wph,self.wal,self.wah,self.wan,self.gain

    # Function returns current filter type:LP,HP,BP or SB
    def filter_is(self):
        return self.name


