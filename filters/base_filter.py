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
    def __init__(self, name=None,Ap=None,Ao=None,wpl=None,wph=None,wal=None,wah=None,gain=None,n=None):
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

    # Function denormalizes the approximation realized previously, it can denormalize to: LP, HP, BP, BS.
    def denormalize(self):
        self.den=np.poly1d([1])
        self.num=np.poly1d([1])

        if self.name=='LowPass': #If required filter is LP
            for i in range(0,len(self.poles)):
                self.den= self.den*np.poly1d([-1/(self.wpl*self.poles[i]),1]) #Filter is denormalized by frequency scaling it S=Sn/wc
            self.poles =[i * self.wpl for i in self.poles]

        elif self.name=='HighPass': #If required filter is HP
            self.zeroes = np.zeros(self.n)
            self.num=np.poly1d(self.zeroes,r=True)
            for i in range(0,len(self.poles)):
                self.den= self.den*np.poly1d([1,-self.wpl/(self.poles[i])]) #Filter is denormalized by frequency scaling it S=wc/Sn
            self.poles =[self.wpl/i for i in self.poles]

        elif self.name=='BandPass': #If required filter is BP
            self.zeroes = np.zeros(self.n)
            self.num=np.poly1d(self.zeroes,r=True)
            wo=np.sqrt(self.wpl*self.wph)
            B=(self.wph-self.wpl)/wo
            for i in range (0,len(self.poles)):
                self.den=self.den*np.poly1d([-1/(wo*B*self.poles[i]), 1, -wo/(B*self.poles[i])])
            self.poles=np.roots(self.den)

        elif self.name=='StopBand': #If required filter is SB
            wo=np.sqrt(self.wpl*self.wph)
            B=(self.wph-self.wpl)/wo
            for i in range (0,len(self.poles)):
                self.den=self.den*np.poly1d([1/wo, B/(-self.poles[i]), wo])
                self.num=self.num*np.poly1d([1/wo,0,wo])
            self.poles=np.roots(self.den)
            self.zeroes=np.roots(self.num)

        else:
            pass
        K=np.power(10,self.gain/20)
        self.denorm_sys = signal.TransferFunction(self.num,self.den) #Denormalized system is obtained

    # Function returns current denormalized filter step response
    @abstractmethod
    def get_step(self):
        pass

    # Function returns current denormalized filter impulse response
    @abstractmethod
    def get_impulse(self):
        pass

    # Function returns current filter frequency response (frec,magnitude,phase)
    @abstractmethod
    def get_bode(self):
        pass

    # Function returns current filter group delay
    @abstractmethod
    def get_group_delay(self):
        pass

    # Function returns current filter zeroes and poles (zeroes, poles)
    @abstractmethod
    def get_zeroes_poles(self):
        pass

    # Function returns current filter template limitations
    @abstractmethod
    def get_template(self):
        pass

    # Function returns current filter type:LP,HP,BP or SB
    @abstractmethod
    def filter_is(self):
        pass


