from . import base_filter
import numpy as np
from scipy import signal

class Invchebyshev(base_filter):

    #Filter initialization with initial parameters received
    def __init__(self, name,Ap,Ao,wpl,wph,wal,wah,n):
        if name:
            self.name = name
            self.Ap=Ap
            self.Ao=Ao
            self.wpl=wpl
            self.wph=wph
            self.wal=wal
            self.wah=wah
            self.n=n
            self.poles=[]
            self.zeroes=[]
            self.den=np.poly1d([1])
            self.num=np.poly1d([1])
            self.normalize()#Normalizes current template
            self.do_approximation()#Does normalized approximation and realizes
            self.denormalize()#Denormalizes the approximation to match desired template
    def do_approximation(self):
       pass
               
    def get_step(self):
        return signal.step(self.denorm_sys)

    def get_impulse(self):
        return signal.impulse(self.denorm_sys)

    def get_bode(self):
        self.w,self.mag,self.phase = signal.bode(self.denorm_sys)
        return self.w, self.mag, self.phase

    def get_group_delay(self):
        return -np.gradient(self.phase)

    def get_zeroes_poles(self):
        return self.zeroes, self.poles

    def get_template(self):
        return self.Ap,self.Ao,self.wpl,self.wph,self.wal,self.wah

    def filter_is(self):
        return self.name











