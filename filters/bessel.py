from . import base_filter
import numpy as np
from scipy import signal
from scipy import special

class Bessel(base_filter):

    #Filter initialization with initial parameters received
    def __init__(self, name,Ap,Ao,wpl,wph,wal,wah,gain,n,tao0=None,wrg=None,palm=None):
        if name:
            self.name = name
            self.Ap=Ap
            self.Ao=Ao
            self.wpl=wpl
            self.wph=wph
            self.wal=wal
            self.wah=wah
            self.n=n
            self.gain=gain
            self.tao0=tao0
            self.wrg=wrg
            self.palm=palm
            self.wan=1
            self.poles=[]
            self.zeroes=[]
            self.den=np.poly1d([1])
            self.num=np.poly1d([1])
            self.normalize()#Normalizes current template
            self.do_approximation()#Does normalized approximation and realizes
            self.denormalize()#Denormalizes the approximation to match desired template

    def do_approximation(self):
        self.n=0    #caso que tampoco se va a cumplir
        Tn_Wrgn=1   #caso que nunca se va a cumplir ya que palm va entre 0 y 1
        while  ((1-self.palm) <= Tn_Wrgn) : #condicion para obtener n
            self.n += 1
            self.n=7
            k = np.linspace(0, self.n, num=(self.n + 1))  #k siempre va desde cero hasta n
            k=np.flip(k)
            w, h = signal.freqs(np.poly1d(self.bessel_coef(0)),np.poly1d(self.bessel_coef(k)))
            group_delay = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
            Tn_Wrgn=3
            #wrgn, Tn_Wrgn = signal.group_delay((self.bessel_coef(0), self.bessel_coef(k)),[self.wrgn])
            # me devuelve en Tn_Wrgn el retardo evaluado en Wrgn. Le pase los parametros de la funcion normalizada
        self.norm_sys = signal.TransferFunction(np.poly1d(self.bessel_coef(0)),np.poly1d(self.bessel_coef(k))) #Filter system is obtained
        self.poles=np.roots(np.poly1d(self.bessel_coef(k)))
            


    def bessel_coef (self,k):   #calculo coeficientes bessel
        b_k = (special.factorial(2*self.n-k))/((special.factorial(k))*(special.factorial(self.n-k))*(np.power(2,(self.n-k))))
        #b_k es un vector con todos los coeficientes
        return b_k

    #def bessel_tf(self,k):
    #    num = [self.bessel_coef(0)] #el numerador siempre esta compuesto por b_0
    #    den = [self.bessel_coef(k)]   #el denominador coincide con los coeficientes de bessel
    #    return signal.TransferFunction(num, den)     #funcion transferencia normalizada
