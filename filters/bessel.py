from . import base_filter
import numpy as np
from scipy import signal

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
            self.poles=[]
            self.zeroes=[]
            self.den=np.poly1d([1])
            self.num=np.poly1d([1])
            self.normalize()#Normalizes current template
            self.do_approximation()#Does normalized approximation and realizes
            self.denormalize()#Denormalizes the approximation to match desired template

    def do_approximation(self):
        self.n=0    #caso que tampoco se va a cumplir
        Tn_Wrgn=0   #caso que nunca se va a cumplir ya que palm va entre 0 y 1
        while Tn_Wrgn >= (1-palm): #condicion para obtener n
            self.n += 1
            k = np.linspace(0, self.n, num=(self.n + 1))  #k siempre va desde cero hasta n
            Wrgn, Tn_Wrgn = signal.group_delay((bessel_coef(self,0), bessel_coef(self,k)),[Wrgn])
            # me devuelve en Tn_Wrgn el retardo evaluado en Wrgn. Le pase los parametros de la funcion normalizada
            


    def bessel_coef (self,k):   #calculo coeficientes bessel
        b_k = (math.factorial(2*self.n-k))/((math.factorial(k))*(math.factorial(self.n-k))*(2^(self.n-k)))
        #b_k es un vector con todos los coeficientes
        return b_k

    def bessel_tf(self,k):
        num = [bessel_coef(self,0)] #el numerador siempre esta compuesto por b_0
        den = [bessel_coef(self,k)]   #el denominador coincide con los coeficientes de bessel
        return signal.TransferFunction(num, den)     #funcion transferencia normalizada
