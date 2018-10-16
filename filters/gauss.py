from . import base_filter
import numpy as np
from scipy import signal

class Gauss(base_filter):

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

    
    def factorial(self,n):
        if n <= 0:
            return 1
        else:
            return n*self.factorial(n-1)

    def do_approximation(self):
        #Me pasan Tao(0) WRG y palmerita
        wrgn = self.tao0*self.wrg
        j = np.complex(0,1)
        NoCumplo=True
        self.n=1
        
        while NoCumplo or self.n<=20:
            Terminos= np.zeros(self.n*2)
            Terminos[len(Terminos)]=1
           
            for i in range (0,self.n+1):
                Terminos[len(Terminos)-2*i] = (-1*(j**i))/self.factorial(i)

            POL=np.poly1d(Terminos)
            den=Terminos
            w,TaoWRGN = signal.group_delay((1,den),[wrgn])

            if TaoWRGN >= (1-self.palm):
                NoCumplo=False
            else:
                NoCumplo=True

            self.poles=np.roots(POL)
            self.zeroes=1
            self.n=self.n+1
                