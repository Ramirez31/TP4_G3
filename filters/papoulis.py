from . import base_filter
import numpy as np
from scipy import signal
from scipy import special

class Papoulis(base_filter):

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
       self.epsilon=np.sqrt(np.power(10,self.Ap/10)-1)
       self.n=0
       Ln2_wan=0
       while (Ln2_wan >= ((np.power(10,self.Ao/10)-1)/np.power(self.epsilon,2)))==False:
           self.n+=1
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
                phi=np.polymul(temp_pol,np.poly1d([1,1]))
                Ln=np.polyint(phi)#Primitive is calculated, Barrow will be applied to this function
           
           Ln2=np.poly1d([0])
           poly2mul=np.poly1d([-2,0,-1])#First (-S^2-1) is replaced in the obtained primitive
           for i in range(0,len(Ln)+1):
               polybeingmul=poly2mul#First (-S^2-1) is replaced in the obtained primitive
               for j in range(0,i-1):
                    polybeingmul=np.polymul(polybeingmul,poly2mul)
               a=Ln[i]*polybeingmul
               Ln2=np.polyadd(Ln2,a)
           Ln2=np.polyadd(Ln2,-np.poly1d([np.polyval(Ln,-1)]))
           if np.mod(self.n,2)==0:#If n is even 
               Ln2=Ln2/np.polyval(Ln2,-1j)
           Ln2_wan=np.polyval(Ln2,self.wan)
       roots=np.roots(np.polyadd(Ln2*np.power(self.epsilon,2),np.poly1d([1])))
       for i in range(0,len(roots)):
           if np.real(roots[i])<=0:
                if roots[i]!=0:
                    pol=np.poly1d([-1/roots[i],1])
                else:
                    pol=np.poly1d([1,0])
                self.den=self.den*pol
       self.norm_sys = signal.TransferFunction(self.num,self.den) #Filter system is obtained
       self.poles=np.roots(self.den)
       self.aprox_gain=1