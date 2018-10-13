#!/usr/apps/Python/bin/python
import numpy as np
import matplotlib, sys
matplotlib.use('TkAgg')
from scipy import signal
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends._backend_tk import NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from tkinter import *
import filters

class TP4:
    #Function creates filter according to user input

    def parse_entry(self):
        error=False
        entries=[]
        if self.is_float(self.ap_entry.get()):
            entries.append(float(self.ap_entry.get()))
        else:
            error=True
        if self.is_float(self.aa_entry.get()):
            entries.append(float(self.aa_entry.get()))
        else:
            error=True
        if self.is_float(self.fpl_entry.get()):
            entries.append(float(self.fpl_entry.get()))
        else:
            error=True
        if self.is_float(self.fph_entry.get()):
            entries.append(float(self.fph_entry.get()))
        else:
            error=True
        if self.is_float(self.fal_entry.get()):
            entries.append(float(self.fal_entry.get()))
        else:
            error=True
        if self.is_float(self.fah_entry.get()):
            entries.append(float(self.fah_entry.get()))
        else:
            error=True
        if self.is_float(self.gain_entry.get()):
            entries.append(float(self.gain_entry.get()))
        else:
            error=True
        entries.append(1)
        self.ap_entry.delete(0,END)
        self.aa_entry.delete(0,END)
        self.fal_entry.delete(0,END)
        self.fah_entry.delete(0,END)
        self.fpl_entry.delete(0,END)
        self.fph_entry.delete(0,END)
        self.gain_entry.delete(0,END)
        return entries, error

    def create_filter(self):
        entries,error =self.parse_entry()
        filter_instance = filters.create('butterworth', name='LowPass',Ap=float(entries[0]),Ao=float(entries[1]),wpl=float(entries[2]),wph=float(entries[3]),wal=float(entries[4]),wah=float(entries[5]),gain=float(entries[6]),n=int(entries[7]))
        
        self.w,self.mag,self.phase = filter_instance.get_bode()
        self.Ap,self.Ao,self.wpl,self.wph,self.wal,self.wah = filter_instance.get_template()
        self.filter_type = filter_instance.filter_is()
        
        self.atenua = -(self.mag)
        self.stepT,self.step_mag = filter_instance.get_step()
        self.impT,self.imp_mag = filter_instance.get_impulse()
        self.zeroes, self.poles = filter_instance.get_zeroes_poles()
        self.group_delay = filter_instance.get_group_delay()

    #Function plots current filter's phase in current subplot
    def plot_phase(self):
        self.axis.clear()
        self.axis.semilogx(self.w,self.phase)
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("Radian Frequency [1/rad]$")
        self.axis.set_ylabel("$Phase [Deegres]$")
        self.data_plot.draw()
       
    #Function plots current filter's magnitude in atenuation in current subplot    
    def plot_atenuacion(self):
        if self.filter_type == 'LowPass':
            xl=-100
            yl=self.Ap
            widthl=np.absolute(xl)+self.wpl
            heightl=100
            xr=self.wal
            yr=-100
            widthr=100
            heightr=100+self.Ao
            xc=0
            yc=0
            widthc=0
            heightc=0
        elif self.filter_type == 'HighPass':
            xl=-100
            yl=-100
            widthl=-xl + self.wal
            heightl=np.absolute(yl) + self.Ao
            xr=self.wpl
            yr=self.Ap
            widthr=100
            heightr=100
            xc=0
            yc=0
            widthc=0
            heightc=0
        elif self.filter_type == 'BandPass':
            xl=-100
            yl=-100
            widthl=np.absolute(xl)+self.wal
            heightl=np.absolute(yl)+self.Ao
            xr=self.wah
            yr=-100
            widthr=100
            heightr=np.absolute(yr)+self.Ao
            xc=self.wpl
            yc=self.Ap
            widthc=self.wph-self.wpl
            heightc=100
        elif self.filter_type == 'StopBand':
            xl=-100
            yl=self.Ap
            widthl=np.absolute(xl)+self.wpl
            heightl=100
            xr=self.wph
            yr=self.Ap
            widthr=100
            heightr=100
            xc=self.wal
            yc=-100
            widthc=self.wah-self.wal
            heightc=np.absolute(yc)+self.Ao
        c_rect = matplotlib.patches.Rectangle( (xc,yc), width= widthc, height=heightc, fill=False,color='red')#template is plotted
        l_rect = matplotlib.patches.Rectangle( (xl,yl), width= widthl, height=heightl, fill=False,color='red')#template is plotted
        r_rect = matplotlib.patches.Rectangle( (xr,yr), width= widthr, height=heightr, fill=False,color='red')#template is plotted
        self.axis.clear()
        self.axis.add_patch(l_rect)
        self.axis.add_patch(r_rect)
        self.axis.add_patch(c_rect)
        self.axis.semilogx(self.w,self.atenua)
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("$Radian Frequency [1/rad]$")
        self.axis.set_ylabel("$Attenuation [dB]$")
        self.data_plot.draw()

    def plot_ganancia(self):
        if self.filter_type == 'LowPass':
            xl=-100
            yl=-100
            widthl=np.absolute(xl)+self.wpl
            heightl=np.absolute(yl)-self.Ap
            xr=self.wal
            yr=-self.Ao
            widthr=100
            heightr=100
            xc=0
            yc=0
            widthc=0
            heightc=0
        elif self.filter_type == 'HighPass':
            xl=-100
            yl=-self.Ao
            widthl=nphase.absolute(xl)+self.wal
            heightl=100
            xr=self.wpl
            yr=-100
            widthr=100
            heightr=np.absolute(yr)-self.Ap
            xc=0
            yc=0
            widthc=0
            heightc=0
        elif self.filter_type == 'BandPass':
            xl=-100
            yl=-self.Ao
            widthl=np.absolute(xl)+self.wal
            heightl=100
            xr=self.wah
            yr=-self.Ao
            widthr=100
            heightr=100
            xc=self.wpl
            yc=-100
            widthc=self.wph-self.wpl
            heightc=np.absolute(yc)-self.Ap
        elif self.filter_type == 'StopBand':
            xl=-100
            yl=-100
            widthl=np.absolute(xl)+self.wpl
            heightl=np.absolute(yl)-self.Ap
            xr=self.wph
            yr=-100
            widthr=100
            heightr=np.absolute(yr)-self.Ap
            xc=self.wal
            yc=-self.Ao
            widthc=self.wah-self.wal
            heightc=100
        c_rect = matplotlib.patches.Rectangle( (xc,yc), width= widthc, height=heightc, fill=False,color='red')#template is plotted
        l_rect = matplotlib.patches.Rectangle( (xl,yl), width= widthl, height=heightl, fill=False,color='red')#template is plotted
        r_rect = matplotlib.patches.Rectangle( (xr,yr), width= widthr, height=heightr, fill=False,color='red')#template is plotted
        self.axis.clear()
        self.axis.add_patch(l_rect)
        self.axis.add_patch(r_rect)
        self.axis.add_patch(c_rect)
        self.axis.semilogx(self.w,self.mag)
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("$Radian Frequency [1/rad]$")
        self.axis.set_ylabel("$Gain [dB]$")
        self.data_plot.draw()

    #Function plots current filter's step response in current subplot
    def plot_step(self):
        self.axis.clear()
        self.axis.plot(self.stepT,self.step_mag)
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("$Time [s]$")
        self.axis.set_ylabel("$V_{out} [Volts]$")
        self.data_plot.draw()

    #Function plots current filter's impulse response in current subplot
    def plot_imp(self):
        self.axis.clear()
        self.axis.plot(self.impT,self.imp_mag)
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("$Time [s]$")
        self.axis.set_ylabel("$V_{out} [Volts]$")
        self.data_plot.draw()

    #Function plots current filter's zeroes and poles
    def plot_zeroes_and_poles(self):
        self.axis.clear()
        self.axis.scatter(np.real(self.poles),np.imag(self.poles),marker="x")
        self.axis.scatter(np.real(self.zeroes),np.imag(self.zeroes),marker="o")
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("$Sigma$")
        self.axis.set_ylabel("$jw$")
        self.data_plot.draw()

    def plot_group_delay(self):
        self.axis.clear()
        self.axis.semilogx(self.w,self.group_delay)
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("Radian Frequency [1/rad]$")
        self.axis.set_ylabel("$Group Delay [s]$")
        self.data_plot.draw()

    def is_float(self,value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    def __init__(self):
        self.root = Tk()
        self.root.title("Tc Example")
        #------------------------------------------------------------------------
        side_toolbar=Frame(self.root)
        side_toolbar.pack(side=LEFT,fill=BOTH,expand=True,padx=2,pady=4)

        filter_list = ('LowPass', 'HighPass', 'BandPass','StopBand')
        v = StringVar()
        v.set(filter_list[0])
        dropdown_menu=OptionMenu(side_toolbar,v, *filter_list).grid(row=0,column=0)

        gain_label = Label( side_toolbar, text="Gain:").grid(row=1,column=0)
        self.gain_entry = Entry(side_toolbar,width=5)
        self.gain_entry.grid(row=1,column=1)
        gain_unit = Label( side_toolbar, text="[dB]").grid(row=1,column=2)

        fpl_label = Label( side_toolbar, text="Passband Freq(Fp-):").grid(row=2,column=0)
        self.fpl_entry = Entry(side_toolbar,width=5)
        self.fpl_entry.grid(row=2,column=1)
        fpl_unit = Label( side_toolbar, text="[Hz]").grid(row=2,column=2)

        fph_label = Label( side_toolbar, text="Passband Freq(Fp+):").grid(row=3,column=0)
        self.fph_entry = Entry(side_toolbar,width=5)
        self.fph_entry.grid(row=3,column=1)
        fph_unit = Label( side_toolbar, text="[Hz]").grid(row=3,column=2)

        fal_label = Label( side_toolbar, text="Attenuation Freq(Fa-):").grid(row=4,column=0)
        self.fal_entry = Entry(side_toolbar,width=5)
        self.fal_entry.grid(row=4,column=1)
        fal_unit = Label( side_toolbar, text="[Hz]").grid(row=4,column=2)

        fah_label = Label( side_toolbar, text="Attenuation Freq(Fa+):").grid(row=5,column=0)
        self.fah_entry = Entry(side_toolbar,width=5)
        self.fah_entry.grid(row=5,column=1)
        fah_unit = Label( side_toolbar, text="[Hz]").grid(row=5,column=2)

        ap_label = Label( side_toolbar, text="Attenuation Atten.(Ap):").grid(row=6,column=0)
        self.ap_entry = Entry(side_toolbar,width=5)
        self.ap_entry.grid(row=6,column=1)
        ap_unit = Label( side_toolbar, text="[dB]").grid(row=6,column=2)

        aa_label = Label( side_toolbar, text="Stopband Atten(Aa):").grid(row=7,column=0)
        self.aa_entry = Entry(side_toolbar,width=5)
        self.aa_entry.grid(row=7,column=1)
        aa_unit = Label( side_toolbar, text="[dB]").grid(row=7,column=2)

        button_create_filter = Button(side_toolbar,text="Create Filter",command=self.create_filter).grid(row=8)

        graph_and_buttons = Frame(self.root)
        graph_and_buttons.pack(side=LEFT)
        graph = Canvas(graph_and_buttons)
        graph.pack(side=TOP,fill=BOTH,expand=True,padx=2,pady=4)
        toolbar = Frame(graph_and_buttons)
        button_phase = Button(toolbar,text="Bode Fase",command=self.plot_phase)
        button_phase.pack(side=LEFT,padx=2,pady=2)
        button_mag = Button(toolbar,text="Bode Ganancia",command=self.plot_ganancia)
        button_mag.pack(side=LEFT,padx=2,pady=2)
        button_aten = Button(toolbar,text="Bode Atenuación",command=self.plot_atenuacion)
        button_aten.pack(side=LEFT,padx=2,pady=2)
        button_step = Button(toolbar,text="Escalón",command=self.plot_step)
        button_step.pack(side=LEFT,padx=2,pady=2)
        button_imp = Button(toolbar,text="Impulso",command=self.plot_imp)
        button_imp.pack(side=LEFT,padx=2,pady=4)
        button_zeros_and_poles = Button(toolbar,text="Zeroes and Poles",command=self.plot_zeroes_and_poles)
        button_zeros_and_poles.pack(side=LEFT,padx=2,pady=4)
        button_group_delay = Button(toolbar,text="Group Delay",command=self.plot_group_delay)
        button_group_delay.pack(side=LEFT,padx=2,pady=4)
        toolbar.pack(side=TOP,fill=X)
        
        #-------------------------------------------------------------------------------

        f = Figure()
        self.axis = f.add_subplot(111)

        self.data_plot = FigureCanvasTkAgg(f, master=graph)
        self.data_plot.draw()
        self.data_plot.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        nav = NavigationToolbar2Tk(self.data_plot, graph_and_buttons)
        nav.pack(side=TOP,fill=BOTH,expand=True,padx=2,pady=4)
        nav.update()
        self.data_plot._tkcanvas.pack(side=LEFT, fill=X, expand=True)
        #-------------------------------------------------------------------------------
        self.root.mainloop()

if __name__ == "__main__":
    ex = TP4()