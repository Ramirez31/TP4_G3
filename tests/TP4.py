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
    #Function parses user input returning an error if input is incorrect
    def parse_entry(self):
        error=False
        entries=[]
        aprox=''
        if self.aprox_string.get()=='Butterworth':
            aprox='butterworth'
        elif self.aprox_string.get()=='Chebyshev':
            aprox='chebyshev'
        elif self.aprox_string.get()=='Inverse Chebyshev':
            aprox='invchebyshev'
        elif self.aprox_string.get()=='Legendre':
            aprox='legendre'
        elif self.aprox_string.get()=='Papoulis':
            aprox='papoulis'
        elif self.aprox_string.get()=='Cauer':
            aprox='cauer'
        elif self.aprox_string.get()=='Bessel':
            aprox='bessel'
        elif self.aprox_string.get()=='Gauss':
            aprox='gauss'
        entries.append(self.filter_string.get())
        for entry in self.curr_buttons:
            if self.is_float(entry.get()):
                entries.append(float(entry.get()))
                entry.delete(0,END)
            else:
                error = True
        entries.append(1)#ACA TIENE QUE IR EL BOTONCITO DE SI QUIERO FIJAR UN ORDEN N DE FILTRO

        #error = self.check_template_entry()
        return error,entries,aprox

    #Function checks if user input is correct related to desired template
    def check_template_entry(self):
        error=False
        if float(self.aa_entry.get())>=float(self.ap_entry.get()):
            if  self.filter_string.get()=='LowPass':
                if float(self.fpl_entry.get())>= float(self.fal_entry.get()):
                    error=True
            elif  self.filter_string.get()=='HighPass':
                if float(self.fal_entry.get())>= float(self.fpl_entry.get()):
                    error=True
            elif  self.filter_string.get()=='BandPass':
                if ((float(self.fah_entry.get()) > float(self.fph_entry.get())) and (float(self.fph_entry.get()) > float(self.fpl_entry.get())) and (float(self.fpl_entry.get()) > float(self.fal_entry.get())) ) is False:
                    error=True
            elif  self.filter_string.get()=='StopBand':
                if (((float(self.fph_entry.get()) > float(self.fah_entry.get())) and (float(self.fah_entry.get()) > float(self.fal_entry.get())) and (float(self.fal_entry.get()) > float(self.fpl_entry.get())))) is False:
                    error=True
            elif self.filter_string.get()=='Group Delay':
                pass
        else:
            error=True
        return error

    #Function creates filter according to user input
    def create_filter(self):
        error,entries,aproximation =self.parse_entry()
        if error is False:
            filter_instance = filters.create(aproximation, *entries)
            self.w,self.mag,self.phase = filter_instance.get_bode()
            self.wn,self.magn,self.phasen=filter_instance.get_norm_bode()
            self.template_params = filter_instance.get_template()
            self.filter_type = filter_instance.filter_is()
        
            self.atenua = -(self.mag)
            self.stepT,self.step_mag = filter_instance.get_step()
            self.impT,self.imp_mag = filter_instance.get_impulse()
            self.zeroes, self.poles = filter_instance.get_zeroes_poles()
            self.group_delay = filter_instance.get_group_delay()
        elif True:
            pass

    #Function plots current filter's phase in current subplot
    def plot_phase(self):
        self.axis.clear()
        self.axis.semilogx(self.w,self.phase)
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("Radian Frequency [1/rad]$")
        self.axis.set_ylabel("$Phase [Deegres]$")
        self.data_plot.draw()
    
    #Function plots current normalized filter's magnitude in atenuation in current subplot    
    def plot_norm_atten(self):
        if self.filter_type == 'Group Delay':
            xl=0
            yl=0
            widthl=0
            heightl=0
            xr=0
            yr=0
            widthr=0
            heightr=0
        else:
            xl=-100
            yl=self.template_params[0]
            widthl=np.absolute(xl)+1
            heightl=100
            xr=self.template_params[2]
            yr=-100
            widthr=100
            heightr=100+self.template_params[1]
        l_rect = matplotlib.patches.Rectangle( (xl,yl), width= widthl, height=heightl, fill=False,color='red')#template is plotted
        r_rect = matplotlib.patches.Rectangle( (xr,yr), width= widthr, height=heightr, fill=False,color='red')#template is plotted
        self.axis.clear()
        self.axis.add_patch(l_rect)
        self.axis.add_patch(r_rect)
        self.axis.semilogx(self.wn,-self.magn)
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("$Radian Frequency [1/rad]$")
        self.axis.set_ylabel("$Attenuation [dB]$")
        self.data_plot.draw()    

    #Function plots current filter's magnitude in atenuation in current subplot    
    def plot_atten(self):
        if self.filter_type == 'LowPass':
            xl=-100
            yl=self.template_params[0]-self.template_params[5]
            widthl=np.absolute(xl)+self.template_params[3]
            heightl=100
            xr=self.template_params[4]
            yr=-100
            widthr=100
            heightr=100+self.template_params[1]-self.template_params[5]
            xc=0
            yc=0
            widthc=0
            heightc=0
        elif self.filter_type == 'HighPass':
            xl=-100
            yl=-100
            widthl=-xl + self.template_params[4]
            heightl=np.absolute(yl) + self.template_params[1]-self.template_params[5]
            xr=self.template_params[3]
            yr=self.template_params[0]-self.template_params[5]
            widthr=100
            heightr=100
            xc=0
            yc=0
            widthc=0
            heightc=0
        elif self.filter_type == 'BandPass':
            xl=-100
            yl=-100
            widthl=np.absolute(xl)+self.template_params[5]
            heightl=np.absolute(yl)+self.template_params[1]-self.template_params[7]
            xr=self.template_params[6]
            yr=-100
            widthr=100
            heightr=np.absolute(yr)+self.template_params[1]-self.template_params[7]
            xc=self.template_params[3]
            yc=self.template_params[0]-self.template_params[7]
            widthc=self.template_params[4]-self.template_params[3]
            heightc=100
        elif self.filter_type == 'StopBand':
            xl=-100
            yl=self.template_params[0]-self.template_params[7]
            widthl=np.absolute(xl)+self.template_params[3]
            heightl=100
            xr=self.template_params[4]
            yr=self.template_params[0]-self.template_params[7]
            widthr=100
            heightr=100
            xc=self.template_params[5]
            yc=-100
            widthc=self.template_params[6]-self.template_params[5]
            heightc=np.absolute(yc)+self.template_params[1]-self.template_params[7]
        elif self.filter_type == 'Group Delay':
            xl=0
            yl=0
            widthl=0
            heightl=0
            xr=0
            yr=0
            widthr=0
            heightr=0
            xc=0
            yc=0
            widthc=0
            heightc=0
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

    #Function creates filter according to user input
    def plot_gain(self):
        if self.filter_type == 'LowPass':
            xl=-100
            yl=-100
            widthl=np.absolute(xl)+self.template_params[3]
            heightl=np.absolute(yl)-self.template_params[0]+self.template_params[5]
            xr=self.template_params[4]
            yr=-self.template_params[1]+self.template_params[5]
            widthr=100
            heightr=100
            xc=0
            yc=0
            widthc=0
            heightc=0
        elif self.filter_type == 'HighPass':
            xl=-100
            yl=-self.template_params[1]+self.template_params[5]
            widthl=np.absolute(xl)+self.template_params[4]
            heightl=100
            xr=self.template_params[3]
            yr=-100
            widthr=100
            heightr=np.absolute(yr)-self.template_params[0]+self.template_params[5]
            xc=0
            yc=0
            widthc=0
            heightc=0
        elif self.filter_type == 'BandPass':
            xl=-100
            yl=-self.template_params[1]+self.template_params[7]
            widthl=np.absolute(xl)+self.template_params[5]
            heightl=100
            xr=self.template_params[6]
            yr=-self.template_params[1]+self.template_params[7]
            widthr=100
            heightr=100
            xc=self.template_params[3]
            yc=-100
            widthc=self.template_params[4]-self.template_params[3]
            heightc=np.absolute(yc)-self.template_params[0]+self.template_params[7]
        elif self.filter_type == 'StopBand':
            xl=-100
            yl=-100
            widthl=np.absolute(xl)+self.template_params[3]
            heightl=np.absolute(yl)-self.template_params[0]+self.template_params[7]
            xr=self.template_params[4]
            yr=-100
            widthr=100
            heightr=np.absolute(yr)-self.template_params[0]+self.template_params[7]
            xc=self.template_params[5]
            yc=-self.template_params[1]+self.template_params[7]
            widthc=self.template_params[6]-self.template_params[5]
            heightc=100
        elif self.filter_type == 'Group Delay':
            xl=0
            yl=0
            widthl=0
            heightl=0
            xr=0
            yr=0
            widthr=0
            heightr=0
            xc=0
            yc=0
            widthc=0
            heightc=0
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
    
    #Function creates filter according to user input
    def plot_group_delay(self):
        self.axis.clear()
        self.axis.semilogx(self.w,self.group_delay*1000)
        self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
        self.axis.set_xlabel("Radian Frequency [1/rad]$")
        self.axis.set_ylabel("$Group Delay [ms]$")
        self.data_plot.draw()

    #Function creates filter according to user input
    def is_float(self,value):
        try:
            float(value)
            return True
        except ValueError:
            return False
    
    def set_aproxs(self,*args):
        self.filter_menu['menu'].delete(0, 'end')
        if (self.aprox_string.get()=='Bessel') | (self.aprox_string.get()=='Gauss'):
            filter='Group Delay'
            self.filter_string.set(filter)
            self.filter_menu=OptionMenu(self.side_toolbar,self.filter_string, filter)
            self.filter_menu.grid(row=0,column=0)
        else:
            filter_list = ('LowPass', 'HighPass', 'BandPass','StopBand')
            self.filter_string.set(filter_list[0])
            self.filter_menu=OptionMenu(self.side_toolbar,self.filter_string, *filter_list)
            self.filter_menu.grid(row=0,column=0)

    #Function creates filter according to user input
    def set_entry_buttons(self,*args):
        self.curr_buttons=[]
        self.aprox_figure.delete("all")
        for widget in self.side_toolbar.grid_slaves():
            if int(widget.grid_info()["row"]) > 1:
                widget.grid_forget()
        if  self.filter_string.get()=='LowPass':
            self.aprox_figure.create_image(0,0,image=self.photoLP,anchor='nw')
            self.entry_buttons[0][0].grid(row=2,column=0)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])
            
            self.entry_buttons[1][0].grid(row=3,column=0)
            self.entry_buttons[1][1].grid(row=3,column=1)
            self.entry_buttons[1][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[1][1])

            self.entry_buttons[3][0].grid(row=4,column=0)
            self.entry_buttons[3][1].grid(row=4,column=1)
            self.entry_buttons[3][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[3][1])

            self.entry_buttons[5][0].grid(row=5,column=0)
            self.entry_buttons[5][1].grid(row=5,column=1)
            self.entry_buttons[5][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[5][1])

            self.entry_buttons[6][0].grid(row=6,column=0)
            self.entry_buttons[6][1].grid(row=6,column=1)
            self.entry_buttons[6][2].grid(row=6,column=2)
            self.curr_buttons.append(self.entry_buttons[6][1])
            self.button_create_filter.grid(row=7)
        elif  self.filter_string.get()=='HighPass':
            self.aprox_figure.create_image(0,0,image=self.photoHP,anchor='nw')
            self.entry_buttons[0][0].grid(row=2,column=0)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])
            
            self.entry_buttons[1][0].grid(row=3,column=0)
            self.entry_buttons[1][1].grid(row=3,column=1)
            self.entry_buttons[1][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[1][1])

            self.entry_buttons[3][0].grid(row=4,column=0)
            self.entry_buttons[3][1].grid(row=4,column=1)
            self.entry_buttons[3][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[3][1])

            self.entry_buttons[5][0].grid(row=5,column=0)
            self.entry_buttons[5][1].grid(row=5,column=1)
            self.entry_buttons[5][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[5][1])

            self.entry_buttons[6][0].grid(row=6,column=0)
            self.entry_buttons[6][1].grid(row=6,column=1)
            self.entry_buttons[6][2].grid(row=6,column=2)
            self.curr_buttons.append(self.entry_buttons[6][1])
            self.button_create_filter.grid(row=7)
        elif  self.filter_string.get()=='BandPass':
            self.aprox_figure.create_image(0,0,image=self.photoBP,anchor='nw')
            self.entry_buttons[0][0].grid(row=2,column=0)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])
            
            self.entry_buttons[1][0].grid(row=3,column=0)
            self.entry_buttons[1][1].grid(row=3,column=1)
            self.entry_buttons[1][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[1][1])

            self.entry_buttons[2][0].grid(row=4,column=0)
            self.entry_buttons[2][1].grid(row=4,column=1)
            self.entry_buttons[2][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[2][1])

            self.entry_buttons[3][0].grid(row=5,column=0)
            self.entry_buttons[3][1].grid(row=5,column=1)
            self.entry_buttons[3][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[3][1])

            self.entry_buttons[4][0].grid(row=6,column=0)
            self.entry_buttons[4][1].grid(row=6,column=1)
            self.entry_buttons[4][2].grid(row=6,column=2)
            self.curr_buttons.append(self.entry_buttons[4][1])

            self.entry_buttons[5][0].grid(row=7,column=0)
            self.entry_buttons[5][1].grid(row=7,column=1)
            self.entry_buttons[5][2].grid(row=7,column=2)
            self.curr_buttons.append(self.entry_buttons[5][1])

            self.entry_buttons[6][0].grid(row=8,column=0)
            self.entry_buttons[6][1].grid(row=8,column=1)
            self.entry_buttons[6][2].grid(row=8,column=2)
            self.curr_buttons.append(self.entry_buttons[6][1])
            self.button_create_filter.grid(row=9)
        elif  self.filter_string.get()=='StopBand':
            self.aprox_figure.create_image(0,0,image=self.photoSB,anchor='nw')
            self.entry_buttons[0][0].grid(row=2,column=0)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])

            self.entry_buttons[1][0].grid(row=3,column=0)
            self.entry_buttons[1][1].grid(row=3,column=1)
            self.entry_buttons[1][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[1][1])

            self.entry_buttons[2][0].grid(row=4,column=0)
            self.entry_buttons[2][1].grid(row=4,column=1)
            self.entry_buttons[2][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[2][1])

            self.entry_buttons[3][0].grid(row=5,column=0)
            self.entry_buttons[3][1].grid(row=5,column=1)
            self.entry_buttons[3][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[3][1])

            self.entry_buttons[4][0].grid(row=6,column=0)
            self.entry_buttons[4][1].grid(row=6,column=1)
            self.entry_buttons[4][2].grid(row=6,column=2)
            self.curr_buttons.append(self.entry_buttons[4][1])

            self.entry_buttons[5][0].grid(row=7,column=0)
            self.entry_buttons[5][1].grid(row=7,column=1)
            self.entry_buttons[5][2].grid(row=7,column=2)
            self.curr_buttons.append(self.entry_buttons[5][1])

            self.entry_buttons[6][0].grid(row=8,column=0)
            self.entry_buttons[6][1].grid(row=8,column=1)
            self.entry_buttons[6][2].grid(row=8,column=2)
            self.curr_buttons.append(self.entry_buttons[6][1])
            self.button_create_filter.grid(row=9)
        elif  self.filter_string.get()=='Group Delay':
            self.entry_buttons[0][0].grid(row=2,column=0)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])

            self.entry_buttons[7][0].grid(row=3,column=0)
            self.entry_buttons[7][1].grid(row=3,column=1)
            self.entry_buttons[7][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[7][1])

            self.entry_buttons[8][0].grid(row=4,column=0)
            self.entry_buttons[8][1].grid(row=4,column=1)
            self.entry_buttons[8][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[8][1])

            self.entry_buttons[9][0].grid(row=5,column=0)
            self.entry_buttons[9][1].grid(row=5,column=1)
            self.entry_buttons[9][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[9][1])
            self.button_create_filter.grid(row=6)
    #Function creates filter according to user input
    def __init__(self):
        self.root = Tk()
        self.root.title("Tc Example")
        #------------------------------------------------------------------------
        self.side_toolbar=Frame(self.root,width=300)
        self.side_toolbar.pack(side=LEFT,fill=BOTH,expand=True,padx=2,pady=4)
        self.side_toolbar.grid_propagate(0)

        approximation_list=('Butterworth','Chebyshev','Inverse Chebyshev','Legendre','Papoulis','Cauer','Bessel','Gauss')
        self.aprox_string = StringVar()
        self.aprox_string.set(approximation_list[0])
        aproximation_menu=OptionMenu(self.side_toolbar,self.aprox_string, *approximation_list)
        aproximation_menu.grid(row=0,column=1,columnspan=2)

        filter_list = ('LowPass', 'HighPass', 'BandPass','StopBand')
        self.filter_string = StringVar()
        self.filter_string.set(filter_list[0])
        self.filter_menu=OptionMenu(self.side_toolbar,self.filter_string, *filter_list)
        self.filter_menu.grid(row=0,column=0)

        self.aprox_figure=Canvas(self.side_toolbar,width=222,height=124)
        self.aprox_figure.grid(row=1,columnspan=3)
        self.photoLP=PhotoImage(file="Images\\LowpassImg2.png")
        self.photoHP=PhotoImage(file="Images\\HighpassImg2.png")
        self.photoBP=PhotoImage(file='Images\\BandpassImg2.png')
        self.photoSB=PhotoImage(file='Images\\StopbandImg2.png')

        self.aprox_figure.create_image(0,0,image=self.photoLP,anchor='nw')
        self.entry_buttons=[]
        self.curr_buttons=[]

        self.gain_label = Label( self.side_toolbar, text="Gain:")
        self.gain_label.grid(row=2,column=0)
        self.gain_entry = Entry(self.side_toolbar,width=5)
        self.gain_entry.grid(row=2,column=1)
        self.gain_unit = Label( self.side_toolbar, text="[dB]")
        self.gain_unit.grid(row=2,column=2)
        self.entry_buttons.append([self.gain_label,self.gain_entry,self.gain_unit])
        self.curr_buttons.append(self.gain_entry)

        self.fpl_label = Label( self.side_toolbar, text="Passband Freq(Fp-):")
        self.fpl_label.grid(row=3,column=0)
        self.fpl_entry = Entry(self.side_toolbar,width=5)
        self.fpl_entry.grid(row=3,column=1)
        self.fpl_unit = Label( self.side_toolbar, text="[Hz]")
        self.fpl_unit.grid(row=3,column=2)
        self.entry_buttons.append([self.fpl_label,self.fpl_entry,self.fpl_unit])
        self.curr_buttons.append(self.fpl_entry)

        self.fph_label = Label( self.side_toolbar, text="Passband Freq(Fp+):")
        self.fph_entry = Entry(self.side_toolbar,width=5)
        self.fph_unit = Label( self.side_toolbar, text="[Hz]")
        self.entry_buttons.append([self.fph_label,self.fph_entry,self.fph_unit])

        self.fal_label = Label( self.side_toolbar, text="Attenuation Freq(Fa-):")
        self.fal_label.grid(row=5,column=0)
        self.fal_entry = Entry(self.side_toolbar,width=5)
        self.fal_entry.grid(row=5,column=1)
        self.fal_unit = Label( self.side_toolbar, text="[Hz]")
        self.fal_unit.grid(row=5,column=2)
        self.entry_buttons.append([self.fal_label,self.fal_entry,self.fal_unit])
        self.curr_buttons.append(self.fal_entry)

        self.fah_label = Label( self.side_toolbar, text="Attenuation Freq(Fa+):")
        self.fah_entry = Entry(self.side_toolbar,width=5)
        self.fah_unit = Label( self.side_toolbar, text="[Hz]")
        self.entry_buttons.append([self.fah_label,self.fah_entry,self.fah_unit])

        self.ap_label = Label( self.side_toolbar, text="Attenuation Atten.(Ap):")
        self.ap_label.grid(row=7,column=0)
        self.ap_entry = Entry(self.side_toolbar,width=5)
        self.ap_entry.grid(row=7,column=1)
        self.ap_unit = Label( self.side_toolbar, text="[dB]")
        self.ap_unit.grid(row=7,column=2)
        self.entry_buttons.append([self.ap_label,self.ap_entry,self.ap_unit])
        self.curr_buttons.append(self.ap_entry)

        self.aa_label = Label( self.side_toolbar, text="Stopband Atten(Aa):")
        self.aa_label.grid(row=8,column=0)
        self.aa_entry = Entry(self.side_toolbar,width=5)
        self.aa_entry.grid(row=8,column=1)
        self.aa_unit = Label( self.side_toolbar, text="[dB]")
        self.aa_unit.grid(row=8,column=2)
        self.entry_buttons.append([self.aa_label,self.aa_entry,self.aa_unit])
        self.curr_buttons.append(self.aa_entry)

        self.gp_label = Label( self.side_toolbar, text="Group delay at 0 Hz:")
        self.gp_entry = Entry(self.side_toolbar,width=5)
        self.gp_unit = Label( self.side_toolbar, text="[ms]")
        self.entry_buttons.append([self.gp_label,self.gp_entry,self.gp_unit])

        self.wrg_label = Label( self.side_toolbar, text="(wrg):")
        self.wrg_entry = Entry(self.side_toolbar,width=5)
        self.wrg_unit = Label( self.side_toolbar, text="[rad/s]")
        self.entry_buttons.append([self.wrg_label,self.wrg_entry,self.wrg_unit])

        self.tol_label = Label( self.side_toolbar, text="Group delay Tolerance:")
        self.tol_entry = Entry(self.side_toolbar,width=5)
        self.tol_unit = Label( self.side_toolbar, text="[%]")
        self.entry_buttons.append([self.tol_label,self.tol_entry,self.tol_unit])

        self.button_create_filter = Button(self.side_toolbar,text="Create Filter",command=self.create_filter)
        self.button_create_filter.grid(row=9)

        graph_and_buttons = Frame(self.root)
        graph_and_buttons.pack(side=LEFT)
        graph = Canvas(graph_and_buttons)
        graph.pack(side=TOP,fill=BOTH,expand=True,padx=2,pady=4)
        toolbar = Frame(graph_and_buttons)
        button_phase = Button(toolbar,text="Bode Phase",command=self.plot_phase)
        button_phase.pack(side=LEFT,padx=2,pady=2)
        button_mag = Button(toolbar,text="Bode Denorm Gain",command=self.plot_gain)
        button_mag.pack(side=LEFT,padx=2,pady=2)
        button_aten = Button(toolbar,text="Bode Denorm Atten.",command=self.plot_atten)
        button_aten.pack(side=LEFT,padx=2,pady=2)
        button_aten_norm = Button(toolbar,text="Bode Norm Atten.",command=self.plot_norm_atten)
        button_aten_norm.pack(side=LEFT,padx=2,pady=2)
        button_step = Button(toolbar,text="Step Response",command=self.plot_step)
        button_step.pack(side=LEFT,padx=2,pady=2)
        button_imp = Button(toolbar,text="Impulse response",command=self.plot_imp)
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

        self.aprox_string.trace_add('write',self.set_aproxs)
        self.filter_string.trace_add('write',self.set_entry_buttons)
        #-------------------------------------------------------------------------------
        self.root.mainloop()

if __name__ == "__main__":
    ex = TP4()