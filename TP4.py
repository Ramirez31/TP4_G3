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
from tkinter import ttk
from tkinter import messagebox
import filters 

class TP4:
    #Function parses user input returning an error if input is not numeric
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
            else:
                error = True
            entry.delete(0,END)
        a=self.nvar.get()
        if (self.nvar.get() == 1) and (self.qvar.get() == 0):
            if self.is_int(self.n_entry.get()):
                entries.append(int(self.n_entry.get()))
            else:
                error = True
        else:
            entries.append(None)
        self.n_entry.delete(0,END)
        if (self.qvar.get() == 1) and (self.nvar.get() == 0):
            if self.is_float(self.q_entry.get()):
                entries.append(float(self.q_entry.get()))
            else:
                error = True
        else:
            entries.append(None)
        self.q_entry.delete(0,END)
        if self.is_float(self.denorm_entry.get()):
            entries.append(float(self.denorm_entry.get()))
        else:
            error = True
        return error,entries,aprox

    #Function creates filter according to user input
    def create_filter(self):
        error,entries,aproximation =self.parse_entry()
        if error is False:
            self.filter_instance = filters.create(aproximation, *entries)
            if self.filter_instance.error_was() == '':
                self.filter_ready=True
                self.w,self.mag,self.phase = self.filter_instance.get_bode()
                self.wn,self.magn,self.phasen=self.filter_instance.get_norm_bode()
                self.template_params = self.filter_instance.get_template()
                self.filter_type = self.filter_instance.filter_is()
                self.n,self.q=self.filter_instance.n_and_q_is()
        
                self.atenua = -(self.mag)
                self.stepT,self.step_mag = self.filter_instance.get_step()
                self.impT,self.imp_mag = self.filter_instance.get_impulse()
                self.zeroes, self.poles = self.filter_instance.get_zeroes_poles()
                self.group_delay = self.filter_instance.get_group_delay()

                for widget in self.filter_data.grid_slaves():
                    widget.grid_forget()

                self.filter_type_label=Label( self.filter_data, text="Filter Type:")
                self.filter_type_label.grid(row=0,column=0)

                self.filter_type_data=Label( self.filter_data, text=self.filter_type)
                self.filter_type_data.grid(row=0,column=2)

                self.n_filter_label=Label( self.filter_data, text="Filter Order:")
                self.n_filter_label.grid(row=1,column=0)

                self.n_filter_data=Label( self.filter_data, text=self.n)
                self.n_filter_data.grid(row=1,column=2)

                self.q_filter_label=Label( self.filter_data, text="Pole's maximum Q:")
                self.q_filter_label.grid(row=2,column=0)

                self.q_filter_data=Label( self.filter_data, text=str(round(self.q,3)))
                self.q_filter_data.grid(row=2,column=2)

                self.filter_data.grid(row=14,columnspan=3,sticky=W)

            else:
                messagebox.showerror("Input Error", self.filter_instance.error_was())
        else:
            messagebox.showerror("Input Error", "Check if any active entry box is empty. Input has to be numeric")

    def create_stages(self):
        root2=Toplevel(self.root)
        myGUI=DesignFilter(root2,self.filter_instance.get_poles(), self.filter_instance.get_zeroes(),self.filter_instance.get_gain())


    #Function plots current filter's phase in current subplot
    def plot_phase(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.semilogx(self.w,self.phase)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("Radian Frequency [1/rad]$")
            self.axis.set_ylabel("$Phase [Deegres]$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")
    
    #Function plots current normalized filter's magnitude in atenuation in current subplot    
    def plot_norm_atten(self):
        if self.filter_ready is True:
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
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function plots current filter's magnitude in atenuation in current subplot    
    def plot_atten(self):
        if self.filter_ready is True:
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
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function creates filter according to user input
    def plot_gain(self):
        if self.filter_ready is True:
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
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function plots current filter's step response in current subplot
    def plot_step(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.plot(self.stepT,self.step_mag)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("$Time [s]$")
            self.axis.set_ylabel("$V_{out} [Volts]$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function plots current filter's impulse response in current subplot
    def plot_imp(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.plot(self.impT,self.imp_mag)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("$Time [s]$")
            self.axis.set_ylabel("$V_{out} [Volts]$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function plots current filter's zeroes and poles
    def plot_zeroes_and_poles(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.scatter(np.real(self.poles),np.imag(self.poles),marker="x")
            self.axis.scatter(np.real(self.zeroes),np.imag(self.zeroes),marker="o")
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("$Sigma$")
            self.axis.set_ylabel("$jw$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")
    
    #Function creates filter according to user input
    def plot_group_delay(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.semilogx(self.w,self.group_delay*1000)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("Radian Frequency [1/rad]$")
            self.axis.set_ylabel("$Group Delay [ms]$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function creates filter according to user input
    def is_float(self,value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    def is_int(self,value):
        try:
            int(value)
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
            self.entry_buttons[0][0].grid(row=2,column=0,sticky=W)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])
            
            self.entry_buttons[1][0].grid(row=3,column=0,sticky=W)
            self.entry_buttons[1][1].grid(row=3,column=1)
            self.entry_buttons[1][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[1][1])

            self.entry_buttons[3][0].grid(row=4,column=0,sticky=W)
            self.entry_buttons[3][1].grid(row=4,column=1)
            self.entry_buttons[3][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[3][1])

            self.entry_buttons[5][0].grid(row=5,column=0,sticky=W)
            self.entry_buttons[5][1].grid(row=5,column=1)
            self.entry_buttons[5][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[5][1])

            self.entry_buttons[6][0].grid(row=6,column=0,sticky=W)
            self.entry_buttons[6][1].grid(row=6,column=1)
            self.entry_buttons[6][2].grid(row=6,column=2)
            self.curr_buttons.append(self.entry_buttons[6][1])

            self.entry_buttons[10][0].grid(row=7,column=0,sticky=W)
            self.entry_buttons[10][2].grid(row=7,column=1)

            self.entry_buttons[11][0].grid(row=8,column=0,sticky=W)
            self.entry_buttons[11][2].grid(row=8,column=1)

            self.entry_buttons[12][0].grid(row=9,column=0,sticky=W)
            self.entry_buttons[12][1].grid(row=9,column=1)
            self.entry_buttons[12][2].grid(row=9,column=2)

            self.button_create_filter.grid(row=10,column=0,sticky=W)

            if self.filter_ready:
                self.filter_data.grid(row=11,columnspan=3,sticky=W)

        elif  self.filter_string.get()=='HighPass':
            self.aprox_figure.create_image(0,0,image=self.photoHP,anchor='nw')
            self.entry_buttons[0][0].grid(row=2,column=0,sticky=W)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])
            
            self.entry_buttons[1][0].grid(row=3,column=0,sticky=W)
            self.entry_buttons[1][1].grid(row=3,column=1)
            self.entry_buttons[1][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[1][1])

            self.entry_buttons[3][0].grid(row=4,column=0,sticky=W)
            self.entry_buttons[3][1].grid(row=4,column=1)
            self.entry_buttons[3][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[3][1])

            self.entry_buttons[5][0].grid(row=5,column=0,sticky=W)
            self.entry_buttons[5][1].grid(row=5,column=1)
            self.entry_buttons[5][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[5][1])

            self.entry_buttons[6][0].grid(row=6,column=0,sticky=W)
            self.entry_buttons[6][1].grid(row=6,column=1)
            self.entry_buttons[6][2].grid(row=6,column=2)
            self.curr_buttons.append(self.entry_buttons[6][1])

            self.entry_buttons[10][0].grid(row=7,column=0,sticky=W)
            self.entry_buttons[10][2].grid(row=7,column=1)

            self.entry_buttons[11][0].grid(row=8,column=0,sticky=W)
            self.entry_buttons[11][2].grid(row=8,column=1)

            self.entry_buttons[12][0].grid(row=9,column=0,sticky=W)
            self.entry_buttons[12][1].grid(row=9,column=1)
            self.entry_buttons[12][2].grid(row=9,column=2)

            self.button_create_filter.grid(row=10,column=0,sticky=W)

            if self.filter_ready:
                self.filter_data.grid(row=11,columnspan=3,sticky=W)

        elif  self.filter_string.get()=='BandPass':
            self.aprox_figure.create_image(0,0,image=self.photoBP,anchor='nw')
            self.entry_buttons[0][0].grid(row=2,column=0,sticky=W)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])
            
            self.entry_buttons[1][0].grid(row=3,column=0,sticky=W)
            self.entry_buttons[1][1].grid(row=3,column=1)
            self.entry_buttons[1][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[1][1])

            self.entry_buttons[2][0].grid(row=4,column=0,sticky=W)
            self.entry_buttons[2][1].grid(row=4,column=1)
            self.entry_buttons[2][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[2][1])

            self.entry_buttons[3][0].grid(row=5,column=0,sticky=W)
            self.entry_buttons[3][1].grid(row=5,column=1)
            self.entry_buttons[3][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[3][1])

            self.entry_buttons[4][0].grid(row=6,column=0,sticky=W)
            self.entry_buttons[4][1].grid(row=6,column=1)
            self.entry_buttons[4][2].grid(row=6,column=2)
            self.curr_buttons.append(self.entry_buttons[4][1])

            self.entry_buttons[5][0].grid(row=7,column=0,sticky=W)
            self.entry_buttons[5][1].grid(row=7,column=1)
            self.entry_buttons[5][2].grid(row=7,column=2)
            self.curr_buttons.append(self.entry_buttons[5][1])

            self.entry_buttons[6][0].grid(row=8,column=0,sticky=W)
            self.entry_buttons[6][1].grid(row=8,column=1)
            self.entry_buttons[6][2].grid(row=8,column=2)
            self.curr_buttons.append(self.entry_buttons[6][1])

            self.entry_buttons[10][0].grid(row=9,column=0,sticky=W)
            self.entry_buttons[10][2].grid(row=9,column=1)

            self.entry_buttons[11][0].grid(row=10,column=0,sticky=W)
            self.entry_buttons[11][2].grid(row=10,column=1)

            self.entry_buttons[12][0].grid(row=11,column=0,sticky=W)
            self.entry_buttons[12][1].grid(row=11,column=1)
            self.entry_buttons[12][2].grid(row=11,column=2)

            self.button_create_filter.grid(row=12,column=0,sticky=W)

            if self.filter_ready:
                self.filter_data.grid(row=13,columnspan=3,sticky=W)

        elif  self.filter_string.get()=='StopBand':
            self.aprox_figure.create_image(0,0,image=self.photoSB,anchor='nw')
            self.entry_buttons[0][0].grid(row=2,column=0,sticky=W)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])

            self.entry_buttons[1][0].grid(row=3,column=0,sticky=W)
            self.entry_buttons[1][1].grid(row=3,column=1)
            self.entry_buttons[1][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[1][1])

            self.entry_buttons[2][0].grid(row=4,column=0,sticky=W)
            self.entry_buttons[2][1].grid(row=4,column=1)
            self.entry_buttons[2][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[2][1])

            self.entry_buttons[3][0].grid(row=5,column=0,sticky=W)
            self.entry_buttons[3][1].grid(row=5,column=1)
            self.entry_buttons[3][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[3][1])

            self.entry_buttons[4][0].grid(row=6,column=0,sticky=W)
            self.entry_buttons[4][1].grid(row=6,column=1)
            self.entry_buttons[4][2].grid(row=6,column=2)
            self.curr_buttons.append(self.entry_buttons[4][1])

            self.entry_buttons[5][0].grid(row=7,column=0,sticky=W)
            self.entry_buttons[5][1].grid(row=7,column=1)
            self.entry_buttons[5][2].grid(row=7,column=2)
            self.curr_buttons.append(self.entry_buttons[5][1])

            self.entry_buttons[6][0].grid(row=8,column=0,sticky=W)
            self.entry_buttons[6][1].grid(row=8,column=1)
            self.entry_buttons[6][2].grid(row=8,column=2)
            self.curr_buttons.append(self.entry_buttons[6][1])

            self.entry_buttons[10][0].grid(row=9,column=0,sticky=W)
            self.entry_buttons[10][2].grid(row=9,column=1)

            self.entry_buttons[11][0].grid(row=10,column=0,sticky=W)
            self.entry_buttons[11][2].grid(row=10,column=1)

            self.entry_buttons[12][0].grid(row=11,column=0,sticky=W)
            self.entry_buttons[12][1].grid(row=11,column=1)
            self.entry_buttons[12][2].grid(row=11,column=2)

            self.button_create_filter.grid(row=12,column=0,sticky=W)

            if self.filter_ready:
                self.filter_data.grid(row=13,columnspan=3,sticky=W)

        elif  self.filter_string.get()=='Group Delay':
            self.entry_buttons[0][0].grid(row=2,column=0,sticky=W)
            self.entry_buttons[0][1].grid(row=2,column=1)
            self.entry_buttons[0][2].grid(row=2,column=2)
            self.curr_buttons.append(self.entry_buttons[0][1])

            self.entry_buttons[7][0].grid(row=3,column=0,sticky=W)
            self.entry_buttons[7][1].grid(row=3,column=1)
            self.entry_buttons[7][2].grid(row=3,column=2)
            self.curr_buttons.append(self.entry_buttons[7][1])

            self.entry_buttons[8][0].grid(row=4,column=0,sticky=W)
            self.entry_buttons[8][1].grid(row=4,column=1)
            self.entry_buttons[8][2].grid(row=4,column=2)
            self.curr_buttons.append(self.entry_buttons[8][1])

            self.entry_buttons[9][0].grid(row=5,column=0,sticky=W)
            self.entry_buttons[9][1].grid(row=5,column=1)
            self.entry_buttons[9][2].grid(row=5,column=2)
            self.curr_buttons.append(self.entry_buttons[9][1])

            self.entry_buttons[10][0].grid(row=6,column=0,sticky=W)
            self.entry_buttons[10][2].grid(row=6,column=1)

            self.entry_buttons[11][0].grid(row=7,column=0,sticky=W)
            self.entry_buttons[11][2].grid(row=7,column=1)

            self.button_create_filter.grid(row=8,column=0,sticky=W)

            if self.filter_ready:
                self.filter_data.grid(row=9,columnspan=3,sticky=W)

    #Function creates filter according to user input
    def __init__(self):
        self.root = Tk()
        self.root.title("Tc Example")
        self.root.resizable(False, False)
        #------------------------------------------------------------------------
        self.side_toolbar=Frame(self.root,width=270)
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
        self.gain_label.grid(row=2,column=0,sticky=W)
        self.gain_entry = Entry(self.side_toolbar,width=5)
        self.gain_entry.grid(row=2,column=1)
        self.gain_unit = Label( self.side_toolbar, text="[dB]")
        self.gain_unit.grid(row=2,column=2)
        self.entry_buttons.append([self.gain_label,self.gain_entry,self.gain_unit])
        self.curr_buttons.append(self.gain_entry)

        self.fpl_label = Label( self.side_toolbar, text="Passband Freq(Fp-):")
        self.fpl_label.grid(row=3,column=0,sticky=W)
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
        self.fal_label.grid(row=5,column=0,sticky=W)
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
        self.ap_label.grid(row=7,column=0,sticky=W)
        self.ap_entry = Entry(self.side_toolbar,width=5)
        self.ap_entry.grid(row=7,column=1)
        self.ap_unit = Label( self.side_toolbar, text="[dB]")
        self.ap_unit.grid(row=7,column=2)
        self.entry_buttons.append([self.ap_label,self.ap_entry,self.ap_unit])
        self.curr_buttons.append(self.ap_entry)

        self.aa_label = Label( self.side_toolbar, text="Stopband Atten(Aa):")
        self.aa_label.grid(row=8,column=0,sticky=W)
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

        self.nvar = IntVar()
        self.n_entry= Entry(self.side_toolbar,width=5)
        self.n_entry.grid(row=9,column=1)
        self.n_check= Checkbutton(self.side_toolbar, text="Fixed  n order ", variable=self.nvar)
        self.n_check.grid(row=9,column=0,sticky=W)
        self.entry_buttons.append([self.n_check,0,self.n_entry])

        self.qvar = IntVar()
        self.q_entry= Entry(self.side_toolbar,width=5)
        self.q_entry.grid(row=10,column=1)
        self.q_check=Checkbutton(self.side_toolbar, text="Maximum Q", variable=self.qvar)
        self.q_check.grid(row=10,column=0,sticky=W)
        self.entry_buttons.append([self.q_check,0,self.q_entry])

        self.denorm_label = Label( self.side_toolbar, text="Denormalization Range:")
        self.denorm_label.grid(row=11,column=0)
        self.denorm_entry = Entry(self.side_toolbar,width=5)
        self.denorm_entry.grid(row=11,column=1)
        self.denorm_unit = Label( self.side_toolbar, text="%")
        self.denorm_unit.grid(row=11,column=2)
        self.entry_buttons.append([self.denorm_label,self.denorm_entry,self.denorm_unit])


        self.button_create_filter = Button(self.side_toolbar,text="Create Filter",command=self.create_filter)
        self.button_create_filter.grid(row=12,column=0,sticky=W)
        #button_create_filter = Button(side_toolbar,text="Create Filter",command=self.create_filter).grid(row=9)
        self.button_create_Stages = Button(self.side_toolbar,text="Create Stages",command=self.create_stages)
        self.button_create_Stages.grid(row=13,column=0,sticky=W)
        
        
        self.filter_data=Frame(self.side_toolbar,width=270)

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
        self.filter_ready=False
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

class DesignFilter:
     def __init__(self,master,poles,zeros,gain):

        self.master=master 
        self.master.title("Design Filters")
       
        self.side_toolbar=Frame(self.master,width=200)
        self.side_toolbar.pack(side=LEFT,fill=BOTH,padx=2,pady=4)
        self.side_toolbar.grid_propagate(0)
        #self.side_toolbar.pack(side=LEFT)

        self.texto=StringVar()
        self.texto.set("Poles & Zeroes Selected:")
       
        self.etiqueta=Label(self.side_toolbar,textvariable=self.texto).place(x=0,y=70)
       
        self.Gain = gain
       #-------Seleccion de Polos--------------------------------
        self.comboPolos = ttk.Combobox(self.side_toolbar,width=10)
        self.comboPolos.place(x=0,y=0)
        self.poles = np.array(np.around(poles, decimals=5)).tolist() #Aca cargo cada polo tmb OJO DEBEN SER ARREGLOS IGUALES EN ORDEN Y TAMAÑO
        self.polesAux = []
        self.polesAux.extend(self.poles)
        self.comboPolos['values'] = self.poles #Aca cargo cada polo
         
        self.comboPolos.current(0)
        self.PolosSeleccionados = [] #Aca guardo polos para hacer etapa

        self.AddPoleButton = Button(self.side_toolbar,text="Add Pole",command=self.AddPole)
        self.AddPoleButton.place(x=0,y=40)
        
        self.SelectedPoles = Listbox(self.side_toolbar,width=10) #Aca muestro polos para hacer estapa
        self.SelectedPoles.place(x=0,y=100)

        
        
        #-----Seleccion de Ceros-----------------------------------
        self.comboZeros = ttk.Combobox(self.side_toolbar,width=10)
        self.comboZeros.place(x=100,y=0)
        
        self.Zeros = np.array(np.around(zeros, decimals=5)).tolist()
        self.ZerosAux = [] #Aca cargo cada polo tmb OJO DEBEN SER ARREGLOS IGUALES EN ORDEN Y TAMAÑOself.ZerosAux.extend(self.Zeros)
        self.comboZeros['values'] = self.Zeros #Aca cargo cada polo
        self.comboZeros.current(0)
        self.ZerosSeleccionados = [] #Aca guardo polos para hacer etapa

        self.AddZeroButton = Button(self.side_toolbar,text="Add Zero",command=self.AddZero)
        self.AddZeroButton.place(x=100,y=40)
       

        self.SelectedZeros = Listbox(self.side_toolbar,width=10) #Aca muestro polos para hacer estapa
        self.SelectedZeros.place(x=100,y=100)


        self.RemoveSelected = Button(self.side_toolbar,text="Del Selection", command=self.Remove)
        self.RemoveSelected.place(x=0,y=270)

        print(self.poles)
        print(self.Zeros)
        print(self.Gain)
       
        #-----------------stages---------------
       
        
        self.forStages=Frame(self.master,width=100)
        self.forStages.pack(side=LEFT,fill=BOTH,padx=2,pady=4)
        self.forStages.grid_propagate(0)
        
        self.texto2=StringVar()
        self.texto2.set("Etapas:")
        self.etiqueta=Label(self.forStages,textvariable=self.texto2).place(x=0,y=2)

        self.SelectedStage = Listbox(self.forStages,width=15,height=25) #Aca muestro polos para hacer estapa
        self.SelectedStage.place(x=0,y=25)

        #-------------Botones---------
        graph_and_buttons = Frame(self.master)
        graph_and_buttons.pack(side=LEFT,fill = BOTH)
        graph = Canvas(graph_and_buttons)
        graph.pack(side=TOP,fill=BOTH,expand=True,padx=2,pady=4)
        toolbar = Frame(graph_and_buttons)
        button_phase = Button(toolbar,text="Bode Phase", command = self.TransferOfStage)
        button_phase.pack(side=LEFT,padx=2,pady=2)
        button_mag = Button(toolbar,text="Bode Module")
        button_mag.pack(side=LEFT,padx=2,pady=2)
        button_step = Button(toolbar,text="Step Response")
        button_step.pack(side=LEFT,padx=2,pady=2)
        button_imp = Button(toolbar,text="Impulse response")
        button_imp.pack(side=LEFT,padx=2,pady=4)
        button_zeros_and_poles = Button(toolbar,text="Zeroes and Poles")
        button_zeros_and_poles.pack(side=LEFT,padx=2,pady=4)
        button_group_delay = Button(toolbar,text="Group Delay")
        button_group_delay.pack(side=LEFT,padx=2,pady=4)
        button_aten = Button(toolbar,text="Individual")
        button_aten.pack(side=LEFT,padx=2,pady=2)
        button_aten_norm = Button(toolbar,text="Accumulative")
        button_aten_norm.pack(side=LEFT,padx=2,pady=2)

        toolbar.pack(side=TOP,fill=X)
        
        f = Figure()
        
        self.axis = f.add_subplot(111)
        self.data_plot = FigureCanvasTkAgg(f, master=graph)
        self.data_plot.draw()
        self.data_plot.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        nav = NavigationToolbar2Tk(self.data_plot, graph_and_buttons)
        nav.pack(side=TOP,fill=BOTH,expand=True,padx=2,pady=4)
        nav.update()
        self.data_plot._tkcanvas.pack(side=LEFT, fill=X, expand=True)

        #----------Generate Stage-----------------------------------------

        self.GenerateStage = Button(self.side_toolbar,text="Generate Stage",command=self.GenerateStage)
        self.GenerateStage.place(x=0,y=300)
        self.TransferList = [] #Voy a ir agregando Stages ej [Polos[] Zeros[], Polos[] Zeros[]] siendo estos poly1d

        #---------Remove Stages--------------------------------------------

        self.DeleteStages = Button(self.side_toolbar,text="Delete Stages",command=self.DeleteStages)
        self.DeleteStages.place(x=0,y=330)

        #--------Automatic Cascade Stages-------------------------------------------

        self.AutoStages = Button(self.side_toolbar,text="Automatic Cascade Stages",command=self.AutoStages)
        self.AutoStages.place(x=0,y=360)
    
     def AddPole(self):

        self.val=self.comboPolos.get()
        
        self.val2=complex(self.val)
        if(self.val2.imag==0):
            self.val2=self.val2.real
        print(self.val2) 
        if(self.val2 not in self.polesAux):
            print("Polo ya Seleccionado")
            return
        if len(self.PolosSeleccionados) == 0 :
            self.PolosSeleccionados.append(self.val2)
            self.SelectedPoles.insert(END,self.val2)
            self.polesAux.remove(self.val2)
            self.comboPolos['values']=self.polesAux
            print(self.polesAux)
            
            if(np.iscomplexobj(self.val2)):
                self.PolosSeleccionados.append(np.conjugate(self.val2))
                self.SelectedPoles.insert(END,np.conjugate(self.val2))
                self.polesAux.remove(np.conjugate(self.val2))
                self.comboPolos['values']=self.polesAux
                print(self.polesAux)
            
        else:
            if (len(self.PolosSeleccionados) <2) and  not np.iscomplexobj(self.val2):
                self.PolosSeleccionados.append(self.val2)
                self.SelectedPoles.insert(END,self.val2)
                print(self.polesAux)
                print(self.val2)
                self.polesAux.remove(self.val2)
                self.comboPolos['values']=self.polesAux
                print(self.polesAux)  

            else:
                
                print("Polo ya Seleccionado o orden maximo de Polos alcanzado para la etapa")
                    
        #Si es imaginario agrego tambien el conjugado pero si ya
        
        
     def AddZero(self):

        self.val=self.comboZeros.get()
        self.val2=complex(self.val)
        if(self.val2.imag==0):
            self.val2=self.val2.real
        if(self.val2 not in self.ZerosAux):
            print("Cero ya Seleccionado")
            return    
        if len(self.ZerosSeleccionados) == 0 :
            self.ZerosSeleccionados.append(self.val2)
            self.SelectedZeros.insert(END,self.val2)
            self.ZerosAux.remove(self.val2)
            self.comboZeros['values']=self.ZerosAux
            print(self.ZerosAux)

            if(np.iscomplexobj(self.val2)):
                self.ZerosSeleccionados.append(np.conjugate(self.val2))
                self.SelectedZeros.insert(END,np.conjugate(self.val2))
                self.ZerosAux.remove(np.conjugate(self.val2))
                self.comboZeros['values']=self.ZerosAux
                print(self.ZerosAux)
        else:
                if (len(self.ZerosSeleccionados) <2) and  not np.iscomplexobj(self.val2):
                    self.ZerosSeleccionados.append(self.val2)
                    self.SelectedZeros.insert(END,self.val2)
                    print(self.ZerosAux)
                    print(self.val2)
                    self.ZerosAux.remove(self.val2)
                    self.comboZeros['values']=self.ZerosAux
                    print(self.ZerosAux)  #Creo que me tira el error en el segundo polo xq hace el for 2 veces con mi nuevo valor alto flash
                else:
                    print("Cero ya Seleccionado o Orden maximo de Ceros alcanzado para la etapa")

     def Remove(self):
         self.ZerosSeleccionados.clear()
         self.PolosSeleccionados.clear()
         self.SelectedPoles.delete(0, END)
         self.SelectedZeros.delete(0,END)
         self.comboPolos['values'] = self.poles
         self.comboZeros['values'] = self.Zeros
         self.polesAux = []
         self.polesAux.extend(self.poles)
         self.ZerosAux = []
         self.ZerosAux.extend(self.Zeros)
         print(self.polesAux)
         
     def GenerateStage(self):
           self.den = []
           self.num = []

           if(len(self.PolosSeleccionados)==0 and len(self.ZerosSeleccionados)==0 ):
               return
           if(len(self.PolosSeleccionados)==0):
                self.den=[]
                self.num.extend(self.ZerosSeleccionados)
           if(len(self.ZerosSeleccionados)==0):
                self.den.extend(self.PolosSeleccionados)
                self.num=[]
           else:                  
                self.den.extend(self.PolosSeleccionados)
                self.num.extend(self.ZerosSeleccionados)
           

           self.HdeStage = []
           self.HdeStage.append(self.den)
           self.HdeStage.append(self.num)
           self.TransferList.append(self.HdeStage)
           print(self.TransferList)
           
           #Elimino los valores de la lsita y queda el combobox sin los valores de los polos seleccionados
           #Hago Clear de los polos y ceros seleccionados para que los proximos no acumulen los anteiores
           self.ZerosSeleccionados.clear()
           self.PolosSeleccionados.clear() 
           self.SelectedPoles.delete(0,END)
           self.SelectedZeros.delete(0,END) 

           self.SelectedStage.insert(END, ["Stage", len(self.TransferList)])
            
         
     def DeleteStages(self):
       
         self.ZerosSeleccionados.clear()
         self.PolosSeleccionados.clear()
         self.SelectedPoles.delete(0,END)
         self.SelectedZeros.delete(0,END)

         self.comboPolos['values'] = self.poles
         self.comboZeros['values'] = self.Zeros
         self.polesAux = []
         self.polesAux.extend(self.poles)
         self.ZerosAux = []
         self.ZerosAux.extend(self.Zeros)

         self.HdeStage = []
         self.TransferList.clear()
         self.SelectedStage.delete(0, END)

     def AutoStages(self):
         #tengo que limpiar el combo box
         self.comboPolos['values'] = ['']
         self.comboZeros['values'] = ['']
         #Tengo que limpiar mi lista
         self.SelectedPoles.delete(0, END)
         self.SelectedZeros.delete(0,END)
         #Tengo que limpiar mis selecciones
         self.ZerosSeleccionados.clear()
         self.PolosSeleccionados.clear()
         #Agarro y genero automaticamente TransferList con mis Aux eliminando los valores de ellos y poniendolos
         self.polesAux = []
         self.polesAux.extend(self.poles)
         self.ZerosAux = []
         self.ZerosAux.extend(self.Zeros)

         print(self.polesAux) 
         for i in self.polesAux:
            if np.iscomplexobj(i):
                self.polesAux.remove(np.conjugate(i))
         
         self.MisPolos = sorted(self.polesAux,key=lambda x:np.sqrt(x.imag**2+x.real**2),reverse=True)       
         

         for i in self.ZerosAux:
            if np.iscomplexobj(i):
                self.ZerosAux.remove(np.conjugate(i))
         
         self.MisCeros = sorted(self.ZerosAux,key=lambda x:np.sqrt(x.imag**2+x.real**2),reverse=True)       
            
         self.HdeStage = []
         self.TransferList.clear() 
         
         print(self.MisPolos)
         print(self.MisCeros)

         while (len(self.MisPolos) != 0) or (len(self.MisCeros) != 0):  
             if len(self.MisPolos)==0 :
                 self.den=[]
             if len(self.MisCeros)==0:
                 self.num=[]   
             for i in self.MisPolos :
                if np.iscomplexobj(i):
                    self.den=[i,np.conjugate(i)]
                    self.MisPolos.remove(i)
                    break
                else:
                    self.den=[i]
                    self.MisPolos.remove(i)
                    break
                break    
             for i in self.MisCeros :
                if np.iscomplexobj(i):
                    self.num=[i,np.conjugate(i)]
                    self.MisCeros.remove(i)
                    break
                else:
                    self.num=[i]
                    self.MisCeros.remove(i)
                    break
                break    
             
             self.HdeStage.append(self.den) 
             self.HdeStage.append(self.num)        
       
             self.TransferList.append(self.HdeStage)
             self.HdeStage = []

         print(self.TransferList)
        #----------------Buttons functions----------------------------
     
     def TransferOfStage(self):
        self.detectStage = self.SelectedStage.curselection()
        for i in self.detectStage:
            self.stageIs = self.TransferList[i]
            #polos de la stage
        self.polesOfStage = self.stageIs[0] #siempre el primer elemento son los polos
        #ceros de la stage
        self.zerosOfStage = self.stageIs[1] #siempre el segundo elemento son los ceros
        #return self.zerosOfStage, self.polesOfStage
        self.den=np.poly1d([1])
        self.num=np.poly1d([1])
        for i in range(0,len(self.polesOfStage)):
            pol = np.poly1d([-1/(self.polesOfStage[i]), 1])
            self.den= self.den*pol
        for i in range(0,len(self.zerosOfStage)):
            pol = np.poly1d([-1/(self.zerosOfStage[i]), 1])
            self.num= self.num*pol
        self.TFofStage = signal.TransferFunction(self.num,self.den)
        print(self.TFofStage)
     

if __name__ == "__main__":
    ex = TP4()