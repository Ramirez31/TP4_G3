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
        elif self.aprox_string.get()=='Papoulis':
            aprox='papoulis'
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
        a=self.nvar.get()
        if (self.nvar.get() == 1) and (self.qvar.get() == 0):
            if self.is_int(self.n_entry.get()):
                entries.append(int(self.n_entry.get()))
            else:
                error = True
        else:
            entries.append(None)
        if (self.qvar.get() == 1) and (self.nvar.get() == 0):
            if self.is_float(self.q_entry.get()):
                entries.append(float(self.q_entry.get()))
            else:
                error = True
        else:
            entries.append(None)
        if self.filter_string.get()!='Group Delay':
            if self.denormvar.get()==1:
                if self.is_float(self.denorm_entry.get()):
                    entries.append(float(self.denorm_entry.get()))
                else:
                    error = True
            else:
                entries.append(0)
        return error,entries,aprox

    #Function creates filter according to user input
    def create_filter(self):
        error,entries,aproximation =self.parse_entry()
        if error is False:
            self.filter_instance = filters.create(aproximation, *entries)
            if self.filter_instance.error_was() == '':
                self.filter_ready=True
                self.w,self.mag,self.phase = self.filter_instance.get_bode()
                self.w=self.w/(2*np.pi)
                self.wn,self.magn,self.phasen=self.filter_instance.get_norm_bode()
                self.wn=self.wn/(2*np.pi)
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

                self.aprox_type_label=Label( self.filter_data, text="Aproximation Type:",background='firebrick4', font='Helvetica 9 bold')
                self.aprox_type_label.grid(row=0,column=0)

                self.aprox_type_data=Label( self.filter_data, text=aproximation,background='firebrick4', font='Helvetica 9 bold')
                self.aprox_type_data.grid(row=0,column=2)

                self.filter_type_label=Label( self.filter_data, text="Filter Type:",background='firebrick4', font='Helvetica 9 bold')
                self.filter_type_label.grid(row=1,column=0)

                self.filter_type_data=Label( self.filter_data, text=self.filter_type,background='firebrick4', font='Helvetica 9 bold')
                self.filter_type_data.grid(row=1,column=2)

                self.n_filter_label=Label( self.filter_data, text="Filter Order:",background='firebrick4', font='Helvetica 9 bold')
                self.n_filter_label.grid(row=2,column=0)

                self.n_filter_data=Label( self.filter_data, text=self.n,background='firebrick4', font='Helvetica 9 bold')
                self.n_filter_data.grid(row=2,column=2)

                self.q_filter_label=Label( self.filter_data, text="Pole's maximum Q:",background='firebrick4', font='Helvetica 9 bold')
                self.q_filter_label.grid(row=3,column=0)

                self.q_filter_data=Label( self.filter_data, text=str(round(self.q,3)),background='firebrick4', font='Helvetica 9 bold')
                self.q_filter_data.grid(row=3,column=2)

                self.filter_data.grid(row=14,columnspan=3,sticky=W)
                if self.filter_type == 'Group Delay':
                    self.plot_group_delay()
                else:
                    self.plot_atten()
            else:
                messagebox.showerror("Input Error", self.filter_instance.error_was())
        else:
            messagebox.showerror("Input Error", "Check if any active entry box is empty. Input has to be numeric")

    def create_stages(self):
        if self.filter_ready is True:
            root2=Toplevel(self.root)
            myGUI=DesignFilter(root2,self.filter_instance.get_poles(), self.filter_instance.get_zeroes(),self.filter_instance.get_gain())
        else:
            messagebox.showerror("Error", "No filter was created, stages cannot be created")

    #Function plots current filter's phase in current subplot
    def plot_phase(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.set_title('Bode Diagram phase plot')
            self.axis.semilogx(self.w,self.phase)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("Frequency [Hz]$")
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
                widthl=np.absolute(xl)+1/(2*np.pi)
                heightl=100
                xr=self.template_params[2]
                yr=-100
                widthr=100
                heightr=100+self.template_params[1]
            l_rect = matplotlib.patches.Rectangle( (xl,yl), width= widthl, height=heightl, fill=True,color='orange',alpha=0.5)#template is plotted
            r_rect = matplotlib.patches.Rectangle( (xr,yr), width= widthr, height=heightr, fill=True,color='orange',alpha=0.5)#template is plotted
            self.axis.clear()
            self.axis.set_title('Bode Diagram Normalized Attenuation plot')
            self.axis.add_patch(l_rect)
            self.axis.add_patch(r_rect)
            self.axis.set_xlim(1/10,self.template_params[2]+5)
            self.axis.set_ylim(-5,self.template_params[1]+10)
            self.axis.semilogx(self.wn,-self.magn)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("$Frequency [Hz]$")
            self.axis.set_ylabel("$Attenuation [dB]$")
            self.data_plot.draw()  
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function plots current filter's magnitude in atenuation in current subplot    
    def plot_atten(self):
        if self.filter_ready is True:
            if self.filter_type == 'LowPass':
                xl=-100000
                yl=self.template_params[0]-self.template_params[5]
                widthl=np.absolute(xl)+self.template_params[3]
                heightl=100000
                xr=self.template_params[4]
                yr=-100000
                widthr=100000
                heightr=100000+self.template_params[1]-self.template_params[5]
                xc=0
                yc=0
                widthc=0
                heightc=0
                self.axis.clear()
                self.axis.set_xlim(self.template_params[3]/10,self.template_params[4]*2)
                self.axis.set_ylim(-5-self.template_params[0],self.template_params[1]+self.template_params[5]+5)
            elif self.filter_type == 'HighPass':
                xl=-100000
                yl=-100000
                widthl=-xl + self.template_params[4]
                heightl=np.absolute(yl) + self.template_params[1]-self.template_params[5]
                xr=self.template_params[3]
                yr=self.template_params[0]-self.template_params[5]
                widthr=100000
                heightr=100000
                xc=0
                yc=0
                widthc=0
                heightc=0
                self.axis.clear()
                self.axis.set_xlim(self.template_params[4]/10,self.template_params[3]*2)
                self.axis.set_ylim(-5-self.template_params[0],self.template_params[1]+self.template_params[5]+5)
            elif self.filter_type == 'BandPass':
                xl=-100000
                yl=-100000
                widthl=np.absolute(xl)+self.template_params[5]
                heightl=np.absolute(yl)+self.template_params[1]-self.template_params[7]
                xr=self.template_params[6]
                yr=-100000
                widthr=100000
                heightr=np.absolute(yr)+self.template_params[1]-self.template_params[7]
                xc=self.template_params[3]
                yc=self.template_params[0]-self.template_params[7]
                widthc=self.template_params[4]-self.template_params[3]
                heightc=100000
                self.axis.clear()
                self.axis.set_xlim(self.template_params[5]/5,self.template_params[6]*5)
                self.axis.set_ylim(-5-self.template_params[0],self.template_params[1]+self.template_params[7]+5)
            elif self.filter_type == 'StopBand':
                xl=-100000
                yl=self.template_params[0]-self.template_params[7]
                widthl=np.absolute(xl)+self.template_params[3]
                heightl=100000
                xr=self.template_params[4]
                yr=self.template_params[0]-self.template_params[7]
                widthr=100000
                heightr=100000
                xc=self.template_params[5]
                yc=-100000
                widthc=self.template_params[6]-self.template_params[5]
                heightc=np.absolute(yc)+self.template_params[1]-self.template_params[7]
                self.axis.clear()
                self.axis.set_xlim(self.template_params[3]/5,self.template_params[4]*5)
                self.axis.set_ylim(-5-self.template_params[0],self.template_params[1]+self.template_params[7]+5)
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
                self.axis.clear()
            c_rect = matplotlib.patches.Rectangle( (xc,yc), width= widthc, height=heightc, fill=True,color='orange',alpha=0.5)#template is plotted
            l_rect = matplotlib.patches.Rectangle( (xl,yl), width= widthl, height=heightl, fill=True,color='orange',alpha=0.5)#template is plotted
            r_rect = matplotlib.patches.Rectangle( (xr,yr), width= widthr, height=heightr, fill=True,color='orange',alpha=0.5)#template is plotted
            self.axis.set_title('Bode Diagram Denormalized Attenuation plot')
            self.axis.add_patch(l_rect)
            self.axis.add_patch(r_rect)
            self.axis.add_patch(c_rect)
            self.axis.semilogx(self.w,self.atenua)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("$Frequency [Hz]$")
            self.axis.set_ylabel("$Attenuation [dB]$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function creates filter according to user input
    def plot_gain(self):
        if self.filter_ready is True:
            if self.filter_type == 'LowPass':
                xl=-100000
                yl=-100000
                widthl=np.absolute(xl)+self.template_params[3]
                heightl=np.absolute(yl)-self.template_params[0]+self.template_params[5]
                xr=self.template_params[4]
                yr=-self.template_params[1]+self.template_params[5]
                widthr=self.template_params[4]*1000
                heightr=100000
                xc=0
                yc=0
                widthc=0
                heightc=0
                self.axis.clear()
                self.axis.set_xlim(self.template_params[3]/100,self.template_params[4]*100)
                self.axis.set_ylim(-5-self.template_params[1],self.template_params[0]+self.template_params[5]+5)
            elif self.filter_type == 'HighPass':
                xl=-100000
                yl=-self.template_params[1]+self.template_params[5]
                widthl=np.absolute(xl)+self.template_params[4]
                heightl=100000
                xr=self.template_params[3]
                yr=-100000
                widthr=self.template_params[3]*1000
                heightr=np.absolute(yr)-self.template_params[0]+self.template_params[5]
                xc=0
                yc=0
                widthc=0
                heightc=0
                self.axis.clear()
                self.axis.set_xlim(self.template_params[4]/100,self.template_params[3]*100)
                self.axis.set_ylim(-5-self.template_params[1]+self.template_params[5],self.template_params[0]+self.template_params[5]+5)
            elif self.filter_type == 'BandPass':
                xl=-100000
                yl=-self.template_params[1]+self.template_params[7]
                widthl=np.absolute(xl)+self.template_params[5]
                heightl=100000
                xr=self.template_params[6]
                yr=-self.template_params[1]+self.template_params[7]
                widthr=self.template_params[6]*1000
                heightr=100000
                xc=self.template_params[3]
                yc=-100000
                widthc=self.template_params[4]-self.template_params[3]
                heightc=np.absolute(yc)-self.template_params[0]+self.template_params[7]
                self.axis.clear()
                self.axis.set_xlim(self.template_params[5]/100,self.template_params[6]*100)
                self.axis.set_ylim(-5-self.template_params[1]+self.template_params[7],self.template_params[0]+self.template_params[7]+5)
            elif self.filter_type == 'StopBand':
                xl=-100000
                yl=-100000
                widthl=np.absolute(xl)+self.template_params[3]
                heightl=np.absolute(yl)-self.template_params[0]+self.template_params[7]
                xr=self.template_params[4]
                yr=-100000
                widthr=self.template_params[4]*1000
                heightr=np.absolute(yr)-self.template_params[0]+self.template_params[7]
                xc=self.template_params[5]
                yc=-self.template_params[1]+self.template_params[7]
                widthc=self.template_params[6]-self.template_params[5]
                heightc=100000
                self.axis.clear()
                self.axis.set_xlim(self.template_params[3]/100,self.template_params[4]*100)
                self.axis.set_ylim(-5-self.template_params[1]+self.template_params[7],self.template_params[0]+self.template_params[7]+5)
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
                self.axis.clear()
            self.axis.set_title('Bode Diagram Denormalized Gain plot')
            c_rect = matplotlib.patches.Rectangle( (xc,yc), width= widthc, height=heightc, fill=True,color='orange',alpha=0.5)#template is plotted
            l_rect = matplotlib.patches.Rectangle( (xl,yl), width= widthl, height=heightl, fill=True,color='orange',alpha=0.5)#template is plotted
            r_rect = matplotlib.patches.Rectangle( (xr,yr), width= widthr, height=heightr, fill=True,color='orange',alpha=0.5)#template is plotted
            self.axis.add_patch(l_rect)
            self.axis.add_patch(r_rect)
            self.axis.add_patch(c_rect)
            self.axis.semilogx(self.w,self.mag)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("$Frequency [Hz]$")
            self.axis.set_ylabel("$Gain [dB]$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function plots current filter's step response in current subplot
    def plot_step(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.set_title('Filter\'s Step Response')
            self.axis.set_aspect('auto',adjustable='box')
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
            self.axis.set_title('Filter\'s Impulse Response')
            self.axis.set_aspect('auto',adjustable='box')
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
            self.axis.set_title('Filter\'s Zeroes and Poles')
            self.axis.set_aspect(aspect='equal')
            self.axis.axhline(0, color='black')
            self.axis.axvline(0, color='black')
            maxmod=0
            for pole in self.poles:
                if np.absolute(pole)>maxmod:
                    maxmod=np.absolute(pole)
            for zero in self.zeroes:
                if np.absolute(zero)>maxmod:
                    maxmod=np.absolute(zero)
            self.axis.set_xlim(-maxmod-2,maxmod+2)
            self.axis.set_ylim(-maxmod-2,maxmod+2)
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
            arrows=[]
            self.axis.clear()
            group_delay=self.group_delay
            if self.filter_type =='Group Delay':
                x=-100
                y=0
                width=np.abs(x)+self.template_params[1]
                height=np.abs(y)+1000*self.template_params[0]*(1-self.template_params[2])
            else:
                x=0
                y=0
                width=0
                height=0
            rect = matplotlib.patches.Rectangle( (x,y), width= width, height=height, fill=True,color='orange',alpha=0.5)#template is plotted
            self.axis.add_patch(rect)
            self.axis.set_title('Filter\'s Group Delay')
           #for i in range(0,len(group_delay)-1):
           #    if (np.absolute(group_delay[i]-group_delay[i+1])*1000>40):
           #        arrows.append(matplotlib.patches.Arrow(self.w[i],group_delay[i]*1000,0,-10,width=12))
           #        j=i
           #        while (np.absolute(group_delay[j]-group_delay[i])<0.01) and (j<(len(group_delay)-1)):
           #            group_delay[j]=group_delay[i-1]*0.9
           #            j=j+1
           #        group_delay[j]=group_delay[i-1]
           #        
           #    elif ((group_delay[i]-group_delay[i+1])*1000<-40):
           #        arrows.append(matplotlib.patches.Arrow(self.w[i],group_delay[i]*1000,0,10,width=15))
           #        j=i
           #        while (np.absolute(group_delay[j]-group_delay[i])<0.01) and (j<(len(group_delay)-1)):
           #            group_delay[j]=group_delay[i-1]
           #            j=j+1
           #        group_delay[j]=group_delay[i-1]
           #for arr in arrows:
           #    self.axis.add_patch(arr)
            self.axis.semilogx(self.w,group_delay*1000)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("Frequency [Hz]$")
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
        self.filter_menu.config(highlightbackground='firebrick4')
        self.filter_menu.configure(activebackground = 'firebrick1')
        self.filter_menu.config(bg = 'firebrick1')
        self.filter_menu["menu"].config(bg='firebrick1')

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

            self.button_create_Stages.grid(row=17,column=0,sticky=W)

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

            self.button_create_Stages.grid(row=17,column=0,sticky=W)

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

            self.button_create_Stages.grid(row=17,column=0,sticky=W)

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

            self.button_create_Stages.grid(row=17,column=0,sticky=W)

        elif  self.filter_string.get()=='Group Delay':
            self.aprox_figure.create_image(0,0,image=self.photoLP,anchor='nw')
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

            self.button_create_Stages.grid(row=17,column=0,sticky=W)

    #Function creates filter according to user input
    def __init__(self):
        self.root = Tk()
        self.root.configure(background='firebrick4')
        self.root.title("Tc Example")
        self.root.resizable(False, False)
        #------------------------------------------------------------------------
        self.side_toolbar=Frame(self.root,width=330,borderwidth=7,relief=RAISED,background='firebrick3')
        self.side_toolbar.pack(side=LEFT,fill=BOTH,expand=True,padx=2,pady=4)
        self.side_toolbar.grid_propagate(0)

        approximation_list=('Butterworth','Chebyshev','Inverse Chebyshev','Papoulis','Bessel','Gauss')
        self.aprox_string = StringVar()
        self.aprox_string.set(approximation_list[0])
        aproximation_menu=OptionMenu(self.side_toolbar,self.aprox_string, *approximation_list)
        aproximation_menu.grid(row=0,column=1,columnspan=2,padx=5,pady=5)
        aproximation_menu.config(highlightbackground='firebrick4')
        aproximation_menu.configure(activebackground = 'firebrick1')
        aproximation_menu.config(bg = 'firebrick1')
        aproximation_menu["menu"].config(bg='firebrick1')

        filter_list = ('LowPass', 'HighPass', 'BandPass','StopBand')
        self.filter_string = StringVar()
        self.filter_string.set(filter_list[0])
        self.filter_menu=OptionMenu(self.side_toolbar,self.filter_string, *filter_list)
        self.filter_menu.grid(row=0,column=0,padx=5,pady=5)
        self.filter_menu.config(highlightbackground='firebrick4')
        self.filter_menu.configure(activebackground = 'firebrick1')
        self.filter_menu.config(bg = 'firebrick1')
        self.filter_menu["menu"].config(bg='firebrick1')

        self.aprox_figure=Canvas(self.side_toolbar,width=222,height=124)
        self.aprox_figure.grid(row=1,columnspan=3,padx=5,pady=5)
        self.aprox_figure.config(highlightbackground='firebrick4')
        self.photoLP=PhotoImage(file="Images\\LowpassImg2.png")
        self.photoHP=PhotoImage(file="Images\\HighpassImg2.png")
        self.photoBP=PhotoImage(file='Images\\BandpassImg2.png')
        self.photoSB=PhotoImage(file='Images\\StopbandImg2.png')

        self.aprox_figure.create_image(0,0,image=self.photoLP,anchor='nw')
        self.entry_buttons=[]
        self.curr_buttons=[]

        self.gain_label = Label( self.side_toolbar, text="Gain:",background='firebrick3', font='Helvetica 9 bold')
        self.gain_label.grid(row=2,column=0,sticky=W)
        self.gain_entry = Entry(self.side_toolbar,width=10)
        self.gain_entry.grid(row=2,column=1)
        self.gain_unit = Label( self.side_toolbar, text="[dB]",background='firebrick3', font='Helvetica 9 bold')
        self.gain_unit.grid(row=2,column=2)
        self.entry_buttons.append([self.gain_label,self.gain_entry,self.gain_unit])
        self.curr_buttons.append(self.gain_entry)

        self.fpl_label = Label( self.side_toolbar, text="Passband Freq(Fp-):",background='firebrick3', font='Helvetica 9 bold')
        self.fpl_label.grid(row=3,column=0,sticky=W)
        self.fpl_entry = Entry(self.side_toolbar,width=10)
        self.fpl_entry.grid(row=3,column=1)
        self.fpl_unit = Label( self.side_toolbar, text="[Hz]",background='firebrick3', font='Helvetica 9 bold')
        self.fpl_unit.grid(row=3,column=2)
        self.entry_buttons.append([self.fpl_label,self.fpl_entry,self.fpl_unit])
        self.curr_buttons.append(self.fpl_entry)

        self.fph_label = Label( self.side_toolbar, text="Passband Freq(Fp+):",background='firebrick3', font='Helvetica 9 bold')
        self.fph_entry = Entry(self.side_toolbar,width=10)
        self.fph_unit = Label( self.side_toolbar, text="[Hz]",background='firebrick3', font='Helvetica 9 bold')
        self.entry_buttons.append([self.fph_label,self.fph_entry,self.fph_unit])

        self.fal_label = Label( self.side_toolbar, text="Attenuation Freq(Fa-):",background='firebrick3', font='Helvetica 9 bold')
        self.fal_label.grid(row=5,column=0,sticky=W)
        self.fal_entry = Entry(self.side_toolbar,width=10)
        self.fal_entry.grid(row=5,column=1)
        self.fal_unit = Label( self.side_toolbar, text="[Hz]",background='firebrick3', font='Helvetica 9 bold')
        self.fal_unit.grid(row=5,column=2)
        self.entry_buttons.append([self.fal_label,self.fal_entry,self.fal_unit])
        self.curr_buttons.append(self.fal_entry)

        self.fah_label = Label( self.side_toolbar, text="Attenuation Freq(Fa+):",background='firebrick3', font='Helvetica 9 bold')
        self.fah_entry = Entry(self.side_toolbar,width=10)
        self.fah_unit = Label( self.side_toolbar, text="[Hz]",background='firebrick3', font='Helvetica 9 bold')
        self.entry_buttons.append([self.fah_label,self.fah_entry,self.fah_unit])

        self.ap_label = Label( self.side_toolbar, text="Attenuation Atten.(Ap):",background='firebrick3', font='Helvetica 9 bold')
        self.ap_label.grid(row=7,column=0,sticky=W)
        self.ap_entry = Entry(self.side_toolbar,width=10)
        self.ap_entry.grid(row=7,column=1)
        self.ap_unit = Label( self.side_toolbar, text="[dB]",background='firebrick3', font='Helvetica 9 bold')
        self.ap_unit.grid(row=7,column=2)
        self.entry_buttons.append([self.ap_label,self.ap_entry,self.ap_unit])
        self.curr_buttons.append(self.ap_entry)

        self.aa_label = Label( self.side_toolbar, text="Stopband Atten(Aa):",background='firebrick3', font='Helvetica 9 bold')
        self.aa_label.grid(row=8,column=0,sticky=W)
        self.aa_entry = Entry(self.side_toolbar,width=10)
        self.aa_entry.grid(row=8,column=1)
        self.aa_unit = Label( self.side_toolbar, text="[dB]",background='firebrick3', font='Helvetica 9 bold')
        self.aa_unit.grid(row=8,column=2)
        self.entry_buttons.append([self.aa_label,self.aa_entry,self.aa_unit])
        self.curr_buttons.append(self.aa_entry)

        self.gp_label = Label( self.side_toolbar, text="Group delay at 0 Hz:",background='firebrick3', font='Helvetica 9 bold')
        self.gp_entry = Entry(self.side_toolbar,width=10)
        self.gp_unit = Label( self.side_toolbar, text="[ms]",background='firebrick3', font='Helvetica 9 bold')
        self.entry_buttons.append([self.gp_label,self.gp_entry,self.gp_unit])

        self.wrg_label = Label( self.side_toolbar, text="Maximum Tol Freq:",background='firebrick3', font='Helvetica 9 bold')
        self.wrg_entry = Entry(self.side_toolbar,width=10)
        self.wrg_unit = Label( self.side_toolbar, text="[Hz]",background='firebrick3', font='Helvetica 9 bold')
        self.entry_buttons.append([self.wrg_label,self.wrg_entry,self.wrg_unit])

        self.tol_label = Label( self.side_toolbar, text="Group delay Tolerance:",background='firebrick3', font='Helvetica 9 bold')
        self.tol_entry = Entry(self.side_toolbar,width=10)
        self.tol_unit = Label( self.side_toolbar, text="[%]",background='firebrick3', font='Helvetica 9 bold')
        self.entry_buttons.append([self.tol_label,self.tol_entry,self.tol_unit])

        self.nvar = IntVar()
        self.n_entry= Entry(self.side_toolbar,width=10)
        self.n_entry.grid(row=9,column=1)
        self.n_check= Checkbutton(self.side_toolbar, text="Fixed  n order ", variable=self.nvar,background='firebrick3', font='Helvetica 9 bold')
        self.n_check.configure(activebackground = 'firebrick3')
        self.n_check.grid(row=9,column=0,sticky=W)
        self.entry_buttons.append([self.n_check,0,self.n_entry])

        self.qvar = IntVar()
        self.q_entry= Entry(self.side_toolbar,width=10)
        self.q_entry.grid(row=10,column=1)
        self.q_check=Checkbutton(self.side_toolbar, text="Maximum Q", variable=self.qvar,background='firebrick3', font='Helvetica 9 bold')
        self.q_check.configure(activebackground = 'firebrick3')
        self.q_check.grid(row=10,column=0,sticky=W)
        self.entry_buttons.append([self.q_check,0,self.q_entry])

        self.denormvar=IntVar()
        self.denorm_check = Checkbutton( self.side_toolbar, text="Denormalization Range:", variable=self.denormvar,background='firebrick3', font='Helvetica 9 bold')
        self.denorm_check.grid(row=11,column=0)
        self.denorm_check.configure(activebackground = 'firebrick3')
        self.denorm_entry = Entry(self.side_toolbar,width=10)
        self.denorm_entry.grid(row=11,column=1)
        self.denorm_unit = Label( self.side_toolbar, text="[%]",background='firebrick3', font='Helvetica 9 bold')
        self.denorm_unit.grid(row=11,column=2)
        self.entry_buttons.append([self.denorm_check,self.denorm_entry,self.denorm_unit])


        self.button_create_filter = Button(self.side_toolbar,text="Create Filter",command=self.create_filter)
        self.button_create_filter.grid(row=12,column=0,sticky=W,padx=5,pady=5)
        self.button_create_filter.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')
        #button_create_filter = Button(side_toolbar,text="Create Filter",command=self.create_filter).grid(row=9)
        self.button_create_Stages = Button(self.side_toolbar,text="Create Stages",command=self.create_stages)
        self.button_create_Stages.grid(row=17,column=0,sticky=W,padx=5,pady=5)
        self.button_create_Stages.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')

        
        self.filter_data=Frame(self.side_toolbar,width=270,borderwidth=2,relief=RAISED,background='firebrick4')

        graph_and_buttons = Frame(self.root,borderwidth=2,relief=RAISED,background='firebrick3')
        graph_and_buttons.pack(side=LEFT)
        graph = Canvas(graph_and_buttons)
        graph.pack(side=TOP,fill=BOTH,expand=True,padx=2,pady=4)
        toolbar = Frame(graph_and_buttons,borderwidth=2,relief=RAISED,background='firebrick3')
        button_phase = Button(toolbar,text="Bode Phase",command=self.plot_phase)
        button_phase.pack(side=LEFT,padx=2,pady=2)
        button_phase.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')
        button_mag = Button(toolbar,text="Bode Denorm Gain",command=self.plot_gain)
        button_mag.pack(side=LEFT,padx=2,pady=2)
        button_mag.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')
        button_aten = Button(toolbar,text="Bode Denorm Atten.",command=self.plot_atten)
        button_aten.pack(side=LEFT,padx=2,pady=2)
        button_aten.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')
        button_norm_aten = Button(toolbar,text="Bode Norm Atten.",command=self.plot_norm_atten)
        button_norm_aten.pack(side=LEFT,padx=2,pady=2)
        button_norm_aten.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')
        button_step = Button(toolbar,text="Step Response",command=self.plot_step)
        button_step.pack(side=LEFT,padx=2,pady=2)
        button_step.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')
        button_imp = Button(toolbar,text="Impulse response",command=self.plot_imp)
        button_imp.pack(side=LEFT,padx=2,pady=4)
        button_imp.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')
        button_zeros_and_poles = Button(toolbar,text="Zeroes and Poles",command=self.plot_zeroes_and_poles)
        button_zeros_and_poles.pack(side=LEFT,padx=2,pady=4)
        button_zeros_and_poles.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')
        button_group_delay = Button(toolbar,text="Group Delay",command=self.plot_group_delay)
        button_group_delay.pack(side=LEFT,padx=2,pady=4)
        button_group_delay.configure(highlightbackground='firebrick4',activebackground = 'brown1',bg = 'firebrick1')
        toolbar.pack(side=TOP,fill=X)
        
        #-------------------------------------------------------------------------------
        self.f = Figure()
        self.filter_ready=False
        self.axis = self.f.add_subplot(111)
        self.data_plot = FigureCanvasTkAgg(self.f, master=graph)
        self.data_plot.draw()
        self.data_plot.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        nav = NavigationToolbar2Tk(self.data_plot, graph_and_buttons)
        nav.pack(side=TOP,fill=BOTH,expand=True,padx=2,pady=4)
        nav.update()
        nav.configure(background='firebrick4')
        nav._message_label.config(background='firebrick4')
        self.data_plot._tkcanvas.pack(side=LEFT, fill=X, expand=True)


        self.aprox_string.trace_add('write',self.set_aproxs)
        self.filter_string.trace_add('write',self.set_entry_buttons)
        #-------------------------------------------------------------------------------
        self.root.mainloop()

class DesignFilter:
     def __init__(self,master,poles,zeros,gain):

        self.master=master 
        self.master.configure(background='skyblue4')
        self.master.title("Design Filters")
        self.master.resizable(False, False)
       
        self.side_toolbar=Frame(self.master,width=300,borderwidth=7,relief=RAISED,background='skyblue3')
        self.side_toolbar.pack(side=LEFT,fill=BOTH,padx=2,pady=4)
        #self.side_toolbar.grid_propagate(0)
        #self.side_toolbar.pack(side=LEFT)

        self.texto=StringVar()
        self.texto.set("Poles & Zeroes Selected: ")
       
        self.etiqueta=Label(self.side_toolbar,textvariable=self.texto,background='skyblue3',font='Helvetica 9 bold')
        self.etiqueta.grid(row=2,column=0,columnspan=2)
        self.Gain = gain
       #-------Seleccion de Polos--------------------------------
        self.comboPolos = ttk.Combobox(self.side_toolbar)
        self.comboPolos.grid(row=0,column=0)

        self.poles = np.array(np.around(poles, decimals=5)).tolist() #Aca cargo cada polo tmb OJO DEBEN SER ARREGLOS IGUALES EN ORDEN Y TAMAÑO
        self.polesAux = []
        self.polesAux.extend(self.poles)
        self.comboPolos['values'] = self.poles #Aca cargo cada polo
         
        self.PolosSeleccionados = [] #Aca guardo polos para hacer etapa

        self.AddPoleButton = Button(self.side_toolbar,text="Add Pole",command=self.AddPole)
        self.AddPoleButton.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        self.AddPoleButton.grid(row=1,column=0)
        
        self.SelectedPoles = Listbox(self.side_toolbar,height=15) #Aca muestro polos para hacer estapa
        self.SelectedPoles.grid(row=3,column=0)
  
        #-----Seleccion de Ceros-----------------------------------
        self.comboZeros = ttk.Combobox(self.side_toolbar)
        self.comboZeros.grid(row=0,column=1)
        
        self.Zeros = np.array(np.around(zeros, decimals=5)).tolist()
        self.ZerosAux = [] #Aca cargo cada polo tmb OJO DEBEN SER ARREGLOS IGUALES EN ORDEN Y TAMAÑOself.ZerosAux.extend(self.Zeros)
        self.ZerosAux.extend(self.Zeros)
        self.comboZeros['values'] = self.Zeros #Aca cargo cada polo
        self.ZerosSeleccionados = [] #Aca guardo polos para hacer etapa

        self.AddZeroButton = Button(self.side_toolbar,text="Add Zero",command=self.AddZero)
        self.AddZeroButton.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        self.AddZeroButton.grid(row=1,column=1)

        self.SelectedZeros = Listbox(self.side_toolbar,height=15) #Aca muestro polos para hacer estapa
        self.SelectedZeros.grid(row=3,column=1)

        self.RemoveSelected = Button(self.side_toolbar,text="Del Selection", command=self.Remove)
        self.RemoveSelected.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        self.RemoveSelected.grid(row=5,column=0, sticky=W, padx=4)

        print(self.poles)
        print(self.Zeros)
        print(self.Gain)
       
        #-----------------stages---------------
       
        self.forStages=Frame(self.master,width=100,borderwidth=5,relief=RAISED,background='skyblue3')
        self.forStages.pack(side=LEFT,fill=BOTH,padx=2,pady=4)
        
        self.texto2=StringVar()
        self.texto2.set("Stages:")
        self.etiqueta=Label(self.forStages,textvariable=self.texto2,background='skyblue3',font='Helvetica 9 bold')
        self.etiqueta.grid(row=0,column=0,columnspan=3)

        self.SelectedStage = Listbox(self.forStages,width=22,height=25) #Aca muestro polos para hacer estapa
        self.SelectedStage.grid(row=1,column=0,columnspan=3)
        self.SelectedStage.bind('<<ListboxSelect>>', self.onselect)

        self.gain_label = Label( self.forStages, text="Gain:",background='skyblue3',font='Helvetica 9 bold')
        self.gain_label.grid(row=2,column=0,sticky=W)
        self.gain_entry = Entry(self.forStages,width=5)
        self.gain_entry.grid(row=2,column=1)
        self.gain_unit = Label( self.forStages, text="[dB]",background='skyblue3',font='Helvetica 9 bold')
        self.gain_unit.grid(row=2,column=2)
        

        #-------------Botones---------
        graph_and_buttons = Frame(self.master,borderwidth=7,relief=RAISED,background='skyblue3')
        graph_and_buttons.pack(side=LEFT,fill = BOTH)
        graph = Canvas(graph_and_buttons)
        graph.pack(side=TOP,fill=BOTH,expand=True,padx=2,pady=4)
        toolbar = Frame(graph_and_buttons,borderwidth=2,relief=RAISED,background='skyblue3')
        button_phase = Button(toolbar,text="Bode Phase", command = self.plot_phase)
        button_phase.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        button_phase.pack(side=LEFT,padx=2,pady=2)
        button_mag = Button(toolbar,text="Bode Module", command = self.plot_atten)
        button_mag.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        button_mag.pack(side=LEFT,padx=2,pady=2)
        button_step = Button(toolbar,text="Step Response", command = self.plot_step)
        button_step.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        button_step.pack(side=LEFT,padx=2,pady=2)
        button_imp = Button(toolbar,text="Impulse response",command = self.plot_imp)
        button_imp.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        button_imp.pack(side=LEFT,padx=2,pady=4)
        button_zeros_and_poles = Button(toolbar,text="Zeroes and Poles", command = self.plot_zeroes_and_poles)
        button_zeros_and_poles.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        button_zeros_and_poles.pack(side=LEFT,padx=2,pady=4)
        button_group_delay = Button(toolbar,text="Group Delay", command = self.plot_group_delay)
        button_group_delay.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        button_group_delay.pack(side=LEFT,padx=2,pady=4)
        button_accumul = Button(toolbar,text="Accumulative",command = self.plot_accumul)
        button_accumul.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        button_accumul.pack(side=LEFT,padx=2,pady=2)

        toolbar.pack(side=TOP,fill=X)
        
        f = Figure()
        
        self.axis = f.add_subplot(111)
        self.data_plot = FigureCanvasTkAgg(f, master=graph)
        self.data_plot.draw()
        self.data_plot.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        nav = NavigationToolbar2Tk(self.data_plot, graph_and_buttons)
        nav.pack(side=TOP,fill=BOTH,expand=True,padx=2,pady=4)
        nav.configure(background='skyblue4')
        nav._message_label.config(background='skyblue4')
        nav.update()
        self.data_plot._tkcanvas.pack(side=LEFT, fill=X, expand=True)
        self.filter_ready=False

        #----------Generate Stage-----------------------------------------

        self.GenerateStage = Button(self.side_toolbar,text="Generate Stage",command=self.GenerateStage)
        self.GenerateStage.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        self.GenerateStage.grid(row=4,column=1, sticky=W, padx=4, pady=2)
        self.TransferList = [] #Voy a ir agregando Stages ej [Polos[] Zeros[], Polos[] Zeros[]] siendo estos poly1d
        self.GainOfStages = [] #guardo las ganancias de cada etapa. que me ingresa el usuario
        self.GainOfStagesVeces = []
        self.ListaDePolosPasados = []
        self.ListaDeCerosPasados = []
        #---------Remove Stages--------------------------------------------

        self.DeleteStages = Button(self.side_toolbar,text="Delete Stages",command=self.DeleteStages)
        self.DeleteStages.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        self.DeleteStages.grid(row=4,column=0, sticky=W, padx=4, pady=2)

        #--------Automatic Cascade Stages-------------------------------------------

        self.AutoStages = Button(self.side_toolbar,text="Automatic Cascade Stages",command=self.AutoStages)
        self.AutoStages.configure(highlightbackground='skyblue3',activebackground = 'lightskyblue1',bg = 'lightskyblue2')
        self.AutoStages.grid(row=5,column=1, sticky=W, padx=4)
  
        #-------Calculo de Rango Dinamico-------------------------------------------
        self.RDT = StringVar()
        self.RDT.set("RD:")
        self.RD_label = Label(self.forStages, textvariable=self.RDT,background='skyblue3',font='Helvetica 9 bold')
        self.RD_label.grid(row=3,column=0,sticky=W)
        
        self.RD_unit = Label(self.forStages, text="[dB]",background='skyblue3',font='Helvetica 9 bold')
        self.RD_unit.grid(row=3,column=1)

        self.RD_label2 = Label(self.forStages, text="Vsat = 14V",background='skyblue3',font='Helvetica 6')
        self.RD_label2.grid(row=4,column=0,sticky=W)
        self.RD_entry = Entry(self.forStages,width=5)

        self.RD_labe3l = Label(self.forStages, text="Vmin = 10mV",background='skyblue3',font='Helvetica 6')
        self.RD_labe3l.grid(row=5,column=0,sticky=W)
        self.RD_entry = Entry(self.forStages,width=5)

    #Function plots current filter's phase in current subplot
     def plot_phase(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.set_title('Stage\'s Phase plot')
            self.axis.semilogx(self.w,self.phase)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("Frequency [Hz]$")
            self.axis.set_ylabel("$Phase [Deegres]$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")
    
    #Function plots current filter's magnitude in atenuation in current subplot    
     def plot_atten(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.set_title('Stage\'s Attenuation plot')
            self.axis.semilogx(self.w,-self.mag)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("$Frequency [Hz]$")
            self.axis.set_ylabel("$Attenuation [dB]$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")

    #Function plots current filter's step response in current subplot
     def plot_step(self):
        if self.filter_ready is True:
            self.axis.clear()
            self.axis.set_title('Stage\'s Step Response')
            self.axis.set_aspect(aspect='auto')
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
            self.axis.set_title('Stage\'s Impulse Response')
            self.axis.clear()
            self.axis.set_aspect(aspect='auto')
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
            self.axis.set_title('Stage\'s Zeroes and Poles')
            self.axis.set_aspect(aspect='equal')
            self.axis.axhline(0, color='black')
            self.axis.axvline(0, color='black')
            maxmod=0
            for pole in self.polesOfStage:
                if np.absolute(pole)>maxmod:
                    maxmod=np.absolute(pole)
            for zero in self.zerosOfStage:
                if np.absolute(zero)>maxmod:
                    maxmod=np.absolute(zero)
            self.axis.set_xlim(-maxmod-2,maxmod+2)
            self.axis.set_ylim(-maxmod-2,maxmod+2)
            self.axis.scatter(np.real(self.polesOfStage),np.imag(self.polesOfStage),marker="x")
            self.axis.scatter(np.real(self.zerosOfStage),np.imag(self.zerosOfStage),marker="o")
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
            self.axis.set_title('Stage\'s Group Delay')
            self.axis.semilogx(self.w,self.group_delay*1000)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("Frequency [Hz]$")
            self.axis.set_ylabel("$Group Delay [ms]$")
            self.data_plot.draw()
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")        

     def plot_accumul(self):
        if self.filter_ready is True:
            self.detectStage = self.SelectedStage.curselection()
            den=np.poly1d([1])
            num=np.poly1d([1])
            for i in range(0,self.detectStage[0]+1):
                self.stageIs = self.TransferList[i]
                temp_poles = self.stageIs[0] #siempre el primer elemento son los polos
                temp_zeros = self.stageIs[1] #siempre el segundo elemento son los ceros
                for pole in temp_poles:
                    if pole !=0:
                        den=den*np.poly1d([-1/pole,1])
                    else:
                        den=den*np.poly1d([1,0])
                for zero in temp_zeros:
                    if zero !=0:
                        num=num*np.poly1d([1/zero,1])
                    else:
                        num=num*np.poly1d([1,0])
            w,mag,phase=signal.bode(signal.TransferFunction(num,den))
            self.axis.clear()
            self.axis.set_title('Accumulative\'s Attenuation plot')
            self.axis.semilogx(w,-mag)
            self.axis.grid(color='grey',linestyle='-',linewidth=0.1)
            self.axis.set_xlabel("$Frequency [Hz]$")
            self.axis.set_ylabel("$Attenuation [dB]$")
            self.data_plot.draw()

            #ceros de la stage

            #return self.zerosOfStage, self.polesOfStage
        else:
            messagebox.showerror("Error", "No filter was created, plot cannot be realized")   

     def AddPole(self):

        self.val=self.comboPolos.get()
        if(len(self.val)==0):
            print("Ingresar Valor")
            return
        self.val2=complex(self.val)
        if(self.val2.imag==0):
            self.val2=self.val2.real
        print(self.val2) 
        if(self.val2 not in self.polesAux):
            print("Polo ya Seleccionado")
            return
        if(self.val2 in self.ListaDePolosPasados):  
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
        if(len(self.val)==0):
            print("Ingresar Valor")
            return
        self.val2=complex(self.val)
        if(self.val2.imag==0):
            self.val2=self.val2.real
        if(self.val2 not in self.ZerosAux):
            print("Cero ya Seleccionado")
            return
        if(self.val2 in self.ListaDeCerosPasados):  
            print("Polo ya Seleccionado")  
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
         #self.PolosNoEnStages = [item for item in temp1 if item not in temp2]   
         
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

           if(len(self.PolosSeleccionados)<len(self.ZerosSeleccionados)) :
               return
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
           self.ListaDePolosPasados = self.ListaDePolosPasados + self.den
           self.ListaDeCerosPasados = self.ListaDeCerosPasados + self.num
           self.TransferList.append(self.HdeStage)
           print(self.TransferList)
           print(self.ListaDePolosPasados)
           print(self.ListaDeCerosPasados)
           
           #Elimino los valores de la lsita y queda el combobox sin los valores de los polos seleccionados
           #Hago Clear de los polos y ceros seleccionados para que los proximos no acumulen los anteiores
           self.ZerosSeleccionados.clear()
           self.PolosSeleccionados.clear() 
           self.SelectedPoles.delete(0,END)
           self.SelectedZeros.delete(0,END) 

           self.GainOfStages.clear()                #borro todo para solamente tener las etapas actuales
           for i in range(0,len(self.TransferList)):
                self.GainOfStages.append(0)         #todas mis etapas arrancan con ganancia 0dB
                self.GainOfStagesVeces.append(1) 
           self.SelectedStage.insert(END, ["Stage", len(self.TransferList), "->", "Gain:", self.GainOfStages[len(self.TransferList)-1], "dB"])
          
         
     def DeleteStages(self):
       
         
         self.ListaDeCerosPasados.clear()
         self.ListaDePolosPasados.clear()
         
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
         self.filter_ready=False
         self.axis.clear()
         self.data_plot.draw()

     def AutoStages(self):
         
         if(len(self.TransferList)!=0):
            print("AutoGen ya Realizado, Delete Stages antes")
            return
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
         self.ListaDeCerosPasados.extend(self.Zeros)
         self.ListaDePolosPasados.extend(self.poles)

         self.den = []
         self.num = []
         
         print(self.polesAux) 
         
         for i in range(len(self.polesAux)) :
                if self.polesAux[i].imag == 0 :
                    self.polesAux[i]=self.polesAux[i].real
         
         for i in self.polesAux:
            if np.iscomplexobj(i):
                self.polesAux.remove(np.conjugate(i))
         
         self.MisPolos = sorted(self.polesAux,key=lambda x: ((np.sqrt(x.imag**2+x.real**2))/np.absolute(2*x.real)),reverse=True)       
         self.PolosAuxAux = []
         self.PolosAuxAux.extend(self.MisPolos)  

         for i in range(len(self.ZerosAux)):
                if self.ZerosAux[i].imag == 0 :
                    self.ZerosAux[i]=self.ZerosAux[i].real

         for i in self.ZerosAux:
            if np.iscomplexobj(i):
                self.ZerosAux.remove(np.conjugate(i))
         
         self.MisCeros = sorted(self.ZerosAux,key=lambda x:((np.sqrt(x.imag**2+x.real**2))/np.absolute(2*x.real)),reverse=True)       
         self.CerosAuxAux = []
         self.CerosAuxAux.extend(self.MisCeros)   
         self.HdeStage = []
         self.TransferList.clear() 
         
         print(self.MisPolos)
         print(self.MisCeros)
         
         counter = 0
         for i in self.MisCeros :
            if i == 0 :
                counter = counter +1 

         
         for i in range(int(counter/2)):
             if (0 in self.MisPolos) :
                self.den=[0,0]
                self.num=[0,0]
                self.MisCeros.remove(self.num[0])
                self.MisCeros.remove(self.num[0])
                self.MisPolos.remove(self.den[0])
                self.MisPolos.remove(self.den[0]) 
                self.HdeStage.append(self.den) 
                self.HdeStage.append(self.num)        
        
                self.TransferList.append(self.HdeStage)
                self.HdeStage = []
         
         for i in self.CerosAuxAux :
            self.den = []
            self.num = []
            if np.iscomplexobj(i) :
            #Si tengo un cero conjugado lo asocio a un par de polos conjugados
                for k in self.MisPolos :
                    if np.iscomplexobj(k) :
                        self.den=[k,np.conjugate(k)]
                        self.num=[i,np.conjugate(i)]
                        self.MisPolos.remove(k)
                        self.MisCeros.remove(i)
                        break
            else:
            #Si tengo un cero simple lo asocio a un polo simple
                for k in self.MisPolos :
                    if np.isrealobj(k) :
                        self.den=[k]
                        self.num=[i]
                        self.MisPolos.remove(k)
                        self.MisCeros.remove(i)
                        break
            print("Transfers del Primer For")
            print(self.den)
            print(self.num)


            if(len(self.den)!=0 or len(self.num)!=0):
                self.HdeStage.append(self.den) 
                self.HdeStage.append(self.num)        
       
                self.TransferList.append(self.HdeStage)
            self.HdeStage = []

        #Si tengo un par de polos conjugados le asocio dos ceros reales
        #Notar ya no quedn ceros conjugados, y quedan tantos polos simples como PolosSimplesTotales - CerosSimples
        #Notar que si quedan ceros simples no quedan polos simples por lo que asocio esos ceros a conjugados en el codigo este
         for i in self.PolosAuxAux :
            self.den = []
            self.num = []

            if np.iscomplexobj(i) :
            #Si tengo un cero conjugado lo asocio a un par de polos conjugados
                for k in range (len(self.MisCeros)) :
                        if (len(self.MisCeros) > 1) and (i in self.MisPolos) :
                            self.den=[i,np.conjugate(i)]
                            self.num=[self.MisCeros[k],self.MisCeros[k+1]]
                            self.MisCeros.remove(self.num[0])
                            self.MisCeros.remove(self.num[1])
                            self.MisPolos.remove(i)
                            
                            break
                        elif len(self.MisCeros) == 1 and (i in self.MisPolos) :
                            self.den=[i,np.conjugate(i)]
                            self.num=[self.MisCeros[k]]
                            self.MisCeros.remove(self.num[0])
                            self.MisPolos.remove(i)
                            break
                print("Transfers del Segundo For")
                print(self.den)
                print(self.num)
                if(len(self.den)!=0 or len(self.num)!=0):
                    self.HdeStage.append(self.den) 
                    self.HdeStage.append(self.num)        
       
                    self.TransferList.append(self.HdeStage)
                self.HdeStage = []
            #Si llegue hasta aca tengo:
            #Si tengo Ceros conjugados y polos simples tengo que asociar los polos simples al par conjugados
         for i in self.MisCeros :
            self.den = []
            self.num = []
            if np.iscomplexobj(i)  :
                for k in range (len(self.MisPolos)) :
                    if np.isrealobj(k) :
                        self.den=[self.MisPolos[k],self.MisPolos[k+1]]
                        self.num=[i,np.conjugate(i)]
                        self.MisPolos.remove(self.den[0])
                        self.MisPolos.remove(self.den[1])
                        self.MisCeros.remove(i)
                        break    

                if(len(self.den)!=0 or len(self.num)!=0):
                    self.HdeStage.append(self.den) 
                    self.HdeStage.append(self.num)        
       
                    self.TransferList.append(self.HdeStage)
                self.HdeStage = []    

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
             print("Transfers del While")
             print(self.den)
             print(self.num)
             
             self.HdeStage.append(self.den) 
             self.HdeStage.append(self.num)        
       
             self.TransferList.append(self.HdeStage)
             self.HdeStage = []

         print(self.TransferList)
         self.GainOfStages.clear()  
         self.GainOfStagesVeces.clear()
         for i in range(0,len(self.TransferList)): 
             self.GainOfStages.append(0)         #todas mis etapas arrancan con ganancia 0dB
             self.SelectedStage.insert(END, ["Stage", i+1 , "->" ,"Gain:", self.GainOfStages[i], "dB"])

        #----------------Buttons functions----------------------------
     
     def TransferOfStage(self):
        #a,b=self.getMax()
        self.GainOfStagesVeces.clear()
        self.detectStage = self.SelectedStage.curselection()
        for i in self.detectStage:
            self.stageIs = self.TransferList[i]
            #polos de la stage
        aux=self.gain_entry.get()   #aca guardo lo que ingreso el usuario
        self.gain_entry.delete(0, END) #una vez que lo tome lo elimino
        if len(aux) > 0:        #si el usuario introdujo algo
            #aca deberia parsear
            self.GainOfStages[i]=aux
            self.SelectedStage.delete(i)    #elimino la stage con ganancia vieja
            self.SelectedStage.insert(i, ["Stage", i+1 , "->" ,"Gain:", self.GainOfStages[i], "dB"]) #nueva stage con nueva ganancia

        K=np.power(10,float(self.GainOfStages[i])/20) #constante para obtener la ganancia en dB que pide el usuario    
        self.polesOfStage = self.stageIs[0] #siempre el primer elemento son los polos
        #ceros de la stage
        self.zerosOfStage = self.stageIs[1] #siempre el segundo elemento son los ceros
        #return self.zerosOfStage, self.polesOfStage
        self.den=np.poly1d([1])
        self.num=np.poly1d([1])
        for i in range(0,len(self.polesOfStage)):
            if self.polesOfStage[i] !=0:
                pol = np.poly1d([-1/(self.polesOfStage[i]), 1])
            else:
                pol = np.poly1d([1,0])
            self.den= self.den*pol
        for i in range(0,len(self.zerosOfStage)):
            if self.zerosOfStage[i] !=0:
                pol = np.poly1d([1/(self.zerosOfStage[i]), 1])
            else:
                pol = np.poly1d([1,0])
            self.num= self.num*pol
        self.num = self.num*K    
        self.TFofStage = signal.TransferFunction(self.num,self.den)
        self.w,self.mag,self.phase = signal.bode(self.TFofStage)
        self.stepT,self.step_mag=signal.step(self.TFofStage,N=1000)
        self.impT,self.imp_mag=signal.impulse(self.TFofStage,N=1000)
        dphase=np.ediff1d(self.phase)#Phase is in deegres
        dw=np.ediff1d(self.w)
        gd=-dphase/dw
        gd=np.append(gd,gd[len(gd)-1])
        self.group_delay= gd/(180/np.pi)
        self.w=self.w/(2*np.pi)
        self.filter_ready=True
        for i in range(0,len(self.GainOfStages)) : 
            kveces = np.power(10,float(self.GainOfStages[i])/20)
            self.GainOfStagesVeces.append(kveces)   
        self.CalcRD(self.GainOfStagesVeces)
        #Actualizo valor de RD
        print(str(self.RD))
        self.RDT.set("RD: " + str(self.RD))

     def onselect(self,evt):
        # Note here that Tkinter passes an event object to onselect()
        w = evt.widget
        index = int(w.curselection()[0])
        value = w.get(index)
        self.TransferOfStage()
        self.plot_atten()
      
      
     def CalcRD(self,ListaDeConstantes):
        Maximos = self.getMax()
        MaximosVeces = []
        for i in range(len(Maximos)) :
            M = np.power(10,float(Maximos[i])/20)
            MaximosVeces.append(M)
        Vsat = 14
         
        Vmin = 0.01 
        #armo un arreglo1 lleno de Vmax de cantidad dependinedo de la cantidad de etapas
        self.arreglo1= [Vsat for number in range(len(ListaDeConstantes))]
        #armo un arreglo2 lleno de vmin de cantidad dependiendo de la cantidad de etapas
        self.arreglo2 = [Vmin for number in range(len(ListaDeConstantes))]
        ConstantesDeGanancia= []
        ConstantesDeGanancia = ListaDeConstantes
        for i in range(len(ConstantesDeGanancia)) :
           for k in range(0,i+1) : 
               self.arreglo1[i]=self.arreglo1[i]/(ConstantesDeGanancia[k]*MaximosVeces[k])

        for i in range(len(ConstantesDeGanancia)):
            for k in range(0,i+1): 
               self.arreglo2[i]=self.arreglo2[i]/ConstantesDeGanancia[k]

        Vmax= min(self.arreglo1)
        Vmin = max(self.arreglo2)    
        self.RD = round(20*np.log10(Vmax/Vmin),2)
        print("valores Maximos")
        print(MaximosVeces)
        print(ConstantesDeGanancia)
        print("Arreglo de Ganancias")
        print(self.arreglo1)
        print("Arreglo de Atenuaciones")
        print(self.arreglo2)

     def getMax(self):
        maxs=[]
        for transfer in self.TransferList:
            den=np.poly1d([1])
            num=np.poly1d([1])
            temppoles=transfer[0]
            tempzeros=transfer[1]
            for pole in temppoles:
                if pole != 0:
                    den=den*np.poly1d([-1/pole,1])
                else:
                    den=den*np.poly1d([1,0])
            for zero in tempzeros:
                if zero != 0:
                    num=num*np.poly1d([1/zero,1])
                else:
                    num=num*np.poly1d([1,0])
            w,mag,phase=signal.bode(signal.TransferFunction(num,den))
            tempmax=0
            for m in mag:
                if m>tempmax:
                    tempmax=m
            maxs.append(tempmax)
        return maxs

if __name__ == "__main__":
    ex = TP4()