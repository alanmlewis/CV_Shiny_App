from shiny import App, render, ui, reactive
import numpy as np
from descriptions import desc
import matplotlib.pyplot as plt
from ratelimit import debounce

choices = {'o':'Oxidised Species','r':'Reduced Species'}
peak_choices = {'pp':'Peak Positions', 'ph':'Peak Height'}

app_ui = ui.page_fluid(
    ui.div(
    {"style": "text-align: center;"},ui.markdown("# **Simulating Cyclic Voltametry**")),
    ui.row(
        ui.column(3, ui.markdown("### Experimental Parameters"),
        ui.tooltip(ui.input_slider("T", "Temperature (K)", 100, 500, 298, step=1),desc['T']),
        ui.tooltip(ui.input_slider("E1", "Initial Potential (V)", -3.0000, 3.000, 0.0, step=0.01),desc['E1']),
        ui.tooltip(ui.input_slider("E2", "Switching Potential (V)", -3.000, 3.000, 0.8, step=0.01),desc['E2']),
        ui.tooltip(ui.input_slider("Crbulk", "Initial Concentration of Reduced Species (mM/dm\u00b3)", value=1,min=0.1,max=5,step=0.1),desc['Crbulk']),
        ui.tooltip(ui.input_slider("Coxbulk", "Initial Concentration of Oxidised Species (mM/dm\u00b3)", value=0.0,min=0.0,max=5,step=0.1),desc['Crbulk']),
        ui.tooltip(ui.input_text("v", "Scan Rate (mV/s)", 100),desc['v']),
        ui.tooltip(ui.input_checkbox_group('peaks', 'Show:', peak_choices),desc['peaks']),
        ui.input_radio_buttons('C_spe','Show concentration of:',choices)
        ),
        ui.column(7, ui.output_plot("plots",height='1500px')),
        ui.column(2, ui.markdown("### System Parameters"),
        ui.tooltip(ui.input_slider("E0", "Redox Couple Potential (V)", -2.5000, 2.5000, 0.496, step=0.01),desc['E0']),
        ui.tooltip(ui.input_text("k0", "Standard Rate Constant (cm/s)", 1.0),desc['k0']),
        ui.tooltip(ui.input_text("D", "Diffusion Coefficient (cm\u00b2/s)", 0.00005),desc['D']),
        ui.tooltip(ui.input_numeric('a','Transfer coefficient',0.5,min=0.1,max=0.9,step=0.05),desc['a']),
        ui.tooltip(ui.input_numeric('n','Number of electrons transfered',1,min=1,max=3,step=1),desc['n']),
        ui.tooltip(ui.input_numeric('dC','Capacitance',0,min=0,max=1,step=0.01),'test')
        ) 
    )
)

# Simulating the reversible reaction R  <--> O + ne-  using Butler-Volmer equation 
# Diescrete model of the solution adjacent to an electrode

###################### Physical constants (Do not change) #############################################

F = 96485           # Faraday constant, unit: C/mol
R = 8.31451         # ideal gas constant, unit: J/(K·mol)

A = 1              # surface area, unit: cm2
#Coxbulk = 0        # bulk concentration of oxidized species O, unit: mol/cm3
#a = 0.5            # alpha, transfer coefficient
#n = 1              # number of electrons in the redox reaction

nx = 100
nt = 800

####################### Simulation parameters to change ##############################################

def server(input, output, session):
    @debounce(0.5)
    @reactive.Calc
    def simulate():

        E0 = input.E0()
        E1 = input.E1()
        E2 = input.E2()
        v = float(input.v())
        T = input.T()
        k0 = float(input.k0())
        D = float(input.D())
        a = input.a()
        n = input.n()
        dC = input.dC()
        Crbulk = float(input.Crbulk())
        Crbulk *= 1e-9
        Coxbulk = float(input.Coxbulk())
        Coxbulk *= 1e-9

        ttot=abs(E2-E1)*2/v      # total time of full potential sweep (forward and back)
#        nt=int(abs(E2-E1)/dE)    # number of time/potential intervals in one potential sweep
        dE=abs(E2-E1)/nt
        dt=ttot/2*nt+1           # time increment of a diffusion grid

        x= 6*(D*ttot)**0.5       # diffusion distance
        dx= x/nx                 # position increment of a diffusion grid

        tini = 0                 # initial time for sweep
        xini = 0                 # initial position

        ns = 2                   # number of sweeps (forward and reverse)
        lamda = nx*nx/(36*nt)    # lambda, model diffusion coefficient

        # for each sweep, test if the sweep number is odd or even and construct the respective sweep
        E = np.array([])    # potential array
        for i in range(1,ns+1):
            if (i%2 == 1):
                E = np.append(E, np.linspace(E1, E2, nt+1)) # odd
            else:
                E = np.append(E, np.linspace(E2-dE, E1, nt)) # even

                   # redoxfine reaction R  <--> O + ne-
        
        kf = k0*np.exp((1-a)*n*F*(E-E0)/(R*T)) #rate constant for forward oxidation reaction
        kb = k0*np.exp((-a)*n*F*(E-E0)/(R*T)) #rate constant for backward reduction reaction

            #Initial concentration of the reduced species R
        
        M=nx         # number of elements in the diffusion grid for distance
        N=2*nt+1     # number of elements in numpy array E#

        Cr = np.zeros((N, M)) #concentration of oxidized species O is a 2D array
        Cr[0,:] = Crbulk # At t=0, bulk concentrations at all x
        Cr[:,0] = Crbulk # At x=0 (far from electrode), bulk concentration at all t
        
              #Initial concentration of the oxidized species O
        Cox = np.zeros((N, M))
        Cox[0,:] = Coxbulk # At t=0, bulk concentrations at all x
        Cox[:,0] = Coxbulk # At x=0 (far from electrode), bulk concentration at all t
               
            #Fluxes of reduced and oxidized species
        Jr = np.zeros(N)
        Jox = np.zeros(N)

        for k in range (0, N):

            Jr[k]= -(kf[k]*Cr[k,M-2]-kb[k]*Cox[k, M-2])/(1+(kf[k]*dx/D)+(kb[k]*dx)/D)
            Jox[k]= -Jr[k]
            if k>0: Cr[k, M-1]= Cr[k, M-2]+Jr[k]*dx/D # Calculate concentration at electrode for t>0
            if k>0: Cox[k, M-1]=Cox[k, M-2]+Jox[k]*dx/D # Calculate concentration at electrode for t>0
            if k==N-1: continue
        
        # Calculate concentrations at the next time/potential step
            for i in range(1, M-1):
                Cr[k+1, i] = lamda*(Cr[k, i-1] - (2)*Cr[k,i] + Cr[k,i+1])+Cr[k,i]
                Cox[k+1, i] = lamda*(Cox[k, i-1] - (2)*Cox[k,i] + Cox[k,i+1])+Cox[k,i]

             # calculating total faradaic current density
        jfar = -n*F*Jr*A*1e6
        jnf = np.ones(2*nt+1)*A*v*dC
        jnf[nt:] *= -1
        
        jtot = jfar + jnf
        
        return E,jtot,Cr,Cox,dx,dt,nt
        
    @output
    @render.plot
    def plots():
            E, jtot, Cr, Cox, dx, dt, nt = simulate()
            plt.rcParams.update({'font.size': 15}) #change the font size of scale bar
            plt.rcParams["font.family"] = "DejaVu Sans" # change the font
            fig,axs = plt.subplots(3)
#            fig.set_figheight(150)
            axs[0].plot(E, jtot, color = 'green' )
            axs[0].set_xlabel("Potential (V)")
            axs[0].set_ylabel("Current density ($\mu$A/cm$^2$)")
            axs[0].set_xlim([E[0],E[nt+1]])
            #legend = plt.legend(loc='upper left', fontsize= '13') #change the legend position and font size

            xs = np.arange(Cr.shape[1])[::-1]*dx
            ts = np.arange(Cr.shape[0])*dt
            nx = int(len(xs)/2)
            if input.C_spe()=='r':
                C = Cr*1e9
            else:
                C = Cox*1e9
            
            if not input.peaks() == ():
                peak_pos = []
                peak_j = []
                peak_pos.append(E[np.argmax(jtot)])
                peak_j.append(np.max(jtot))
                peak_pos.append(E[np.argmin(jtot)])
                peak_j.append(np.min(jtot))
            if 'pp' in input.peaks():
                for i in range(2):
                    axs[0].plot(peak_pos[i],peak_j[i],'ko')
                    axs[0].annotate('E = '+str(np.round(peak_pos[i],3))+' V',(peak_pos[i]*1.02,peak_j[i]),va = 'center')
            if 'ph' in input.peaks():
                label = ['$j_{pa}$ = ', '$j_{pc}$ = ']
                for i in range(2):
                    r = 1.0
                    slope_len = 0
                    init = i*nt
                    while abs(r) > 0.99:
                        slope_len +=5
                        r = np.corrcoef(E[init:init+slope_len],jtot[init:init+slope_len])[0,1]
                        if slope_len == nt:
                            print('FAIL')
                            break
                    slope_len -=5
                    x = E[init:init+slope_len]
                    m,c = np.linalg.lstsq(np.vstack([x,np.ones(len(x))]).T,jtot[init:init+slope_len],rcond=None)[0]

                    axs[0].plot(E,m*E+c,linestyle='dotted',color='black')
                    axs[0].plot([peak_pos[i],peak_pos[i]],[peak_j[i],m*peak_pos[i]+c],color='black')
                    height = abs(peak_j[i]-(m*peak_pos[i]+c))
                    ypos = min(peak_j[i],m*peak_pos[i]+c) + height/2
                    axs[0].annotate(label[i]+str(np.round(height,2))+' $\mu$A/cm$^2$',(peak_pos[i]*1.02,ypos),va='center')
                    
                
                
            a = axs[1].pcolor(ts[:nt+1],xs[-nx:],C[:nt+1,-nx:].T,cmap='Reds')
            axs[2].pcolor(ts[nt:],xs[-nx:],C[nt:,-nx:].T,cmap='Reds')
            axs[2].set_xlim(axs[2].get_xlim()[::-1])
            
            cb1 = fig.colorbar(a,orientation='vertical',ax=axs[1])
            cb2 = fig.colorbar(a,orientation='vertical',ax=axs[2])
            cb = fig.colorbar(a,orientation='vertical',ax=axs[0],shrink=0.001)
            cb.set_ticks([])
            plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
#            cb.remove()
            
            axs[1].set_xlabel("Time (s)")
            axs[2].set_xlabel("Time (s)")
            axs[1].set_ylabel("Distance from electrode (cm)")
            axs[2].set_ylabel("Distance from electrode (cm)")
            axs[1].set_title("Concentration during forward sweep (mM/dm$^3$)\n (Time increases from left to right)")
            axs[2].set_title("Concentration during reverse sweep (mM/dm$^3$)\n (Time increases from right to left)")


app = App(app_ui, server)
