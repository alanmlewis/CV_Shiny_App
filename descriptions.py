from shiny.ui import markdown
desc = {'T':'Temperature in Kelvin.',
'E1' : markdown('This is the starting point for potential scans. Also known as starting potential (E<sub>start</sub> or E<sub>1</sub>).'), 
'E2':  markdown('This is the end point for potential scan before switching the direction of scanning. Also known as end potential (E<sub>end</sub> or E<sub>2</sub>).'),
'Crbulk' : 'Concentration of electroactive species in the solution.',
'v' : 'Rate at which the voltage or potential is swept or scanned per unit time. Usually expressed in V/s or mV/s.',
'E0' :  markdown('The redox potential of a couple, denoted as E<sub>0</sub> or E<sub>1/2</sub>, reflects the characteristic electron-transfer behaviour of a species or couple under specific experimental conditions. It is also known as the formal potential under those conditions.'),
'k0' :  markdown('This represents the rate constant for transferring electrons between the electrode and the redox species (across the electrode-electrolyte interface) and is also known as heterogeneous electron transfer rate constant. It is simply a measure of the facility of a redox couple. A system with large k<sub>0</sub> means faster electron transfer kinetics, whereas a system with smaller k<sub>0</sub> will be sluggish.'),
'D' :  markdown('The rate of diffusion of a redox species under a unit concentration gradient (dC/dx) and is generally expressed in units of cm<sup>2</sup>/s. This parameter is essential in understanding how quickly a species can migrate through a medium when there is a concentration gradient. Various factors influence the diffusion coefficient (D) of a species, with the size and shape of the species being particularly significant.'),
'a' : 'Generally a constant (a fraction) and is a measure of the symmetry of the energy barrier at the transition state for electron transfer within the normal potential range studied. Experimentally it varies anywhere from 0.3 to 0.7. A value of transfer coefficient 0.5 represents a symmetrical barrier.',
'n' : 'Number of electrons involved in a redox process.',
'peaks' : markdown('The Peak Position is the potential at which rate of oxidation or reduction is maximum. The Peak Height is the difference in current from baseline to the peak. This is known as *j<sub>pa</sub>* and *j<sub>pc</sub>* for oxidation and reduction respectively. For a non-zero baseline, the peak current might not be the same as peak height.')
}
