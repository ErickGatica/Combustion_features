# Importing libraries
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot  as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import numpy as np 

from hapi import *

# Some physical constantes
h=6.62607015e-27 # Planck constant in erg*second
c=2.99792458e10 # Speed of light in cm/second
k=1.380649e-16 # Boltzmann constant in erg/K
c2=1.4387769 # Second radiation constant in cm*K
Tref=296 # Reference temperature in K



# Dictionary with the name of molecule and the respective ID
molecule_id_dict = {
    "H2O": 1,
    "CO2": 2,
    "CO": 5,
    "N2": 22,
    "O2": 7,
    "CH4": 6,
    "C2H6": 27,
    "H2": 45,
    "NO":8,
    "NO2":10
}

# Dictionary Data_table
Data_table = {}

# Functions
def Getting_table():
    global Data_table
    # Getting the data
    molecule_name=molecule_id_var.get()
    molecule_id=molecule_id_dict[molecule_name]
    isotopo_id=1
    nu_min=float(numin.get())
    nu_max=float(numax.get())
    # This is to get the table of linestrenghts and data
    # Defining the molecule
    #molecule=moleculeName(molecule_id)
    fetch(molecule_name, molecule_id, isotopo_id, nu_min, nu_max)
    # Selecting data
    select(molecule_name,ParameterNames=('nu','sw','elower'),Conditions=('between','nu',nu_min,nu_max))
    nu_plot,sw_plot,e_lower_plot = getColumns(molecule_name,['nu','sw','elower'])
    # array to save the ID of the linestrengths, same length than sw_plot
    Id_selection=np.arange(0,len(sw_plot),1)
    # Creating dictionary to save the data
    Data_table = {
        "selection_id_sw": Id_selection,
        "nu_plot": nu_plot,
        "sw_plot": sw_plot,
        "e_lower": e_lower_plot
    }
    print('Data obtained',Data_table)
    # Filtering the Data values to avoid linestrengths less than 1E-27, so re defining the Data_table
    value_to_filter=1E-24
    Data_table = {
        "selection_id_sw": Id_selection[sw_plot>value_to_filter],
        "nu_plot": nu_plot[sw_plot>value_to_filter],
        "sw_plot": sw_plot[sw_plot>value_to_filter],
        "e_lower": e_lower_plot[sw_plot>value_to_filter]
    }
    nu_plot=Data_table['nu_plot']
    sw_plot=Data_table['sw_plot']
    e_lower_plot=Data_table['e_lower']
    Id_selection=np.arange(0,len(sw_plot),1)
    Data_table["Id_selection"]=Id_selection
    print('Filtered Data obtained ',Data_table)
    # Creating a table for ID, Nu, sw and e_lower using Id_selection, nu_plot, sw_plot and e_lower_plot
    # Creating the table in the ax1, header in first column and data goes vertically
    ax1.axis('tight')
    ax1.axis('off')
    cellText=[['ID','Nu','sw','e_lower']]+[[Id_selection[i],nu_plot[i],sw_plot[i],e_lower_plot[i]] for i in range(0,len(sw_plot))]
    table=ax1.table(cellText,loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)
    # Redraw the plot
    canvas.draw()
    return Data_table



def Plotting():
    global Data_table
    # Parameters 
    ID_selected_S=int(linestrengthID.get())
    ID_selected_S1=int(linestrengthID1.get())
    molecule_name=molecule_id_var.get()
    molecule_id=molecule_id_dict[molecule_name]
    isotopo_id=1
    T_min=float(Tmin.get())
    T_max=float(Tmax.get())
    P=float(Pr.get())
    length=float(lengthh.get())
    molar_f=float(molar_fraction.get())
    Tref=296 # Reference temperature in K
    R_g=8.314 # Gas constant in J/(mol*K)
    Avogadro=6.022e23 # Avogadro number
    # Getting the wavenumber of the transition and the line strength
    Nu_plot=Data_table['nu_plot'] #cm^-1
    S_plot=Data_table['sw_plot']
    Eprime_plot=Data_table['e_lower']

    # This is fot the first linestrength
    v_ij=Nu_plot[ID_selected_S] #wavenumber of the transition cm^-1
    S_ij_ref=S_plot[ID_selected_S] #Line strength cm^-1/(molecule*cm^-2)
    E_pp=Eprime_plot[ID_selected_S] # Lower state energy of the transition cm^-1

    # This is for the second linestrength
    v_ij1=Nu_plot[ID_selected_S1] #wavenumber of the transition cm^-1
    S_ij_ref1=S_plot[ID_selected_S1] #Line strength cm^-1/(molecule*cm^-2)
    E_pp1=Eprime_plot[ID_selected_S1] # Lower state energy of the transition cm^-1

    # Partion function
    Q_part_ref_cal=partitionSum(molecule_id,isotopo_id, [Tref])  # Getting the partition function at the reference temperature
    Q_part_ref=Q_part_ref_cal[0]

    # Defining a temperature range and a vector to save the S_ij values
    Temps=np.linspace(T_min,T_max,100)
    S_ij_values=np.zeros(len(Temps)) # first linestrenght  
    S_ij_values1=np.zeros(len(Temps)) # second linestrenght

    # Computing the S_ij values
    for i in range(0,len(Temps)):
        # Computing values for the first linestrength 
        Q_part_HT=partitionSum(molecule_id,isotopo_id, [Temps[i]])
        Q_part=Q_part_HT[0]
        S_cal=S_ij_ref*Q_part_ref/Q_part * np.exp(-c2*E_pp/Temps[i])/np.exp(-c2*E_pp/Tref) * (1-np.exp(-c2*v_ij/Temps[i]))/(1-np.exp(-c2*v_ij/Tref))
        S_ij_values[i]=S_cal

        # Computing values for the second linestrength
        Q_part_HT1=partitionSum(molecule_id,isotopo_id, [Temps[i]])
        Q_part1=Q_part_HT1[0]
        S_cal1=S_ij_ref1*Q_part_ref/Q_part1 * np.exp(-c2*E_pp1/Temps[i])/np.exp(-c2*E_pp1/Tref) * (1-np.exp(-c2*v_ij1/Temps[i]))/(1-np.exp(-c2*v_ij1/Tref))
        S_ij_values1[i]=S_cal1
        
    # Converting the S_ij values to absorbance using pressure, molar fraction and path lengths

    for i in range(0,len(Temps)):
        S_ij_values[i]=S_ij_values[i]*molar_f*P*101325/(R_g*Temps[i])/(100**3)*length*Avogadro
        S_ij_values1[i]=S_ij_values1[i]*molar_f*P*101325/(R_g*Temps[i])/(100**3)*length*Avogadro

    # Getting the ratio of the S2/S1
    S_ratios=np.zeros(len(Temps))
    for i in range(0,len(Temps)):
        S_ratios[i]=S_ij_values1[i]/S_ij_values[i]

    # Plotting the S_ij values in function of the temperature
    ax2.plot(Temps, S_ratios, label='S2/S1'+'Nu='+str(v_ij1)+' and '+str(v_ij)+'cm^-1')
    ax2.set_xlabel('Temperature [K]')
    ax2.set_ylabel('Line intensity S_ij')
    ax2.set_title('Temperature dependence of the line intensity')
    ax2.legend()
    # Redraw the plot
    canvas.draw()

    # Plotting the S_ij values in function of the temperature
    ax3.plot(Temps, S_ij_values, label='S1 '+ 'Nu= '+str(v_ij)+' cm^-1')
    ax3.set_xlabel('Temperature [K]')
    ax3.set_ylabel('Line intensity S_ij')
    ax3.set_title('Temperature dependence of the line intensity')
    ax3.legend()

    # Plotting the S_ij values in function of the temperature
    ax4.plot(Temps, S_ij_values1, label='S2 '+ 'Nu= '+str(v_ij1)+' cm^-1')
    ax4.set_xlabel('Temperature [K]')
    ax4.set_ylabel('Line intensity S_ij')
    ax4.set_title('Temperature dependence of the line intensity')
    ax4.legend()
    # Redraw the plot
    canvas2.draw()



    return Temps,S_ratios



def clear_plot():
    # Function to clear plots
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()
    # Redraw the plot
    canvas.draw()

############### INPUTS ############################

# Creating the frame
root=tk.Tk()
root.title('Thermometry')
root.geometry('1080x1080')
root['bg']='white'

# Tab manager and notebook
notebook = ttk.Notebook(root)

# First tab
tab1 = ttk.Frame(notebook)
notebook.add(tab1, text='Thermometry')

# Intro label
intro_label=ttk.Label(tab1, text="Welcome to the GUI to obtain the spectral of the combustion products")
intro_label.grid(sticky='nw', padx=10, pady=10, columnspan=2, row=0, column=0)

# Adding the notebook to the window
notebook.grid(sticky='nsew')

# Creating frame for inputs
inputs_frame = ttk.LabelFrame(tab1, text='Inputs')
inputs_frame.grid(row=1, column=0, columnspan=2, padx=10, pady=10)

# Molecule of interest multipe choice, considering dictionary of combustion products
molecule_id_label = ttk.Label(inputs_frame, text="Molecule of interest:")
molecule_id_label.grid(row=0, column=0)
molecule_id_var = tk.StringVar()
molecule_id_dropdown = ttk.Combobox(inputs_frame, textvariable=molecule_id_var, values=["H2O", "CO2", "CO", "N2", "O2","CH4","H2","NO","NO2"]) 
molecule_id_dropdown.grid(row=0, column=1)
molecule_id_dropdown.set("H2O")  # Set a default value

# Temperature range
Tmin_label = ttk.Label(inputs_frame, text="Minimum temperature [K]:")
Tmin_label.grid(row=1, column=0)
Tmin=tk.Entry(inputs_frame)
Tmin.grid(row=1, column=1)

Tmax_label = ttk.Label(inputs_frame, text="Maximum temperature [K]:")
Tmax_label.grid(row=2, column=0)
Tmax=tk.Entry(inputs_frame)
Tmax.grid(row=2, column=1)

# Pressure 
Pr_label = ttk.Label(inputs_frame, text="Pressure [atm]:")
Pr_label.grid(row=3, column=0)
Pr=tk.Entry(inputs_frame)
Pr.grid(row=3, column=1)

# Molar fraction
molar_fraction_label = ttk.Label(inputs_frame, text="Molar fraction of the molecule:")
molar_fraction_label.grid(row=4, column=0)
molar_fraction=tk.Entry(inputs_frame)
molar_fraction.grid(row=4, column=1)

# Path length
lengthh_label = ttk.Label(inputs_frame, text="Path length [cm]:")
lengthh_label.grid(row=5, column=0)
lengthh=tk.Entry(inputs_frame)
lengthh.grid(row=5, column=1)

# Wavenumber range
numin_label = ttk.Label(inputs_frame, text="Minimum wavenumber [cm-1]:")
numin_label.grid(row=6, column=0)
numin=tk.Entry(inputs_frame)
numin.grid(row=6, column=1)

numax_label = ttk.Label(inputs_frame, text="Maximum wavenumber [cm-1]:")
numax_label.grid(row=7, column=0)
numax=tk.Entry(inputs_frame)
numax.grid(row=7, column=1)

# Button to get the data from HAPI and be able to see which linestrenghts are available
get_data_button = ttk.Button(inputs_frame, text="Get data", command=Getting_table)
get_data_button.grid(row=8, column=0, columnspan=2, padx=10, pady=10)

# Selecting which linestrength want to simulate, has to be an integer
linestrengthID_label = ttk.Label(inputs_frame, text="1st ID of Linestrength to compute S(T):")
linestrengthID_label.grid(row=9, column=0)
linestrengthID=tk.Entry(inputs_frame)
linestrengthID.grid(row=9, column=1)

# Selecting 2 linestrengths to get the ratio of the S2/S1
linestrengthID1_label = ttk.Label(inputs_frame, text="2nd ID of Linestrength to compute S(T):")
linestrengthID1_label.grid(row=10, column=0)
linestrengthID1=tk.Entry(inputs_frame)
linestrengthID1.grid(row=10, column=1)

# Button to plot the linestrength in function of the temperature
plot_button = ttk.Button(inputs_frame, text="Plot", command=Plotting)
plot_button.grid(row=11, column=0, columnspan=2, padx=10, pady=10)

# Button for clear the plots
clear_button = ttk.Button(inputs_frame, text="Clear the plots", command=clear_plot)
clear_button.grid(row=12, column=0, columnspan=2)

# Creating the plots for TAB 1
######### PLOTS #########
fig, ((ax1,ax2)) = plt.subplots(1, 2, figsize=(12,8))
fig.tight_layout()  # Automatically adjust spacing between subplots

# Creating the canvas for embedding the plot
canvas = FigureCanvasTkAgg(fig, master=tab1)  # Display the plot in tab1
canvas_widget = canvas.get_tk_widget()
canvas_widget.grid(row=13, column=0, columnspan=2)

# Create a frame for the Matplotlib toolbar
toolbar_frame = ttk.Frame(tab1)
toolbar_frame.grid(row=14, column=0, columnspan=2, sticky="nsew")

# Add a Matplotlib toolbar for zooming
toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
toolbar.update()
toolbar.pack(side=tk.BOTTOM, fill=tk.X)
canvas_widget.grid(row=15, column=0, columnspan=2, sticky=tk.NSEW)

# Configure grid row and column of plots weights to make them expand with the window
tab1.grid_rowconfigure(15, weight=1)
tab1.grid_columnconfigure(0, weight=1)

# Configure grid row and column of zoombar weights to make them expand with the window
tab1.grid_rowconfigure(15, weight=1)
tab1.grid_columnconfigure(0, weight=1)

#  Display the canvas
canvas.draw()

# Adding bar to go up and down in the plots frame
scrollbar = ttk.Scrollbar(tab1, orient="vertical", command=canvas.get_tk_widget().yview)
scrollbar.grid(row=15, column=2, sticky="ns")
canvas.get_tk_widget().configure(yscrollcommand=scrollbar.set)

# Second tab
tab2 = ttk.Frame(notebook)
notebook.add(tab2, text='Plot S1 and S2')

# Intro label
intro_label=ttk.Label(tab2, text="Plot of S1 and S2")
intro_label.grid(sticky='nw', padx=10, pady=10, columnspan=2, row=0, column=0)

# Creating the plots for TAB 2
######### PLOTS #########
fig2, ((ax3,ax4)) = plt.subplots(1, 2, figsize=(12,5))
fig2.tight_layout()  # Automatically adjust spacing between subplots

# Creating the canvas for embedding the plot
canvas2 = FigureCanvasTkAgg(fig2, master=tab2)  # Display the plot in tab1
canvas_widget2 = canvas2.get_tk_widget()
canvas_widget2.grid(row=1, column=0, columnspan=2)

# Create a frame for the Matplotlib toolbar
toolbar_frame2 = ttk.Frame(tab2)
toolbar_frame2.grid(row=2, column=0, columnspan=2, sticky="nsew")

# Add a Matplotlib toolbar for zooming
toolbar2 = NavigationToolbar2Tk(canvas2, toolbar_frame2)
toolbar2.update()
toolbar2.pack(side=tk.BOTTOM, fill=tk.X)
canvas_widget2.grid(row=3, column=0, columnspan=2, sticky=tk.NSEW)

# Configure grid row and column of plots weights to make them expand with the window
tab2.grid_rowconfigure(3, weight=1)
tab2.grid_columnconfigure(0, weight=1)

# Configure grid row and column of zoombar weights to make them expand with the window
tab2.grid_rowconfigure(3, weight=1)
tab2.grid_columnconfigure(0, weight=1)

#  Display the canvas
canvas2.draw()

# Adding bar to go up and down in the plots frame
scrollbar2 = ttk.Scrollbar(tab2, orient="vertical", command=canvas2.get_tk_widget().yview)
scrollbar2.grid(row=3, column=2, sticky="ns")
canvas2.get_tk_widget().configure(yscrollcommand=scrollbar2.set)



# loop
root.mainloop()