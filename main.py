# This is a code for computing the simulation features in the combustion products of an engine
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot  as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import numpy as np 

from spectrum import spectrum


# Importing defined file for calculations
from combustion import combustion


#Products as global variable
Products ={}

def update_plot():
    # Function to update the plot
    global Products
    # Updating variables
    T1=float(T.get())
    P1=float(P.get())
    XMethane1=float(XMethane.get())
    XEthane1=float(XEthane.get())
    XEthene1=float(XEthene.get())
    XPropane1=float(XPropane.get())
    XButhane1=float(XButhane.get())
    XHydrogen1=float(XHydrogen.get())
    XCO1=float(XCO.get())
    XCO21=float(XCO2.get())
    XN2=1-XMethane1-XEthane1-XEthene1-XPropane1-XButhane1-XHydrogen1-XCO1-XCO21
    Perc1=float(Perc.get())
    H2ble1=float(H2ble.get())
    
    # Defining a dictionary
    Compositions={
        'T':T1,
        'P':P1,
        'Methane':XMethane1,
        'Ethane':XEthane1,
        'Ethene':XEthene1,
        'Propane':XPropane1,
        'Buthane':XButhane1,
        'Hydrogen':XHydrogen1,
        'CO':XCO1,
        'CO2':XCO21,
        'N2':XN2,
        'Perc':Perc1,
        'H2ble':H2ble1
    }

    # Computing the combustion
    Products=combustion(Compositions)
    
    # Creating plot  bar plot with the products of the combustion outlet using the dictionary obtained from combustion function
    # Creating the bar plot with a label of the molar fraction of the products
    bars=ax1.bar(range(len(Products)), Products.values())
    #plt.bar_label(bars, labels=np.round(list(Products.values()), 3), label_type='edge', padding=3)
    ax1.set_title('Products of the combustion')
    ax1.set_xlabel('Products')
    ax1.set_ylabel('Molar fraction')
    ax1.set_xticks(range(len(Products)))  # Add this line
    ax1.set_xticklabels(Products.keys())
    

    # Create a table from the dictionary data
    table_data = list(zip(Products.keys(), ['{:.3e}'.format(val) for val in Products.values()])) #['{:.3e}'.format(val) for val in Products.values()]  np.round(list(Products.values()), 3)
    table = ax2.table(cellText=table_data, colLabels=['Products', 'Molar fraction'], loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.2)
    # Remove the x and y axis
    ax2.axis('off')

    # Redraw the plot
    canvas.draw()

    # Printing plots were updated
    print('Plots were updated')
    return Products


def clear_plot():
    # Function to clear plots
    ax1.clear()
    ax2.clear()
    # Redraw the plot
    canvas.draw()

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

def update_spectrum():
    # Access the global variable
    global Products
    # Function to update the spectrum
    # Getting the molecule
    P1 = float(P.get())
    T1 = float(T.get())
    molecule_id = molecule_id_var.get()
    molecule_id = molecule_id_dict[molecule_id]
    numin1 = float(numin.get())
    numax1 = float(numax.get())
    numin_cm = 1e7/numax1 # Convert to cm-1
    numax_cm = 1e7/numin1 # Convert to cm-1
    wave_step1 = float(wave_step.get())
    method1 = method_var.get()
    length1 = float(length.get())
    R_gas=8.314 # J/(mol K)
    total_mol=P1*101325/(R_gas*T1) # Total mol per volume m3
    Avogadro=6.022e23 # Avogadro number

    # Defining the dictionary of dilution for the computation of the spectrum
    if molecule_id == 1: # H2O
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':Products['X2_H2O'],'O2':Products['X2_O2'],'CO2':Products['X2_CO2'],'CO':Products['X2_CO'],'CH4':Products['X2_Methane'],'C2H6':Products['X2_Ethane'],'H2':Products['X2_H2'],'N2':N2_approx,'NO':Products['X2_NO'],'NO2':Products['X2_NO2']}
        Factor_absorption=Products['X2_H2O']*total_mol/(100**3)*length1*Avogadro
    elif molecule_id == 2: # CO2
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':Products['X2_CO2'],'O2':Products['X2_O2'],'H2O':Products['X2_H2O'],'CO':Products['X2_CO'],'CH4':Products['X2_Methane'],'C2H6':Products['X2_Ethane'],'H2':Products['X2_H2'],'N2':N2_approx,'NO':Products['X2_NO'],'NO2':Products['X2_NO2']}
        Factor_absorption=Products['X2_CO2']*total_mol/(100**3)*length1*Avogadro
    elif molecule_id == 5: # CO
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':Products['X2_CO'],'O2':Products['X2_O2'],'H2O':Products['X2_H2O'],'CO2':Products['X2_CO2'],'CH4':Products['X2_Methane'],'C2H6':Products['X2_Ethane'],'H2':Products['X2_H2'],'N2':N2_approx,'NO':Products['X2_NO'],'NO2':Products['X2_NO2']}
        Factor_absorption=Products['X2_CO']*total_mol/(100**3)*length1*Avogadro
    elif molecule_id == 22: # N2
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':N2_approx,'O2':Products['X2_O2'],'H2O':Products['X2_H2O'],'CO2':Products['X2_CO2'],'CO':Products['X2_CO'],'CH4':Products['X2_Methane'],'C2H6':Products['X2_Ethane'],'H2':Products['X2_H2'],'NO':Products['X2_NO'],'NO2':Products['X2_NO2']}
        Factor_absorption=N2_approx*total_mol/(100**3)*length1*Avogadro
    elif molecule_id == 7: # O2
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':Products['X2_O2'],'N2':N2_approx,'H2O':Products['X2_H2O'],'CO2':Products['X2_CO2'],'CO':Products['X2_CO'],'CH4':Products['X2_Methane'],'C2H6':Products['X2_Ethane'],'H2':Products['X2_H2'],'NO':Products['X2_NO'],'NO2':Products['X2_NO2']}
        Factor_absorption=Products['X2_O2']*total_mol/(100**3)*length1*Avogadro
    elif molecule_id == 6: # CH4
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':Products['X2_Methane'],'N2':N2_approx,'H2O':Products['X2_H2O'],'CO2':Products['X2_CO2'],'CO':Products['X2_CO'],'O2':Products['X2_O2'],'C2H6':Products['X2_Ethane'],'H2':Products['X2_H2'],'NO':Products['X2_NO'],'NO2':Products['X2_NO2']}
        Factor_absorption=Products['X2_Methane']*total_mol/(100**3)*length1*Avogadro
    elif molecule_id == 27: # C2H6
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':Products['X2_Ethane'],'N2':N2_approx,'H2O':Products['X2_H2O'],'CO2':Products['X2_CO2'],'CO':Products['X2_CO'],'O2':Products['X2_O2'],'CH4':Products['X2_Methane'],'H2':Products['X2_H2'],'NO':Products['X2_NO'],'NO2':Products['X2_NO2']}
        Factor_absorption=Products['X2_Ethane']*total_mol/(100**3)*length1*Avogadro
    elif molecule_id == 45: # H2
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':Products['X2_H2'],'N2':N2_approx,'H2O':Products['X2_H2O'],'CO2':Products['X2_CO2'],'CO':Products['X2_CO'],'O2':Products['X2_O2'],'CH4':Products['X2_Methane'],'C2H6':Products['X2_Ethane'],'NO':Products['X2_NO'],'NO2':Products['X2_NO2']}
        Factor_absorption=Products['X2_H2']*total_mol/(100**3)*length1*Avogadro
    elif molecule_id==8: # NO
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':Products['X2_NO'],'N2':N2_approx,'H2O':Products['X2_H2O'],'CO2':Products['X2_CO2'],'CO':Products['X2_CO'],'O2':Products['X2_O2'],'CH4':Products['X2_Methane'],'C2H6':Products['X2_Ethane'],'H2':Products['X2_H2'],'NO2':Products['X2_NO2']}
        Factor_absorption=Products['X2_NO']*total_mol/(100**3)*length1*Avogadro
    elif molecule_id==10: # NO2
        N2_approx=1-Products['X2_H2O']-Products['X2_CO2']-Products['X2_CO']-Products['X2_Methane']-Products['X2_Ethane']-Products['X2_H2']-Products['X2_O2']-Products['X2_NO']-Products['X2_NO2']
        Diluent={'self':Products['X2_NO2'],'N2':N2_approx,'H2O':Products['X2_H2O'],'CO2':Products['X2_CO2'],'CO':Products['X2_CO'],'O2':Products['X2_O2'],'CH4':Products['X2_Methane'],'C2H6':Products['X2_Ethane'],'H2':Products['X2_H2'],'NO':Products['X2_NO']}
        Factor_absorption=Products['X2_NO2']*total_mol/(100**3)*length1*Avogadro
    # Computing the spectrum
    Data=spectrum(P1,T1,length1,numin_cm,numax_cm,molecule_id,1,method1,wave_step1,Diluent)
    # Getting the data for plotting
    nu=10**7/Data.nu # getting the data in nm
    Absorption=np.multiply(Data.coef,Factor_absorption)
    name_isoto=Data.name_isoto
    # Plotting the spectrum
    ax3.plot(nu,Absorption, label=name_isoto+' $T$=' + str(T1) + 'K $P$=' + str(P1) + 'atm')
    ax3.set_title('Absorption Spectrum')
    ax3.set_xlabel('Wavelength [nm]')
    ax3.set_ylabel('Absorption')
    ax3.legend(loc='upper right',bbox_to_anchor=(1.0, 1))
    # Redraw the plot
    canvas2.draw()


    print('Spectrum was updated')

def clear_spectrum():
    # Function to clear the spectrum
    ax3.clear()
    # Redraw the plot
    canvas2.draw()


############### INPUTS ############################

#### TAB 1 ####
    
# Creating the frame 
root=tk.Tk()
root.title("Combustion reactants and products")
root['bg']='white'
root.geometry("1080x1080")
#intro_frame.grid(row=0, column=0,columnspan=2,sticky="snew")

# Creating the notebook (tab manager)
notebook = ttk.Notebook(root)

# Creating the first tab
tab1 = ttk.Frame(notebook)
notebook.add(tab1, text='Reactants and Products')


# Adding the description of the code in the window
intro_label=ttk.Label(tab1, text="Welcome to the GUI to obtain the spectral of the combustion products")
intro_label.grid(sticky='nw', padx=10, pady=10, columnspan=2, row=0, column=0)

# Creating the second tab
tab2 = ttk.Frame(notebook)
notebook.add(tab2, text='Spectrum')

# Adding the description of the code in the second tab
intro_label2 = ttk.Label(tab2, text="Here is for analyzing the spectrum of the combustion products")
intro_label2.grid(sticky='nw')


# Adding the notebook to the window
notebook.grid(sticky='nsew')

# Frame for inputs tab1
input_frame=ttk.LabelFrame(tab1,text="Input initial conditions of the fuel and air")
input_frame.grid(row=1, column=0,padx=10,pady=10, sticky="nw")
#Temperature input, 0
T_label=ttk.Label(input_frame,text="Temperature [K]:")
T_label.grid(row=0, column=0)
T=ttk.Entry(input_frame)
T.grid(row=0, column=1)
#Pressure input, 1
P_label=ttk.Label(input_frame,text="Pressure [atm]:")
P_label.grid(row=1, column=0)
P=ttk.Entry(input_frame)
P.grid(row=1, column=1)
#Methane molar fraction, 2
XMethane_label=ttk.Label(input_frame,text="Fuel Methane molar fraction: ")
XMethane_label.grid(row=2, column=0)
XMethane=ttk.Entry(input_frame)
XMethane.grid(row=2, column=1)
#Ethane molar fraction, 3
XEthane_label=ttk.Label(input_frame,text="Fuel Ethane molar fraction")
XEthane_label.grid(row=3, column=0)
XEthane=ttk.Entry(input_frame)
XEthane.grid(row=3, column=1)
#Ethene molar fraction, 4
XEthene_label=ttk.Label(input_frame,text="Fuel Ethene molar fraction")
XEthene_label.grid(row=4, column=0)
XEthene=ttk.Entry(input_frame)
XEthene.grid(row=4, column=1)
#Propane molar fraction
XPropane_label=ttk.Label(input_frame,text="Fuel Propane molar fraction")
XPropane_label.grid(row=5, column=0)
XPropane=ttk.Entry(input_frame)
XPropane.grid(row=5, column=1)
#Buthane molar fraction
XButhane_label=ttk.Label(input_frame,text="Fuel Buthane or more molar fraction")
XButhane_label.grid(row=6, column=0)
XButhane=ttk.Entry(input_frame)
XButhane.grid(row=6, column=1)
# Hydrogen molar fraction
XHydrogen_label=ttk.Label(input_frame,text="Fuel Hydrogen molar fraction")
XHydrogen_label.grid(row=7, column=0)
XHydrogen=ttk.Entry(input_frame)
XHydrogen.grid(row=7, column=1)
# CO molar fraction
XCO_label=ttk.Label(input_frame,text="Fuel Carbon monoxide molar fraction")
XCO_label.grid(row=8, column=0)
XCO=ttk.Entry(input_frame)
XCO.grid(row=8, column=1)
# CO2 molar fraction
XCO2_label=ttk.Label(input_frame,text="Fuel Carbon dioxide molar fraction")
XCO2_label.grid(row=9, column=0)
XCO2=ttk.Entry(input_frame)
XCO2.grid(row=9, column=1) 
# % of complete combustion
Perc_label=ttk.Label(input_frame,text="Percentage of complete combustion:")
Perc_label.grid(row=10, column=0)
Perc=ttk.Entry(input_frame)
Perc.grid(row=10, column=1)
# % of Hydrogen blending
H2ble_label=ttk.Label(input_frame,text="Percentage of Hydrogen blending:")
H2ble_label.grid(row=11, column=0)
H2ble=ttk.Entry(input_frame)
H2ble.grid(row=11, column=1)


# Creating button to compute 
compute_button=ttk.Button(input_frame,text="Compute",command=update_plot)
compute_button.grid(row=12,column=0, columnspan=2)

# Button for clear the plots
clear_button = ttk.Button(input_frame, text="Clear the plots", command=clear_plot)
clear_button.grid(row=13, column=0, columnspan=2)

# Creating the plots for TAB 1
######### PLOTS #########
fig, ((ax1), (ax2)) = plt.subplots(2, 1, figsize=(8,6 ))
fig.subplots_adjust(hspace=1.0, wspace=0.5)  # Adjust vertical and horizontal spacing between subplots

# Creating the canvas for embedding the plot
canvas = FigureCanvasTkAgg(fig, master=tab1)  # Display the plot in tab1
canvas_widget = canvas.get_tk_widget()
canvas_widget.grid(row=14, column=0, columnspan=2)

# Create a frame for the Matplotlib toolbar
toolbar_frame = ttk.Frame(tab1)
toolbar_frame.grid(row=15, column=0, columnspan=2, sticky="nsew")

# Add a Matplotlib toolbar for zooming
toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
toolbar.update()
toolbar.pack(side=tk.BOTTOM, fill=tk.X)
canvas_widget.grid(row=16, column=0, columnspan=2, sticky=tk.NSEW)

# Configure grid row and column of plots weights to make them expand with the window
tab1.grid_rowconfigure(16, weight=1)
tab1.grid_columnconfigure(0, weight=1)

# Configure grid row and column of zoombar weights to make them expand with the window
tab1.grid_rowconfigure(16, weight=1)
tab1.grid_columnconfigure(0, weight=1)

#  Display the canvas
canvas.draw()

# Adding bar to go up and down in the plots frame
scrollbar = ttk.Scrollbar(tab1, orient="vertical", command=canvas.get_tk_widget().yview)
scrollbar.grid(row=16, column=2, sticky="ns")
canvas.get_tk_widget().configure(yscrollcommand=scrollbar.set)



#### TAB 2 #### 

# Creating the frame for the second tab

# Frame for inputs tab2
input_frame2=ttk.LabelFrame(tab2,text="Inputs for simulating spectrum")
input_frame2.grid(row=1, column=0,padx=10,pady=10, sticky="nw")

# Path length input
length_label=ttk.Label(input_frame2,text="Path length [cm]:")
length_label.grid(row=0, column=0)
length=ttk.Entry(input_frame2)
length.grid(row=0, column=1)

# Max wavelength input
numax_label = ttk.Label(input_frame2, text="Max Wavelength [nm]:")
numax_label.grid(row=1, column=0)
numax = ttk.Entry(input_frame2)
numax.grid(row=1, column=1)

# Min wavelength input
numin_label = ttk.Label(input_frame2, text="Min Wavelength [nm]:")
numin_label.grid(row=2, column=0)
numin = ttk.Entry(input_frame2)
numin.grid(row=2, column=1)

# Wave step [cm-1] input
wave_step_label = ttk.Label(input_frame2, text="Wave step [cm-1]:")
wave_step_label.grid(row=3, column=0)
wave_step = ttk.Entry(input_frame2)
wave_step.grid(row=3, column=1)

# Molecule of interest multipe choice, considering dictionary of combustion products
molecule_id_label = ttk.Label(input_frame2, text="Molecule of interest:")
molecule_id_label.grid(row=4, column=0)
molecule_id_var = tk.StringVar()
molecule_id_dropdown = ttk.Combobox(input_frame2, textvariable=molecule_id_var, values=["H2O", "CO2", "CO", "N2", "O2","CH4","H2","NO","NO2"]) 
molecule_id_dropdown.grid(row=4, column=1)
molecule_id_dropdown.set("H2O")  # Set a default value

# Method of Compute (New variable)
method_label = ttk.Label(input_frame2, text="Method of Compute HT, Voigt, Lorentz or Doppler:")
method_label.grid(row=5, column=0)
method_var = tk.StringVar()
method_dropdown = ttk.Combobox(input_frame2, textvariable=method_var, values=["HT", "V", "L","D"])  # Add more values as needed
method_dropdown.grid(row=5, column=1)
method_dropdown.set("HT")  # Set a default value

# Button for computing the spectrum
compute_spectrum_button = ttk.Button(input_frame2, text="Compute Spectrum", command=update_spectrum)
compute_spectrum_button.grid(row=6, column=0, columnspan=2)

# Button for clear the plots
clear_spectrum_button = ttk.Button(input_frame2, text="Clear the Spectrum", command=clear_spectrum)
clear_spectrum_button.grid(row=7, column=0, columnspan=2)


# Creating the plots for TAB 2
######### PLOTS #########
fig2, ax3 = plt.subplots(1, 1, figsize=(8,6 )) 
fig2.subplots_adjust(hspace=1.0, wspace=0.5)  # Adjust vertical and horizontal spacing between subplots

# Creating the canvas for embedding the plot
canvas2 = FigureCanvasTkAgg(fig2, master=tab2)  # Display the plot in tab1
canvas_widget2 = canvas2.get_tk_widget()
canvas_widget2.grid(row=7, column=0, columnspan=2)

# Create a frame for the Matplotlib toolbar
toolbar_frame2 = ttk.Frame(tab2)
toolbar_frame2.grid(row=8, column=0, columnspan=2, sticky="nsew")

# Add a Matplotlib toolbar for zooming
toolbar2 = NavigationToolbar2Tk(canvas2, toolbar_frame2)
toolbar2.update()
toolbar2.pack(side=tk.BOTTOM, fill=tk.X)
canvas_widget2.grid(row=9, column=0, columnspan=2, sticky=tk.NSEW)

# Configure grid row and column of plots weights to make them expand with the window
tab2.grid_rowconfigure(9, weight=1)
tab2.grid_columnconfigure(0, weight=1)

# Configure grid row and column of zoombar weights to make them expand with the window
tab2.grid_rowconfigure(9, weight=1)
tab2.grid_columnconfigure(0, weight=1)

#  Display the canvas
canvas2.draw()


# Start the Tkinter main loop
root.mainloop()

