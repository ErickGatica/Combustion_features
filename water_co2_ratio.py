
# Code to plot the water to carbon dioxide ratio as a function of Hydrogen blending and % of complete combustion for 1.5[atm] and 873[K]
# Importing libraries
import numpy as np
import matplotlib.pyplot as plt

from combustion import combustion

Compositions_standar={
        'T':873.15,#K
        'P':0.94,#atm
        'Methane':0.899,
        'Ethene':0.0,
        'Ethane': 0.0726,
        'Propane':0.0044,
        'Buthane': 0.0005,
        'Hydrogen':0.0,
        'CO':0.0,
        'CO2':0.0091,
        'N2':1-0.899-0.0726-0.0044-0.0005-0.0091,
        'Perc':0.0,
        'H2ble':0.0
    }

# Defining linspace vector of hydrogen blending and % of complete combustion
res=100
Hydrogen_B=np.linspace(0,25,res)
Perc=np.linspace(85,res,res)
lenght=len(Hydrogen_B)
# Defining vector to save the data of water, co2 and co molar fraction
H2O_outlet=np.zeros((lenght,lenght))
CO2_outlet=np.zeros((lenght,lenght))
CO_outlet=np.zeros((lenght,lenght))
Water_CO2_ratio=np.zeros((lenght,lenght))
# Computing the data

for i in range(0,lenght): # Blending loop
    for j in range(0,lenght): # Complete combustion loop
        Compositions_standar['H2ble']=Hydrogen_B[i]
        Compositions_standar['Perc']=Perc[j]
        Compositions=combustion(Compositions_standar)
        H2O_outlet[i,j]=Compositions['X2_H2O']
        CO2_outlet[i,j]=Compositions['X2_CO2']
        CO_outlet[i,j]=Compositions['X2_CO']

        Water_CO2_ratio[i,j]=Compositions['X2_H2O']/Compositions['X2_CO2']

# Plotting the data in a 3D plot of the ratio in function of blending and complete combustion
# Plot for the ratio of H2O/CO2
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(Perc, Hydrogen_B)
ax.plot_surface(X, Y, Water_CO2_ratio, cmap='viridis', edgecolor='none')
ax.set_ylabel('Hydrogen blending [%]')
ax.set_xlabel('Complete combustion [%]')
ax.set_zlabel('H2O/CO2 ratio')


# Plot for the molar fraction of H2O
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, H2O_outlet, cmap='viridis', edgecolor='none')
ax.set_ylabel('Hydrogen blending [%]')
ax.set_xlabel('Complete combustion [%]')
ax.set_zlabel('H2O molar fraction')


# Plot for the molar fraction of CO2
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, CO2_outlet,  cmap='viridis', edgecolor='none')
ax.set_ylabel('Hydrogen blending [%]')
ax.set_xlabel('Complete combustion [%]')
ax.set_zlabel('CO2 molar fraction')

# Plot for the molar fraction of CO
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, CO_outlet, cmap='viridis', edgecolor='none')
ax.set_ylabel('Hydrogen blending [%]')
ax.set_xlabel('Complete combustion [%]')
ax.set_zlabel('CO molar fraction')


# Plotting the CO2 molar in function of the blending for 100 combustion
#fig = plt.figure()
#plt.plot(Hydrogen_B,CO2_outlet[:,99])
#plt.xlabel('Hydrogen blending [%]')
#plt.ylabel('CO2 molar fraction')

plt.show()