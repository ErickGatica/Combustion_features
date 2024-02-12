import numpy as np

# Defining dictionary of variables

Compositions_standar={
        'T':330,
        'P':1.2,
        'XMethane':0.899,
        'XEthene':0.0,
        'XEthane': 0.0726,
        'XPropane':0.0044,
        'XButhane': 0.0005,
        'XHydrogen':0.0,
        'XCO':0.0,
        'XCO2':0.0091,
        'XN2':1-0.899-0.0726-0.0044-0.0005-0.0091,
        'Perc':1.0
    }


# Defining the molecular weight of the components
MW={
    'Methane':16.04,
    'Ethane':30.07,
    'Ethene':28.05,
    'Propane':44.09,
    'Buthane':58.12,
    'Hydrogen':2.02,
    'CO':28.01,
    'CO2':44.01,
    'N2':28.01,
    'Oxygen':32.0
}

def combustion(Compositions):
    T=Compositions['T']
    P=Compositions['P']
    XMethane=Compositions['Methane']
    XEthane=Compositions['Ethane']
    XEthene=Compositions['Ethene']
    XPropane=Compositions['Propane']
    XButhane=Compositions['Buthane']
    XHydrogen=Compositions['Hydrogen']
    XCO=Compositions['CO']
    XCO2=Compositions['CO2']
    XN2=Compositions['N2']
    Perc=Compositions['Perc'] # Percentage of the fuel that is burned
    H2ble=Compositions['H2ble'] # Hydrogen bleed in percentage

    Perc=Perc # Selectivity for the combustion, how many percent of the fuel is burned completely
    O2excess=0.00 # 0% of excess of oxygen
    Oxy_fraction=0.21 # Fraction of oxygen in the air

    # Defining the flow of the fuel 
    Baseline=100
    Fuel_flow=Baseline*(1-H2ble/100)
    Hydrogen_flow=Baseline*(H2ble/100)+Fuel_flow*XHydrogen
    # Stoichoimetric equation for combustions
    # Complete combustion generating CO2 and H2O
        # CxHy + (x+y/4)O2 -> xCO2 + y/2H2O
    # Incomplete combusstion generating CO and H2O
        # CxHy + (x+y/2)/2 O2 -> xCO + y/2H2O 
    # Combustion Hydrogen
        # 2H2 + O2 -> 2H2O
    # Production of NOx NO and NO2
        # N2 + O2 -> 2NO
        # N2 + 2O2 -> 2NO2


    # Computing how much oxygen would be necessary for the complete combustion
    Oxyg_methane=XMethane*Fuel_flow*2
    Oxyg_ethane=XEthane*Fuel_flow*7/2
    Oxyg_ethene=XEthene*Fuel_flow*3
    Oxyg_propane=XPropane*Fuel_flow*5
    Oxyg_buthane=XButhane*Fuel_flow*13/2

    # Computing the oxygen for hydrogen combustion
    Oxyg_hydrogen=Hydrogen_flow/2

    # Computin the oxygen for the production of NOx
    # Factor of NO production according to fuel
    NO_fa=1E-1/28*16 #from engineeringtoolbox.com 1E-3
    NO2_fa=1E-1/54*16 #from engineeringtoolbox.com 1E-3
    # Computing the oxygen for the production of NOx
    NO_perc_prod=0.9 # 90% of the nitrogen is converted to NO and 10% to NO2
    Oxyg_NO=NO_fa*Fuel_flow*XMethane*NO_perc_prod
    Oxyg_NO2=2*(NO2_fa*Fuel_flow*XMethane*NO_perc_prod)
    
    # Computing airflow necessary for the combustion
    total_oxygen=Oxyg_methane+Oxyg_ethane+Oxyg_ethene+Oxyg_propane+Oxyg_buthane+Oxyg_NO+Oxyg_NO2 #total oxygen for hydrocarb combustion and NOx production
    O2_necessary=total_oxygen*(1+O2excess)+Oxyg_hydrogen
    Air_flow=O2_necessary/Oxy_fraction # Air is 21% oxygen and 79% nitrogen
    # Computing N2 from the air
    N2_air=Air_flow*(1-Oxy_fraction)

    # Computing the total amount of oxygen necessary for incomplete combustion
    Oxyg_methane_I=XMethane*Fuel_flow*3/2
    Oxyg_ethane_I=XEthane*Fuel_flow*5/2
    Oxyg_ethene_I=XEthene*Fuel_flow*2
    Oxyg_propane_I=XPropane*Fuel_flow*7/2
    Oxyg_buthane_I=XButhane*Fuel_flow*9/2
    total_oxygen_I=Oxyg_methane_I+Oxyg_ethane_I+Oxyg_ethene_I+Oxyg_propane_I+Oxyg_buthane_I

    # Computing the amount of nitrogen converted in NOx
    N_convert=Oxyg_NO+1/2*Oxyg_NO2
    # Computing the new amount of nitrogen in the air
    N2_postNO=N2_air-N_convert
   
    # Computing CO2 generated from complete combustion
    CO2_methane=XMethane*Fuel_flow*Perc/100
    CO2_ethane=2*XEthane*Fuel_flow*Perc/100
    CO2_ethene=2*XEthene*Fuel_flow*Perc/100
    CO2_propane=3*XPropane*Fuel_flow*Perc/100
    CO2_buthane=4*XButhane*Fuel_flow*Perc/100

    # Computing CO generated from incomplete combustion
    CO_methane=XMethane*Fuel_flow*(1-Perc/100)
    CO_ethane=2*XEthane*Fuel_flow*(1-Perc/100)
    CO_ethene=2*XEthene*Fuel_flow*(1-Perc/100)
    CO_propane=3*XPropane*Fuel_flow*(1-Perc/100)
    CO_buthane=4*XButhane*Fuel_flow*(1-Perc/100)

    # Computing the amount of water generated from the combustion
    H2O_methane=2*XMethane*Fuel_flow*Perc/100
    H2O_ethane=3*XEthane*Fuel_flow*Perc/100
    H2O_ethene=2*XEthene*Fuel_flow*Perc/100
    H2O_propane=4*XPropane*Fuel_flow*Perc/100
    H2O_buthane=5*XButhane*Fuel_flow*Perc/100

    # Computing the amount of water generated from the incomplete combustion
    H2O_methane_I=2*XMethane*Fuel_flow*(1-Perc/100)
    H2O_ethane_I=3*XEthane*Fuel_flow*(1-Perc/100)
    H2O_ethene_I=2*XEthene*Fuel_flow*(1-Perc/100)
    H2O_propane_I=4*XPropane*Fuel_flow*(1-Perc/100)
    H2O_buthane_I=5*XButhane*Fuel_flow*(1-Perc/100)

    # Computing the amount of water generated from the combustion of hydrogen
    H2O_hydrogen=Hydrogen_flow*Perc/100

    # Computing all flows
    Oxyg_consum=total_oxygen*Perc/100+total_oxygen_I*(1-Perc/100)+Oxyg_hydrogen
    CO2_prod=CO2_methane+CO2_ethane+CO2_ethene+CO2_propane+CO2_buthane
    CO_prod=CO_methane+CO_ethane+CO_ethene+CO_propane+CO_buthane
    H2O_prod=H2O_methane+H2O_ethane+H2O_ethene+H2O_propane+H2O_buthane+H2O_methane_I+H2O_ethane_I+H2O_ethene_I+H2O_propane_I+H2O_buthane_I+H2O_hydrogen

    # Computing mols flow in the outlet of engine
    #Fuels
    Methane_outlet=XMethane*Fuel_flow*(1-Perc/100)
    Ethane_outlet=XEthane*Fuel_flow*(1-Perc/100)
    Ethene_outlet=XEthene*Fuel_flow*(1-Perc/100)
    Propane_outlet=XPropane*Fuel_flow*(1-Perc/100)
    Buthane_outlet=XButhane*Fuel_flow*(1-Perc/100)
    #Oxygen
    Oxyg_outlet=Air_flow*Oxy_fraction-Oxyg_consum
    #Nitrogen
    N2_outlet=Fuel_flow*XN2 + N2_postNO
    #CO
    CO_outlet=Fuel_flow*XCO+CO_prod
    #CO2
    CO2_outlet=Fuel_flow*XCO2+CO2_prod
    #H2O
    H2O_outlet=H2O_prod
    #H2
    H2_outlet=Hydrogen_flow*(1-Perc/100)
    #NO
    NO_outlet=2*Oxyg_NO
    #NO2
    NO2_outlet=Oxyg_NO2

    # Total mole flow in the outlet
    Total_outlet=Methane_outlet+Ethane_outlet+Ethene_outlet+Propane_outlet+Buthane_outlet+Oxyg_outlet+N2_outlet+CO_outlet+CO2_outlet+H2O_outlet+H2_outlet+NO_outlet+NO2_outlet
    # Computing molar fractions
    XMethane_outlet=Methane_outlet/Total_outlet
    XEthane_outlet=Ethane_outlet/Total_outlet
    XEthene_outlet=Ethene_outlet/Total_outlet
    XPropane_outlet=Propane_outlet/Total_outlet
    XButhane_outlet=Buthane_outlet/Total_outlet
    XOxyg_outlet=Oxyg_outlet/Total_outlet
    XN2_outlet=N2_outlet/Total_outlet
    XCO_outlet=CO_outlet/Total_outlet
    XCO2_outlet=CO2_outlet/Total_outlet
    XH2O_outlet=H2O_outlet/Total_outlet
    XH2_outlet=H2_outlet/Total_outlet
    XNO_outlet=NO_outlet/Total_outlet
    XNO2_outlet=NO2_outlet/Total_outlet
    # Saving results in dictionary
    Results={

        # Outlet products
        'X2_Methane':XMethane_outlet,
        'X2_Ethane':XEthane_outlet,
        'X2_Ethene':XEthene_outlet,
        'X2_Propane':XPropane_outlet,
        'X2_Buthane':XButhane_outlet,
        'X2_O2':XOxyg_outlet,
        'X2_N2':XN2_outlet,
        'X2_CO':XCO_outlet,
        'X2_CO2':XCO2_outlet,
        'X2_H2O':XH2O_outlet,
        'X2_H2':XH2_outlet,
        'X2_NO':XNO_outlet,
        'X2_NO2':XNO2_outlet
    }
    return Results

