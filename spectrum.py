from hapi import *
import numpy as np

def spectrum(P,T,length,numin,numax,molecule_id,isotopo_id,method,wavestep,Diluent1):
    # Defining the molecule
    molecule=moleculeName(molecule_id);
    name_isoto=isotopologueName(molecule_id,isotopo_id); #Getting the name of the isotopologue from the id
    
    # Getting the data
    fetch(molecule,molecule_id,isotopo_id,numin,numax) 


# Lets go with the data
    if method=="HT":
        nu,coef=absorptionCoefficient_HT(SourceTables=molecule,WavenumberStep=wavestep,Diluent=Diluent1,Environment={'T':T,'p':P,'l':length})
    elif method=="V":
        nu,coef=absorptionCoefficient_Voigt(SourceTables=molecule,WavenumberStep=wavestep,Diluent=Diluent1,Environment={'T':T,'p':P,'l':length})
    elif method=="L":
        nu,coef=absorptionCoefficient_Lorentz(SourceTables=molecule,WavenumberStep=wavestep,Diluent=Diluent1,Environment={'T':T,'p':P,'l':length})
    elif method=="D":
        nu,coef=absorptionCoefficient_Doppler(SourceTables=molecule,WavenumberStep=wavestep,Diluent=Diluent1,Environment={'T':T,'p':P,'l':length})


    # Lets go with absorption Spectrum
    nu,absorp=absorptionSpectrum(nu,coef)
    # Lets go with transmittance Spectrum
    nu,trans=transmittanceSpectrum(nu,coef)
    # Lets go with radiance Spectrum
    nu,radi = radianceSpectrum(nu,coef)

    # Lets define a class to save the data
    class Data_spectrum:
        def __init__(self, nu, coef, absorp, trans, radi,name_isoto):
            self.nu = nu
            self.coef = coef
            self.absorp = absorp
            self.trans = trans
            self.radi = radi
            self.name_isoto=name_isoto

    return Data_spectrum(nu, coef, absorp, trans, radi, name_isoto)