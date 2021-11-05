from typing import *

class Codon:

    def __init__(self, fullname : str, shortname : str, symbol : str, hydropathy : str, chemical_class : str, charge : str, hydrogen_donor_status : str, polarity : str) -> None:
        self.fullname = fullname
        self.shortname = shortname
        self.symbol = symbol
        self.hydropathy = hydropathy
        self.chemical_class = chemical_class
        self.charge = charge
        self.hydrogen_donor_status = hydrogen_donor_status
        self.polarity = polarity

    def translate(self) -> str:
        return self.symbol
    
    def translate_fullname(self) -> str:
        return self.fullname
    
    def translate_shortname(self) -> str:
        return self.shortname
    
    def get_charge(self) -> str:
        return self.charge
    
    def get_hydropathy(self) -> str:
        return self.hydropathy
    
    def get_chemical_class(self) -> str:
        return self.chemical_class

    def get_donor_status(self) -> str:
        return self.hydrogen_donor_status
    
    def get_polarity(self) -> str:
        return self.polarity
    