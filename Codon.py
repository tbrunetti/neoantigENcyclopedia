from typing import *

class Codon:

    '''
    defines all requirements of a Codon object
    5 attributes have a specified number of choices
    '''
    HYDROPATHY = ('hydrophobic', 'hydrophilic', 'neutral', 'n/a')
    CHEMICAL_CLASS = ('aliphatic', 'aromatic', 'hydroxyl', 'sulfur', 'basic', 'amide', 'acidic', 'n/a')
    CHARGE = ('positive', 'negative', 'uncharged', 'n/a')
    DONOR_STATUS = ('none', 'acceptor', 'donor', 'donor and acceptor', 'n/a')
    POLARITY = ('polar', 'nonpolar', 'n/a')
    
    def __init__(self, fullname : str, shortname : str, symbol : str, hydropathy : str, chemical_class : str, charge : str, 
                 hydrogen_donor_status : str, polarity : str, hydropathy_options = HYDROPATHY, chem_options = CHEMICAL_CLASS,
                 charge_options = CHARGE, donor_options = DONOR_STATUS, polarity_options = POLARITY) -> None:
        
        # check attribute types
        if hydropathy not in hydropathy_options:
            raise ValueError('%s is not an option for the hydropathy attribute.' % hydropathy)
        
        if chemical_class not in chem_options:
            raise ValueError('%s is not an option for the chemical_class attribute.' % chemical_class)
        
        if charge not in charge_options:
            raise ValueError('%s is not an option for the charge attribute.' % charge)
        
        if hydrogen_donor_status not in donor_status:
            raise ValueError('%s is not an option for the hydrogen_donor_status attribute.' % hydrogen_donor_status)
        
        if polarity not in polarity_options:
            raise ValueError('%s is not an option for the polarity attribute.' % polarity)
        
        self.fullname = fullname
        self.shortname = shortname
        self.symbol = symbol
        self.hydropathy = hydropathy
        self.chemical_class = chemical_class
        self.charge = charge
        self.hydrogen_donor_status = hydrogen_donor_status
        self.polarity = polarity


    def translate(self) -> str:
        '''
        input: Codon object (self)
        description: method applied to object Codon that returns the single character amino acid symbol
        return: string of length 1
        '''
        return self.symbol
    
    def translate_fullname(self) -> str:
        '''
        input: Codon object (self) 
        description: method applied to Codon object that returns the full official name of amino acid
        return: string
        '''
        return self.fullname
    
    def translate_shortname(self) -> str:
        '''
        input: Codon object (self)
        description: method applied to Codon object that returns the abbreviated amino acid name
        return: string typically of length 3, with exception of STOP codons
        '''
        return self.shortname
    
    def get_charge(self) -> str:
        '''
        input: Codon object (self)
        description: method applied to Codon object to get the charge status of the amino acid generated from codon
        return: string of one of following options: positive, negative, uncharged, n/a
        '''
        return self.charge
    
    def get_hydropathy(self) -> str:
        '''
        input: Codon object (self)
        description: method applied to Codon object to get the hydropathy status of the amino acid generated from codon
        return: string of one of following options: hydrophobic, hydrophilic, neutral, n/a
        '''
        return self.hydropathy
    
    def get_chemical_class(self) -> str:
        '''
        input: Codon object (self)
        description: method applied to Codon object to get the chemical class of the amino acid generated from codon
        return: string of one of following options: aliphatic, aromatic, hydroxyl, sulfur, basic, amide, acidic, n/a
        '''
        return self.chemical_class

    def get_donor_status(self) -> str:
        '''
        input: Codon object (self)
        description: method applied to Codon object to get the hydrogen donor status of the amino acid generated from codon
        return: string of one of following options: none, acceptor, donor, donor and acceptor, n/a
        '''
        return self.hydrogen_donor_status
    
    def get_polarity(self) -> str:
        '''
        input: Codon object (self)
        description: method applied to Codon object to get the polarity status of the amino acid generated from codon
        return: string of one of following options: polar, nonpolar, n/a
        '''
        return self.polarity
    