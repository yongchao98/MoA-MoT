import rdkit
from rdkit import Chem
from rdkit.Chem import Draw

def get_compound_name():
    """
    Identifies and returns the name of the compound that reacts with potassium hydroxide
    to produce 1-methyl-4,4a,5,6,7,8-hexahydronaphthalen-2(3H)-one.

    The reaction is an intramolecular aldol condensation, a key step in the
    Robinson annulation. The required compound is the 1,5-diketone precursor.
    """
    # The systematic name of the diketone precursor
    compound_name = "2-(1-methyl-3-oxobutyl)cyclohexan-1-one"
    
    print("The compound that reacted with potassium hydroxide is:")
    print(compound_name)

if __name__ == "__main__":
    get_compound_name()
