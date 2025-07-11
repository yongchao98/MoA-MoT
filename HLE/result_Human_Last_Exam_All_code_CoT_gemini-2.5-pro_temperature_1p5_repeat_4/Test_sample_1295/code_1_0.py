import sys

def solve_molecular_puzzle():
    """
    This function constructs and prints the SMILES representation of a molecule
    based on a series of complex chemical constraints.

    The final deduced molecule is 1,6-bis(dimethylamino)-3,4-diazahexane-2,5-diimine.
    This structure is consistent with the key constraints provided:
    - Molecular Formula: C8H18N6
    - Valence Electrons: (8*4) + (18*1) + (6*5) = 32 + 18 + 30 = 80
    - Molecular Weight: ~198.159 g/mol
    - Heavy Atoms: 14 (8 Carbon, 6 Nitrogen)
    - Heteroatoms: 6 (all Nitrogen)
    - Functional Groups: 2 amidine groups, 1 azo group
    - Nitrogen Types: 2 primary, 2 secondary, 2 tertiary (by connectivity)
    - Acyclic: Yes
    - Formal Charge: 0
    - Rotatable Bonds: 4 (assuming amidine C-N bond rigidity)
    """

    # SMILES representation for (CH3)2N-C(=NH)-CH2-N=N-CH2-C(=NH)-N(CH3)2
    # Broken down:
    # CN(C)   - a dimethylamino group (tert-N)
    # C(=N)   - an imine-amidine part (prim-N)
    # C       - a methylene linker
    # N=N     - the azo group (sec-N x2)
    # C       - a second methylene linker
    # C(=N)   - a second imine-amidine part (prim-N)
    # N(C)C   - a second dimethylamino group (tert-N)
    # The SMILES parser will add implicit hydrogens.
    smiles_string = "CN(C)C(=N)CN=NCC(=N)N(C)C"

    # The problem asks to output numbers from the "final equation".
    # This likely refers to showing the work for the key quantitative values.
    # Note: Python floating point math might have small precision differences.
    print("Derived Molecular Formula: C8H18N6")
    print("Valence Electrons Calculation: (8 * 4) + (18 * 1) + (6 * 5) = 32 + 18 + 30 = 80")
    print("Approximate MW Calculation: (8 * 12.01) + (18 * 1.008) + (6 * 14.007) = 96.08 + 18.144 + 84.042 = 198.266")
    print("\n--- Constructed SMILES Representation ---")
    print(smiles_string)

solve_molecular_puzzle()
<<<CN(C)C(=N)CN=NCC(=N)N(C)C>>>