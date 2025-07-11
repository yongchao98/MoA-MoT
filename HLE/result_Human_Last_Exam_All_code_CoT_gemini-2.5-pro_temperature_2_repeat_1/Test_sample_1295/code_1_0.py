import sys

def solve_molecular_puzzle():
    """
    This function constructs and provides the SMILES string for a molecule
    based on a complex set of constraints.

    The final molecule is N,N''-(diazene-1,2-diylbis(methylene))bis(N'-ethylcarboximidamide).
    SMILES: CCNC(=N)CN=NCC(=N)NCC

    Let's verify its properties based on the prompt:
    - Total Valence Electrons (80): C8H18N6 -> (8*4 + 18*1 + 6*5) = 32 + 18 + 30 = 80. (Correct)
    - Formal Charge (0): The molecule is neutral. (Correct)
    - Molecular Weight (~198.159): C8H18N6 gives ~198.27 g/mol. (Correct)
    - Heavy Atoms (14): 8 Carbons + 6 Nitrogens = 14. (Correct)
    - Heteroatoms (6): 6 Nitrogens. (Correct)
    - Total NH or OH groups (6): There are 2 C=NH and 2 C-NH groups. This totals 4 groups, not 6. This is a primary conflict in the prompt I had to resolve.
    - Hydrogen Bond Acceptors (4): All 6 nitrogens can act as acceptors. This is another conflict I've noted.
    - Hydrogen Bond Donors (4): The molecule has 4 N-H bonds (2 from C=NH, 2 from C-NH). (Correct)
    - Amine Types: The prompt's requirement for 2 of each (primary, secondary, tertiary) conflicts with the functional group and H-bond donor counts. My molecule contains secondary amine-like nitrogens within the amidine groups.
    - Amidine Groups (2): The structure contains two C(=N)N groups. (Correct)
    - Azo Group (1): The structure contains the central -N=N- azo group. (Correct)
    - No Rings: The molecule is acyclic. (Correct)
    - Rotatable Bonds (4): This specific structure has more than 4 rotatable bonds. It is extremely difficult to satisfy this constraint simultaneously with all others.
    """

    # The SMILES representation for the constructed molecule.
    # The molecule is C_H3-C_H2-N_H-C(=N_H)-C_H2-N=N-C_H2-C(=N_H)-N_H-C_H2-C_H3
    smiles_string = "CCNC(=N)CN=NCC(=N)NCC"

    print("Generated SMILES string based on the provided constraints:")
    print(smiles_string)

solve_molecular_puzzle()