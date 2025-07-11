import numpy as np
from rdkit import Chem

def calculate_indices_ratio():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol.

    The function performs the following steps:
    1. Defines the molecule perylene-3-thiol using its SMILES string.
    2. Adds explicit hydrogen atoms to the molecule, as required by the problem.
    3. Computes the distance matrix (shortest paths between all atom pairs).
    4. Calculates the Wiener index (W) by summing all pairwise distances.
    5. Calculates the Szeged index (Sz) by iterating over each bond and applying its definition.
    6. Computes the ratio Sz / W and prints the results.
    """
    # SMILES string for perylene-3-thiol, the major reduction product of di(perylene-3-yl) disulfide.
    smiles = 'SC1=CC2=C3C=CC4=CC=CC5=C4C3=C(C=C2)C2=C1C=CC=C25'

    # Create an RDKit molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    # Add explicit hydrogen atoms to the molecule graph
    mol_h = Chem.AddHs(mol)

    # Get the total number of atoms (vertices in the graph)
    num_atoms = mol_h.GetNumAtoms()

    # Calculate the all-pairs shortest path distance matrix
    dist_matrix = Chem.GetDistanceMatrix(mol_h)

    # Calculate the Wiener Index (W)
    # W is half the sum of all elements in the distance matrix.
    wiener_index = np.sum(dist_matrix) / 2.0

    # Calculate the Szeged Index (Sz)
    szeged_index = 0
    # Iterate over all bonds (edges) in the molecule
    for bond in mol_h.GetBonds():
        # Get the indices of the two atoms in the bond
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()

        n_u = 0  # Counter for atoms closer to u
        n_v = 0  # Counter for atoms closer to v

        # For the current bond(u,v), partition all atoms in the graph
        for w in range(num_atoms):
            dist_wu = dist_matrix[w, u]
            dist_wv = dist_matrix[w, v]
            if dist_wu < dist_wv:
                n_u += 1
            elif dist_wv < dist_wu:
                n_v += 1
        
        # Add the product of the partition sizes to the index
        szeged_index += n_u * n_v

    # Calculate the final ratio
    if wiener_index == 0:
        ratio = 0.0 # Avoid division by zero for single-atom molecules
    else:
        ratio = szeged_index / wiener_index

    # Print the results in a clear format, including the numbers for the final equation
    print(f"Molecule: Perylene-3-thiol (C20H11SH)")
    print(f"Number of atoms (including H): {num_atoms}")
    print(f"Szeged Index (Sz): {int(szeged_index)}")
    print(f"Wiener Index (W): {int(wiener_index)}")
    print(f"Szeged/Wiener Ratio = {int(szeged_index)} / {int(wiener_index)} = {ratio}")

if __name__ == '__main__':
    # Execute the calculation
    # Note: RDKit is required. You can install it via pip: pip install rdkit
    try:
        calculate_indices_ratio()
    except ImportError:
        print("RDKit library not found. Please install it using 'pip install rdkit-pypi'")
