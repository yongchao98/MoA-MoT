# First, you may need to install the required libraries. You can do this by running:
# pip install rdkit numpy scipy

import rdkit
from rdkit import Chem
import numpy as np
from scipy.sparse.csgraph import floyd_warshall

def solve():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol.
    """
    # Step 1 & 2: Define the molecule (perylene-3-thiol) from its SMILES string
    # and add explicit hydrogen atoms to create the full molecular graph.
    # The SMILES for Perylene-3-thiol is C1=CC=C2C3=C4C=CC=C5C6=C(C(=C1)C2=C6C=C(C=C4)S)C=C5
    smiles = "C1=CC=C2C3=C4C=CC=C5C6=C(C(=C1)C2=C6C=C(C=C4)S)C=C5"
    mol = Chem.MolFromSmiles(smiles)
    mol_with_hs = Chem.AddHs(mol)
    num_atoms = mol_with_hs.GetNumAtoms()

    # Get the adjacency matrix, which represents the connections (bonds) in the molecule.
    adj_matrix = Chem.GetAdjacencyMatrix(mol_with_hs)

    # Step 3 (Part 1): Calculate the Wiener Index (W)
    # First, compute the all-pairs shortest path matrix (distance matrix) using the Floyd-Warshall algorithm.
    dist_matrix = floyd_warshall(csgraph=adj_matrix, directed=False)

    # The Wiener index is the sum of distances between all unique pairs of atoms.
    # This is equivalent to summing the upper (or lower) triangle of the distance matrix.
    wiener_index = np.sum(np.triu(dist_matrix, k=1))

    # Step 3 (Part 2): Calculate the Szeged Index (Sz)
    szeged_index = 0
    # Iterate over each bond (edge) in the molecule.
    for bond in mol_with_hs.GetBonds():
        # Get the indices of the atoms connected by the bond.
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()

        n_u = 0  # Counter for atoms closer to u
        n_v = 0  # Counter for atoms closer to v

        # For the current bond (u,v), iterate over all atoms (w) in the molecule.
        for w in range(num_atoms):
            dist_w_u = dist_matrix[w, u]
            dist_w_v = dist_matrix[w, v]

            # Partition the atoms based on their distance from u and v.
            if dist_w_u < dist_w_v:
                n_u += 1
            elif dist_w_v < dist_w_u:
                n_v += 1
        
        # Add the contribution of this bond to the total Szeged index.
        szeged_index += n_u * n_v

    # Step 4: Compute the final ratio
    if wiener_index == 0:
        ratio = 0.0
    else:
        ratio = szeged_index / wiener_index

    # Print the final equation with all components
    print("Molecule: Perylene-3-thiol (C20H12S)")
    print(f"Total atoms (vertices), N: {num_atoms}")
    print(f"Calculated Wiener Index (W): {int(wiener_index)}")
    print(f"Calculated Szeged Index (Sz): {int(szeged_index)}")
    print("\n--- Final Ratio Calculation ---")
    print(f"Sz / W = {int(szeged_index)} / {int(wiener_index)} = {ratio}")

if __name__ == '__main__':
    solve()