import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def solve():
    """
    Calculates the Szeged/Wiener index ratio for perylene-3-thiol, including hydrogen atoms.
    """
    # Step 1 & 2: Define and represent the molecule, perylene-3-thiol.
    # The major reduction product of di(perylene-3-yl) disulfide is perylene-3-thiol.
    # SMILES string for perylene-3-thiol.
    smiles = "Sc1cc2ccc3c4ccc5ccc6c(c12)c3c4c56"
    mol = Chem.MolFromSmiles(smiles)

    # Step 3: Add hydrogens to create the full molecular graph as specified.
    mol_h = Chem.AddHs(mol)
    num_atoms = mol_h.GetNumAtoms()

    # Calculate the all-pairs shortest path distance matrix for the graph.
    dist_matrix = Chem.GetDistanceMatrix(mol_h)

    # Step 4: Calculate the Wiener Index (W).
    # The Wiener index is the sum of distances between all unique pairs of vertices,
    # which is half the sum of all elements in the symmetric distance matrix.
    wiener_index = np.sum(dist_matrix) / 2.0

    # Step 5: Calculate the Szeged Index (Sz).
    szeged_index = 0
    # Iterate over every bond (edge) in the molecule.
    for bond in mol_h.GetBonds():
        # Get the indices of the two atoms connected by the bond.
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()

        n_u = 0
        n_v = 0
        # For the current bond (u,v), iterate over all atoms (w) in the molecule
        # to partition them based on their distance to u and v.
        for w in range(num_atoms):
            dist_wu = dist_matrix[w, u]
            dist_wv = dist_matrix[w, v]
            if dist_wu < dist_wv:
                n_u += 1
            elif dist_wv < dist_wu:
                n_v += 1
        
        # Add the product of the partition sizes to the Szeged index.
        szeged_index += n_u * n_v

    # Step 6: Compute the final ratio.
    # Ensure wiener_index is not zero to avoid division errors.
    if wiener_index > 0:
        ratio = szeged_index / wiener_index
    else:
        ratio = 0

    # Print the results in a clear, step-by-step format.
    print(f"Molecule: Perylene-3-thiol")
    print(f"Formula: {rdMolDescriptors.CalcMolFormula(mol_h)}")
    print(f"Total number of atoms (including H): {num_atoms}")
    print("-" * 30)
    print(f"Calculated Wiener Index (W): {int(wiener_index)}")
    print(f"Calculated Szeged Index (Sz): {int(szeged_index)}")
    print("-" * 30)
    print(f"Szeged/Wiener Ratio = {int(szeged_index)} / {int(wiener_index)} = {ratio}")

solve()
<<<11.633484162895928>>>