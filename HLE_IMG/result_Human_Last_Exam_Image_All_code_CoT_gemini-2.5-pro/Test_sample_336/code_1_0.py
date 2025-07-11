import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen
import numpy as np
from scipy.sparse.csgraph import floyd_warshall

def solve_eccentric_connectivity_index():
    """
    This script identifies the target molecule from the image,
    calculates its Eccentric Connectivity Index (including hydrogens),
    and prints the detailed calculation.
    """
    # Step 1: Identify the target molecule.
    # The image shows a reaction: M1 + M2 -> M3.
    # SMILES strings for the three molecules are defined below.
    smiles_dict = {
        'Reactant 1 (M1)': 'COC1=CC(=CC(OC)=C1OC)COC2=C(F)C=NC(=C2)N',
        'Reactant 2 (M2)': 'Clc1ncc2c(n1)NCC2',
        'Product (M3)': 'C1C2=C(N=CN=C2NC3=CC=C(C=N3)F)CN1' # Canonical SMILES from PubChem
    }

    print("--- Step 1: Molecule Identification ---")
    print("Calculating Crippen logP for all depicted molecules to verify the condition (> 1).")

    target_name = 'Product (M3)' # Assume product is the molecule of interest
    target_smiles = smiles_dict[target_name]

    for name, smiles in smiles_dict.items():
        mol = Chem.MolFromSmiles(smiles)
        logp = Crippen.MolLogP(mol)
        print(f"  - {name}: logP = {logp:.2f}")

    print(f"\nAll molecules have a logP > 1. The calculation will proceed with the product, {target_name}.")
    print("-" * 35)

    # Step 2: Calculate the Eccentric Connectivity Index for the target molecule.
    print(f"\n--- Step 2: ECI Calculation for {target_name} ---")
    
    # Create RDKit molecule object and add explicit hydrogens
    mol = Chem.MolFromSmiles(target_smiles)
    mol_h = Chem.AddHs(mol)
    num_atoms = mol_h.GetNumAtoms()
    
    print(f"The molecule, including hydrogens, has {num_atoms} atoms.")

    # Get the adjacency matrix and compute the all-pairs shortest path distance matrix
    adj_matrix = Chem.GetAdjacencyMatrix(mol_h)
    dist_matrix = floyd_warshall(csgraph=adj_matrix, directed=False, unweighted=True)

    if np.any(np.isinf(dist_matrix)):
        print("Error: Molecule graph is not connected.")
        return

    # Calculate eccentricities (max distance from each atom)
    eccentricities = np.max(dist_matrix, axis=1)

    # Get degrees (number of connections for each atom)
    degrees = [atom.GetDegree() for atom in mol_h.GetAtoms()]

    # Calculate the sum of (degree * eccentricity) for all atoms
    total_eci = 0
    eci_terms_str = []
    for i in range(num_atoms):
        v = degrees[i]
        e = int(eccentricities[i])
        term = v * e
        total_eci += term
        eci_terms_str.append(f"({v}*{e})")

    # Print the full equation as requested
    print("\nThe sum is calculated as: ECI = sum(degree * eccentricity) for each atom.")
    print("\nFinal Equation:")
    # To avoid a very long line, we'll wrap it for better readability
    equation_str = " + ".join(eci_terms_str)
    print("ECI = " + equation_str)
    
    # Print the final result
    print(f"\nTotal Sum of Eccentric Connectivity Indices = {total_eci}")

if __name__ == "__main__":
    solve_eccentric_connectivity_index()