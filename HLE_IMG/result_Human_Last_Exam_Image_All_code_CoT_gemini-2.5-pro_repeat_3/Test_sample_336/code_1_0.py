# To run this script, you need to install rdkit and numpy.
# You can install them using pip:
# pip install rdkit numpy

import collections
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def solve_task():
    """
    This function performs the entire task as described in the plan.
    1. Identifies molecules from their SMILES strings.
    2. Filters them based on Crippen logP > 1.
    3. Calculates the Eccentric Connectivity Index (ECI) for the filtered molecules.
    4. Sums the ECI values and prints the result.
    """
    
    # Step 1: Define the molecules from the image using their SMILES strings.
    smiles_map = {
        "Reactant 1 (M1)": "COC1=C(C=C(C(=C1)OC)COC2=C(C=NC(=C2)N)F)OC",
        "Reactant 2 (M2)": "C1N(C2=C(C=NC=N2)C1)Cl",
        "Product (M3)": "C1N(C2=C(C=NC=N2)C1)NC1=NC=CC(F)=C1"
    }
    print("--- Step 1: Identifying Molecules ---")
    for name, smi in smiles_map.items():
        print(f"{name}: {smi}")
    print("-" * 35)

    # Step 2: Filter molecules with Crippen logP > 1.
    print("\n--- Step 2: Filtering by Crippen logP > 1 ---")
    molecules_to_process = []
    for name, smi in smiles_map.items():
        mol = Chem.MolFromSmiles(smi)
        logp = Descriptors.MolLogP(mol)
        print(f"logP for {name}: {logp:.4f}")
        if logp > 1:
            molecules_to_process.append((name, mol))
            print(f"-> Selected '{name}' for calculation.")
        else:
            print(f"-> Excluded '{name}'.")
    print("-" * 35)

    # Step 3: Calculate ECI for each selected molecule.
    print("\n--- Step 3: Calculating Eccentric Connectivity Index (ECI) ---")
    total_eci = 0
    eci_values = []

    for name, mol in molecules_to_process:
        # Add explicit hydrogens to the molecular graph
        mol_with_hs = Chem.AddHs(mol)
        num_atoms = mol_with_hs.GetNumAtoms()

        # Get adjacency matrix and list for graph traversal
        adj_matrix = Chem.GetAdjacencyMatrix(mol_with_hs)
        adj_list = collections.defaultdict(list)
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                if adj_matrix[i, j] == 1:
                    adj_list[i].append(j)
                    adj_list[j].append(i)
        
        # Calculate all-pairs shortest path (distance matrix) using BFS from each node
        dist_matrix = np.full((num_atoms, num_atoms), -1, dtype=int)
        for start_node in range(num_atoms):
            queue = collections.deque([(start_node, 0)])
            visited = {start_node}
            dist_matrix[start_node, start_node] = 0
            while queue:
                u, d = queue.popleft()
                for v in adj_list[u]:
                    if v not in visited:
                        visited.add(v)
                        dist_matrix[start_node, v] = d + 1
                        queue.append((v, d + 1))
        
        # Calculate degrees and eccentricities
        degrees = np.sum(adj_matrix, axis=1)
        eccentricities = np.max(dist_matrix, axis=1)
        
        # Calculate ECI for the molecule
        eci = int(np.sum(degrees * eccentricities))
        eci_values.append(eci)
        total_eci += eci
        
        print(f"\nCalculation for '{name}':")
        print(f"  - Number of atoms (including H): {num_atoms}")
        print(f"  - Eccentric Connectivity Index (ECI): {eci}")
        
    print("-" * 35)

    # Step 4: Sum the indices and display the final result.
    print("\n--- Step 4: Final Summation ---")
    equation_str = " + ".join(map(str, eci_values))
    print(f"The total Sum of Eccentric Connectivity Indices is:")
    print(f"{equation_str} = {int(total_eci)}")

# Execute the main function
solve_task()