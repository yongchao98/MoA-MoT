import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen, rdmolops, Descriptors
import numpy as np
from collections import deque

def calculate_logp(smiles_dict):
    """Calculates Crippen logP for a dictionary of molecules."""
    qualifying_mols = []
    print("Step 1: Calculate Crippen logP for each depicted molecule.")
    for name, smi in smiles_dict.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            formula = Descriptors.CalcMolFormula(mol)
            logp = Crippen.MolLogP(mol)
            print(f"  - logP for {name} ({formula}): {logp:.4f}")
            if logp > 1:
                print(f"    -> Qualifies (logP > 1).")
                qualifying_mols.append((name, smi, formula))
            else:
                print(f"    -> Does not qualify.")
        else:
            print(f"  - Could not process SMILES for {name}.")
    return qualifying_mols

def calculate_eccentric_connectivity_index(smi):
    """Calculates the Eccentric Connectivity Index (including H) for a molecule."""
    mol = Chem.MolFromSmiles(smi)
    mol_with_h = Chem.AddHs(mol)
    num_atoms = mol_with_h.GetNumAtoms()
    adj_matrix = rdmolops.GetAdjacencyMatrix(mol_with_h)

    # Use Breadth-First Search from each node to find all-pairs shortest paths
    dist_matrix = np.zeros((num_atoms, num_atoms), dtype=int)
    for i in range(num_atoms):
        q = deque([(i, 0)])
        visited = {i}
        while q:
            u, d = q.popleft()
            dist_matrix[i, u] = d
            for v in range(num_atoms):
                if adj_matrix[u, v] and v not in visited:
                    visited.add(v)
                    q.append((v, d + 1))
    
    degrees = [atom.GetDegree() for atom in mol_with_h.GetAtoms()]
    eccentricities = [np.max(dist_matrix[i]) for i in range(num_atoms)]
    
    eci = sum(d * e for d, e in zip(degrees, eccentricities))
    
    return eci, degrees, eccentricities

def main():
    # SMILES strings for the three molecules depicted in the image
    smiles_data = {
        'Molecule 1': 'NC1=C(OCC2=CC(OC)=C(OC)C=C2OC)C=C(F)C=N1',
        'Molecule 2': 'ClN1CC2=CN=CN=C2C1',
        'Molecule 3': 'Fc1ccc(NN2CC3=C(C2)N=CN=C3)nc1',
    }
    
    qualifying_molecules = calculate_logp(smiles_data)
    
    if not qualifying_molecules:
        print("\nNo molecule found with Crippen logP > 1.")
        return

    print("\nStep 2: Calculate the Sum of Eccentric Connectivity Indices for qualifying molecule(s).")
    total_eci = 0
    full_equation = []

    for name, smi, formula in qualifying_molecules:
        print(f"\nProcessing {name} ({formula})...")
        eci, degrees, eccentricities = calculate_eccentric_connectivity_index(smi)
        total_eci += eci
        
        print(f"The molecule has {len(degrees)} atoms (including hydrogens).")
        print("The calculation for the Eccentric Connectivity Index is the sum of (degree * eccentricity) for each atom.")
        
        equation_terms = [f"{d}*{e}" for d, e in zip(degrees, eccentricities)]
        full_equation.extend(equation_terms)
        
    print("\nFinal Calculation:")
    print("Sum = " + " + ".join(full_equation))
    print(f"\nTotal Sum of Eccentric Connectivity Indices = {total_eci}")
    
if __name__ == "__main__":
    main()