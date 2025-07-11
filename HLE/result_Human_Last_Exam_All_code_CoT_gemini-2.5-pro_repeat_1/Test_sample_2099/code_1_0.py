import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops

def get_connected_mol_from_smiles(smiles_string):
    """
    Creates a connected RDKit Mol object from a SMILES string.
    If the SMILES represents disconnected fragments (like an ionic salt),
    it connects them with a single bond to form a single graph.
    This is a necessary heuristic for graph-based calculations.
    """
    mol = Chem.MolFromSmiles(smiles_string, sanitize=False)
    # Basic sanitization is needed for property calculation
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_FINALIZE | Chem.SanitizeFlags.SANITIZE_KEKULIZE)

    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        # It's a disconnected graph; we need to connect the pieces.
        # We'll use an editable molecule to add all atoms and bonds,
        # plus new bonds to connect the fragments.
        rw_mol = Chem.RWMol(fragments[0])
        # Connect subsequent fragments to the first atom of the first fragment
        connection_point_in_main_frag = 0

        for i in range(1, len(fragments)):
            frag = fragments[i]
            # Map old atom indices from the fragment to new indices in the combined molecule
            frag_atom_map = {}
            for atom in frag.GetAtoms():
                new_idx = rw_mol.AddAtom(atom)
                frag_atom_map[atom.GetIdx()] = new_idx
            
            # Copy bonds from the fragment
            for bond in frag.GetBonds():
                begin_idx = frag_atom_map[bond.GetBeginAtomIdx()]
                end_idx = frag_atom_map[bond.GetEndAtomIdx()]
                rw_mol.AddBond(begin_idx, end_idx, bond.GetBondType())
            
            # Add a single bond to connect the new fragment to the main one
            first_atom_in_new_frag = frag_atom_map[0]
            rw_mol.AddBond(connection_point_in_main_frag, first_atom_in_new_frag, Chem.BondType.SINGLE)
        
        mol = rw_mol.GetMol()
        Chem.SanitizeMol(mol)

    # Add explicit hydrogens, which are crucial for the mass-based calculations
    return Chem.AddHs(mol)

def solve_chemoinformatics_puzzle():
    """
    Solves the entire problem by calculating energies and Moran's I values.
    """
    # Based on the decoded riddle, these are the molecules to analyze.
    # Y1, Y2, Y3 are all Phosgene. The others are Y4 to Y10.
    smiles_map = {
        'Phosgene': 'O=C(Cl)Cl',
        'Barium sulfate': '[Ba+2].[O-]S(=O)(=O)[O-]',
        'Hydrogen cyanide': 'N#C',
        'Ethyl pyruvate': 'CCOC(=O)C(=O)C',
        'Chlorpyrifos': 'CCOP(=S)(OCC)OC1=NC(Cl)=C(Cl)C=C1Cl',
        'Tetramethylphosphonium chloride': 'C[P+](C)(C)C.[Cl-]',
        'Formamide': 'C(=O)N',
        'Tetraethyl pyrophosphate': 'CCOP(=O)(OCC)OP(=O)(OCC)OCC'
    }

    energies = {}

    print("Step 1: Calculating Mass-Weighted Barysz Graph Energy for each molecule...")
    for name, smiles in smiles_map.items():
        try:
            mol = get_connected_mol_from_smiles(smiles)
            atoms = mol.GetAtoms()
            n_atoms = len(atoms)
            masses = np.array([atom.GetMass() for atom in atoms])
            adj_matrix = rdmolops.GetAdjacencyMatrix(mol)

            # Build the Mass-Weighted Adjacency Matrix
            mass_weighted_adj = np.zeros_like(adj_matrix, dtype=float)
            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    if adj_matrix[i, j] == 1:
                        # Use bond order 1 for simplicity as per the definition
                        val = 1.0 / np.sqrt(masses[i] * masses[j])
                        mass_weighted_adj[i, j] = val
                        mass_weighted_adj[j, i] = val

            # Calculate the energy as the sum of absolute eigenvalues
            eigenvalues = np.linalg.eigvalsh(mass_weighted_adj)
            barysz_energy = np.sum(np.abs(eigenvalues))
            energies[name] = barysz_energy
            print(f"  - {name}: {barysz_energy:.4f}")

        except Exception as e:
            print(f"An error occurred while processing {name}: {e}")
            energies[name] = float('inf')

    # Find the molecule with the minimum energy
    min_energy_molecule_name = min(energies, key=energies.get)
    min_energy_value = energies[min_energy_molecule_name]
    
    print(f"\nStep 2: Identifying the molecule with the lowest energy.")
    print(f"The molecule with the lowest energy is '{min_energy_molecule_name}'.")
    print(f"Identified Energy (E_min): {min_energy_value}")

    # Calculate Moran's I for the identified molecule
    print(f"\nStep 3: Calculating Mass-Weighted Moran's I for '{min_energy_molecule_name}'.")
    
    mol = get_connected_mol_from_smiles(smiles_map[min_energy_molecule_name])
    atoms = mol.GetAtoms()
    n_atoms = len(atoms)
    masses = np.array([atom.GetMass() for atom in atoms])
    dist_matrix = rdmolops.GetDistanceMatrix(mol)
    
    mean_mass = np.mean(masses)
    mass_deviations = masses - mean_mass
    sum_sq_dev = np.sum(mass_deviations**2)

    moran_I_values = []
    if sum_sq_dev > 1e-9: # Avoid division by zero if all masses are the same
        max_dist = int(np.max(dist_matrix))
        for d in range(1, max_dist + 1):
            numerator_sum = 0
            pair_count = 0
            # Sum over all ordered pairs (i, j) at distance d
            indices = np.where(dist_matrix == d)
            pair_count = len(indices[0])
            if pair_count > 0:
                numerator_sum = np.sum(mass_deviations[indices[0]] * mass_deviations[indices[1]])
                # Moran's I formula: I = (N / S) * (sum of cross-products) / (sum of squares)
                # S is the total number of pairs at distance d, which is pair_count
                moran_I = (n_atoms / pair_count) * (numerator_sum / sum_sq_dev)
                moran_I_values.append(moran_I)
                print(f"  - Moran's I for distance d={d}: {moran_I:.4f}")

    if not moran_I_values:
         print("Could not calculate any Moran's I values.")
         min_moran_I = 0.0
         max_moran_I = 0.0
    else:
        min_moran_I = min(moran_I_values)
        max_moran_I = max(moran_I_values)

    print(f"Minimum Moran's I (I_min): {min_moran_I}")
    print(f"Maximum Moran's I (I_max): {max_moran_I}")

    # Final calculation
    final_product = min_energy_value * min_moran_I * max_moran_I

    print("\nStep 4: Final Calculation.")
    print("Result = E_min * I_min * I_max")
    print(f"Product = {min_energy_value} * {min_moran_I} * {max_moran_I}")
    print(f"Final Answer = {final_product}")
    print(f"<<<{final_product}>>>")


if __name__ == '__main__':
    solve_chemoinformatics_puzzle()