import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops

def solve_chemical_puzzle():
    """
    Solves the puzzle by assuming the molecules are a series of PAHs,
    calculating the specified properties, and finding the final product.
    """
    
    # Step 1: Assume the identity of the molecules as a series of common PAHs
    # This is an educated guess as the cipher is not provided.
    smiles_list = [
        "c1ccccc1",                                  # Y1: Benzene
        "c1ccc2ccccc2c1",                             # Y2: Naphthalene
        "c1cc2ccccc2cc1",                             # Y3: Anthracene
        "c1ccc2c(c1)ccc3ccccc23",                      # Y4: Phenanthrene
        "c1cc2cccc3ccc(c1)cc23",                      # Y5: Pyrene
        "c1ccc2c(c1)cc3ccc4ccccc4c3c2",                # Y6: Chrysene
        "c1cc2cccc3c2c(c1)c4cccc5cccc(c3)c45",       # Y7: Perylene
        "c1ccc2c(c1)c3ccc4cccc5c4c3c(c2)c5",       # Y8: Benzo[a]pyrene
        "c1cc2cc3cc4ccc5cc6cc1c(c2c3)c4c56",           # Y9: Coronene
        "c1cc2cc3cc4cc5c6c7c(cc8cc9cc1c%10c%11c9c8c(c76)c%12c5c4c3c%12%11)c2%10" # Y10: Ovalene
    ]

    energies = []
    molecules = []
    
    # Step 2: Calculate Mass-Weighted Barysz Graph Energy for each molecule
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Skip any invalid SMILES string
            continue
            
        mol_h = Chem.AddHs(mol)
        molecules.append(mol_h)
        
        num_atoms = mol_h.GetNumAtoms()
        masses = np.array([atom.GetMass() for atom in mol_h.GetAtoms()])
        adj_matrix = rdmolops.GetAdjacencyMatrix(mol_h)
        
        barysz_matrix = np.zeros((num_atoms, num_atoms))
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                if adj_matrix[i, j] == 1:
                    val = (masses[i] * masses[j])**(-0.5)
                    barysz_matrix[i, j] = barysz_matrix[j, i] = val
        
        eigenvalues = np.linalg.eigvalsh(barysz_matrix)
        energy = np.sum(np.abs(eigenvalues))
        energies.append(energy)
        
    # Step 3: Identify the molecule with the lowest energy
    min_energy_index = np.argmin(energies)
    lowest_energy = energies[min_energy_index]
    target_molecule = molecules[min_energy_index]

    # Step 4: Calculate Mass-Weighted Moran's I for the identified molecule
    masses = np.array([atom.GetMass() for atom in target_molecule.GetAtoms()])
    num_atoms = target_molecule.GetNumAtoms()
    mean_mass = np.mean(masses)
    z = masses - mean_mass
    z_sq_sum = np.sum(z**2)
    
    dist_matrix = rdmolops.GetDistanceMatrix(target_molecule)
    max_dist = int(np.max(dist_matrix))
    
    moran_i_values = []
    if z_sq_sum > 0:
        for d in range(1, max_dist + 1):
            wij_zizj_sum = 0.0
            sd = 0
            
            # Sum over all pairs (i, j)
            for i in range(num_atoms):
                for j in range(num_atoms):
                    if i == j: continue
                    if dist_matrix[i, j] == d:
                        wij_zizj_sum += z[i] * z[j]
                        sd += 1
            
            if sd > 0:
                i_d = (num_atoms / sd) * (wij_zizj_sum / z_sq_sum)
                moran_i_values.append(i_d)
                
    min_moran_i = min(moran_i_values) if moran_i_values else 0
    max_moran_i = max(moran_i_values) if moran_i_values else 0

    # Step 5: Final Calculation
    final_product = lowest_energy * min_moran_i * max_moran_i
    
    # Step 6: Output the result
    print("The molecule with the lowest Mass-Weighted Barysz Graph Energy was identified.")
    print("The final result is the product of this energy, and the minimum and maximum values of Mass-Weighted Moran's I for this molecule.")
    print("\n--- Calculation Components ---")
    print(f"Lowest Energy (Y) = {lowest_energy}")
    print(f"Minimum Moran's I = {min_moran_i}")
    print(f"Maximum Moran's I = {max_moran_i}")
    print("\n--- Final Equation ---")
    print(f"{lowest_energy} * {min_moran_i} * {max_moran_i} = {final_product}")

    return final_product

if __name__ == '__main__':
    result = solve_chemical_puzzle()
    print(f"\n<<<{result}>>>")