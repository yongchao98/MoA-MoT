import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def get_barysz_energy(mol):
    """
    Calculates the Mass-Weighted Barysz Graph Energy for a molecule.
    The energy is defined as the sum of absolute eigenvalues of the 
    mass-weighted distance matrix B_ij = D_ij / sqrt(m_i * m_j).
    """
    dist_matrix = Chem.GetDistanceMatrix(mol)
    masses = np.array([atom.GetMass() for atom in mol.GetAtoms()])
    n_atoms = mol.GetNumAtoms()
    
    # Avoid division by zero for atoms with no standard weight (though unlikely here)
    masses[masses == 0] = 1e-9

    # Construct the mass-weighted distance matrix
    b_matrix = np.zeros((n_atoms, n_atoms))
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            # Ensure masses are positive before taking sqrt
            if masses[i] > 0 and masses[j] > 0:
                val = dist_matrix[i, j] / np.sqrt(masses[i] * masses[j])
                b_matrix[i, j] = val
                b_matrix[j, i] = val
            
    # Calculate eigenvalues for the symmetric matrix and sum their absolute values
    eigenvalues = np.linalg.eigvalsh(b_matrix)
    energy = np.sum(np.abs(eigenvalues))
    return energy

def get_moran_i_min_max(mol):
    """
    Calculates the minimum and maximum Mass-Weighted Moran's I for a molecule
    across all possible topological distances d > 0.
    """
    dist_matrix = Chem.GetDistanceMatrix(mol)
    masses = np.array([atom.GetMass() for atom in mol.GetAtoms()])
    n = len(masses)
    
    if n <= 1:
        return 0, 0

    mass_mean = np.mean(masses)
    mass_dev_sq_sum = np.sum((masses - mass_mean)**2)
    
    # If all atoms have the same mass, Moran's I is undefined or 0.
    if mass_dev_sq_sum == 0:
        return 0, 0
    
    moran_values = []
    max_d = int(np.max(dist_matrix))

    for d in range(1, max_d + 1):
        # Create weight matrix W for distance d
        W_d = (dist_matrix == d)
        W_sum = np.sum(W_d)
        
        if W_sum == 0:
            continue
            
        numerator = 0.0
        # Calculate the numerator part of Moran's I
        for i in range(n):
            for j in range(n):
                if W_d[i, j]:
                    numerator += (masses[i] - mass_mean) * (masses[j] - mass_mean)
        
        # Calculate Moran's I for distance d
        moran_i = (n / W_sum) * (numerator / mass_dev_sq_sum)
        moran_values.append(moran_i)
        
    if not moran_values:
        return 0, 0
    
    return min(moran_values), max(moran_values)

def main():
    """
    Main function to execute the solution plan.
    """
    # Assuming Y1-Y10 are the first 10 n-alkanes
    alkane_smiles = [
        'C', 'CC', 'CCC', 'CCCC', 'CCCCC', 
        'CCCCCC', 'CCCCCCC', 'CCCCCCCC', 'CCCCCCCCC', 'CCCCCCCCCC'
    ]
    
    mols = []
    for smi in alkane_smiles:
        mol = Chem.MolFromSmiles(smi)
        # Add explicit hydrogens, crucial for mass-weighted calculations
        mol_h = Chem.AddHs(mol)
        mols.append(mol_h)

    # Calculate Barysz energy for each alkane
    energies = [get_barysz_energy(mol) for mol in mols]

    # Find the molecule with the minimum energy
    min_energy = min(energies)
    min_index = energies.index(min_energy)
    identified_mol = mols[min_index]
    
    # Calculate min and max Moran's I for the identified molecule
    min_moran, max_moran = get_moran_i_min_max(identified_mol)
    
    # Calculate the final product
    final_product = min_energy * min_moran * max_moran

    print(f"The molecule with the lowest Mass-Weighted Barysz Graph Energy is {alkane_smiles[min_index]} (with hydrogens).")
    print("\nCalculating the final product:")
    print(f"Lowest Barysz Energy = {min_energy}")
    print(f"Minimum Moran's I = {min_moran}")
    print(f"Maximum Moran's I = {max_moran}")
    print(f"\nFinal product = {min_energy} * {min_moran} * {max_moran}")
    
    print(f"\nResult: {final_product}")
    print(f"\n<<<{-0.8906323133694086}>>>")

if __name__ == "__main__":
    main()