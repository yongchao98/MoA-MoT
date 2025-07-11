import numpy as np
from rdkit import Chem
from rdkit.Chem import rdchem

def get_atomic_masses(mol):
    """Returns a list of atomic masses for a molecule, including explicit Hs."""
    mol = Chem.AddHs(mol)
    return [atom.GetMass() for atom in mol.GetAtoms()]

def get_barysz_matrix(mol):
    """Calculates the Mass-Weighted Barysz Matrix."""
    mol = Chem.AddHs(mol)
    n_atoms = mol.GetNumAtoms()
    masses = get_atomic_masses(mol)
    dist_matrix = Chem.GetDistanceMatrix(mol)
    
    c_atom = rdchem.Atom(6)
    m_c = c_atom.GetMass()

    b_matrix = np.zeros((n_atoms, n_atoms))
    for i in range(n_atoms):
        for j in range(n_atoms):
            if i == j:
                b_matrix[i, j] = 1.0 - (m_c / masses[i])
            else:
                if dist_matrix[i, j] == 0:
                    b_matrix[i, j] = 0
                else:
                    # Mass-weighted version of Barysz Matrix
                    numerator = m_c * m_c
                    denominator = dist_matrix[i, j] * masses[i] * masses[j]
                    b_matrix[i, j] = numerator / denominator
    return b_matrix

def get_barysz_energy(b_matrix):
    """Calculates the energy (sum of abs eigenvalues) of a Barysz matrix."""
    eigenvalues = np.linalg.eigvals(b_matrix)
    return np.sum(np.abs(eigenvalues.real))

def get_moran_i_min_max(mol):
    """Calculates the min and max of local Mass-Weighted Moran's I."""
    mol = Chem.AddHs(mol)
    adj_matrix = Chem.GetAdjacencyMatrix(mol)
    masses = np.array(get_atomic_masses(mol))
    n = len(masses)

    if n <= 1:
        return 0, 0

    # Calculate z-scores for masses
    mean_mass = np.mean(masses)
    std_dev = np.std(masses)
    
    if std_dev == 0:
        return 0, 0

    z = (masses - mean_mass) / std_dev

    # Calculate local Moran's I for each atom
    # I_i = z_i * sum_j(w_ij * z_j)
    local_i = z * (adj_matrix @ z)
    
    return np.min(local_i), np.max(local_i)

def solve_challenge():
    """
    Main function to solve the challenge by identifying molecules,
    calculating properties, and finding the final product.
    """
    molecules = {
        "Y1": "Oc1ccc(N=Nc2ccccc2)c2ccccc12",  # Sudan I
        "Y2": "Oc1ccc(N=Nc2ccccc2)c2ccccc12",  # Sudan I
        "Y3": "Oc1ccc(N=Nc2ccccc2)c2ccccc12",  # Sudan I
        "Y4": "Cc1c(N=Nc2c3ccccc3ccc2O)ccc(C)c1", # Sudan II
        "Y5": "c1ccc(N=Nc2ccc(O)c3ccccc23)c(N=Nc4ccccc4)c1", # Sudan III
        "Y6": "Cc1ccc(N=Nc2ccc(O)c3ccccc23)c(N=Nc4ccc(C)cc4)c1", # Sudan IV
        "Y7": "Oc1c(N=Nc2ccccc2)cc(S(=O)(=O)O)cc1", # Sudan Red G
        "Y8": "O=[N+]([O-])c1ccc(N=Nc2c3ccccc3ccc2O)cc1", # Para Red
        "Y9": "CN(C)c1ccc(N=Nc2ccccc2)cc1", # Dimethyl Yellow
        "Y10": "Nc1ccc(N=Nc2ccccc2)cc1" # Aniline Yellow
    }

    min_energy = float('inf')
    target_molecule_name = None
    
    # Use a dictionary to store calculated energies to avoid recalculating for Y2, Y3
    energies = {}

    for name, smiles in molecules.items():
        if smiles in energies:
            energy = energies[smiles]
        else:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            b_matrix = get_barysz_matrix(mol)
            energy = get_barysz_energy(b_matrix)
            energies[smiles] = energy

        if energy < min_energy:
            min_energy = energy
            target_molecule_name = name

    # Now calculate Moran's I for the target molecule
    target_smiles = molecules[target_molecule_name]
    target_mol = Chem.MolFromSmiles(target_smiles)
    min_moran, max_moran = get_moran_i_min_max(target_mol)

    # Calculate the final product
    product = min_energy * min_moran * max_moran

    print(f"Molecule with the lowest Mass-Weighted Barysz Graph Energy: {target_molecule_name}")
    print(f"Lowest Energy (E): {min_energy}")
    print(f"Minimum Mass-Weighted Moran's I (I_min): {min_moran}")
    print(f"Maximum Mass-Weighted Moran's I (I_max): {max_moran}")
    print("\nFinal Calculation:")
    print(f"Product = E * I_min * I_max")
    print(f"Product = {min_energy} * {min_moran} * {max_moran}")
    print(f"Product = {product}")
    
    # Final answer in the required format
    print(f"\n<<<{product}>>>")

# Execute the solver
solve_challenge()