import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.linalg import eigvals

def get_molecule_from_string(y_str, mapping):
    """Converts a formula string like 'X1X2...' into an RDKit molecule."""
    mol = Chem.RWMol()
    atoms = []
    # First, add all atoms to the molecule
    for i in range(0, len(y_str), 2):
        x_key = y_str[i:i+2]
        if x_key not in mapping:
            # Handle cases like X3, X4 etc.
            x_key = y_str[i:i+1] if i + 1 == len(y_str) or not y_str[i+1].isdigit() else y_str[i:i+2]
            if not y_str[i:i+2].startswith('X'):
                 x_key = y_str[i:i+1]
            else:
                 x_key = y_str[i:i+2]
            
            # Simplified logic for parsing single digit X's
            if len(x_key) == 1 and x_key.startswith('X'):
                pass # Should not happen with given Ys
            elif not x_key.startswith('X'):
                 key_num_str = y_str[i+1]
                 x_key = f"X{key_num_str}"
            
        element_symbol = mapping[x_key]
        atom = Chem.Atom(element_symbol)
        idx = mol.AddAtom(atom)
        atoms.append(idx)

    # Then, create bonds for the linear chain
    for i in range(len(atoms) - 1):
        mol.AddBond(atoms[i], atoms[i+1], Chem.BondType.SINGLE)
    
    # Finalize the molecule
    final_mol = mol.GetMol()
    # It's good practice to sanitize and add hydrogens if needed, but for these structures it's straightforward
    Chem.SanitizeMol(final_mol)
    return final_mol

def get_mass_weighted_barysz_energy(mol):
    """Calculates the Mass-Weighted Barysz Graph Energy for a molecule."""
    adj_matrix = Chem.GetAdjacencyMatrix(mol)
    
    masses = np.array([atom.GetMass() for atom in mol.GetAtoms()])
    if 0 in masses: # Should not happen with our elements
        return float('inf')

    # Laplacian Matrix L = D - A
    degrees = np.sum(adj_matrix, axis=1)
    laplacian_matrix = np.diag(degrees) - adj_matrix
    
    # Weight matrix W and its inverse square root W_inv_sqrt
    w_inv_sqrt = np.diag(1.0 / np.sqrt(masses))
    
    # Barysz Matrix B = W^(-1/2) * L * W^(-1/2)
    barysz_matrix = w_inv_sqrt @ laplacian_matrix @ w_inv_sqrt
    
    # Energy is the sum of the absolute values of the eigenvalues
    e_vals = eigvals(barysz_matrix)
    energy = np.sum(np.abs(e_vals.real)) # Eigenvalues may have small imaginary parts due to precision
    
    return energy

def get_mass_weighted_moran_i_min_max(mol):
    """Calculates Moran's I for different distances and returns the min and max."""
    num_atoms = mol.GetNumAtoms()
    if num_atoms <= 1:
        return 0.0, 0.0

    masses = np.array([atom.GetMass() for atom in mol.GetAtoms()])
    mean_mass = np.mean(masses)
    
    # Deviations from the mean
    deviations = masses - mean_mass
    
    # Denominator of Moran's I formula
    denom = np.sum(deviations**2)
    if denom == 0:
        return 0.0, 0.0

    # Get topological distance matrix
    dist_matrix = Chem.GetDistanceMatrix(mol)
    
    moran_i_values = []
    
    # Calculate Moran's I for each possible lag distance d
    for d in range(1, num_atoms):
        # Weight matrix Wd for distance d
        wd = (dist_matrix == d).astype(float)
        s0 = np.sum(wd)
        
        if s0 == 0:
            continue
            
        # Numerator calculation
        num = 0.0
        for i in range(num_atoms):
            for j in range(num_atoms):
                num += wd[i, j] * deviations[i] * deviations[j]
        
        # Moran's I formula
        moran_i = (num_atoms / s0) * (num / denom)
        moran_i_values.append(moran_i)

    if not moran_i_values:
        return 0.0, 0.0
        
    return min(moran_i_values), max(moran_i_values)


# 1. Decipher the elements
single_letter_elements = sorted(['B', 'C', 'F', 'H', 'I', 'K', 'N', 'O', 'P', 'S', 'U', 'V', 'W', 'Y'])
x_mapping = {f"X{i+1}": single_letter_elements[i] for i in range(len(single_letter_elements))}

# 2. Define Y formulas
y_formulas = {
    "Y1": "X1X2X3X4X5X6",
    "Y2": "X8X9X10X11X5X6X12X8X9X8X12X13",
    "Y3": "X11X3X4X14X3X4X5X6",
    "Y4": "X12X4X13X5X6X3",
    "Y5": "X8X9X10X11X9X14X5X6X3",
    "Y6": "X1X10X5X1X9X4X3",
    "Y7": "X8X9X10X11X12X4X5X6",
    "Y8": "X10X2X5X13X9X4X12X4X3",
    "Y9": "X9X14X5X11X3X4X14X3X4X3",
    "Y10": "X1X12X1X3X10X12X13X12X4X3",
}

# Add single-digit keys to mapping for parsing (e.g., in "X5X6X3", the "3" needs to be handled)
for i in range(1,10):
    if f"X{i}" in x_mapping:
         y_formulas = {k: v.replace(f"X{i}", f"X0{i}") for k, v in y_formulas.items()}
x_mapping = {f"X{str(i+1).zfill(2)}": single_letter_elements[i] for i in range(len(single_letter_elements))}

y_formulas_mod = {}
for name, formula in y_formulas.items():
    new_formula = ""
    i = 0
    while i < len(formula):
        if formula[i] == 'X':
            # Find the end of the number
            j = i + 1
            while j < len(formula) and formula[j].isdigit():
                j += 1
            num_str = formula[i+1:j]
            new_formula += f"X{num_str.zfill(2)}"
            i = j
        else: # Should not happen
            i+=1
    y_formulas_mod[name] = new_formula


# 3 & 4. Calculate energies and find the minimum
min_energy = float('inf')
min_energy_molecule_name = None
results = {}

for name, formula_str in y_formulas_mod.items():
    mol = get_molecule_from_string(formula_str, x_mapping)
    if mol:
        energy = get_mass_weighted_barysz_energy(mol)
        results[name] = {'mol': mol, 'energy': energy}
        if energy < min_energy:
            min_energy = energy
            min_energy_molecule_name = name

# 5. Calculate Moran's I for the identified molecule
target_molecule_name = min_energy_molecule_name
target_molecule = results[target_molecule_name]['mol']
target_energy = results[target_molecule_name]['energy']

min_moran_i, max_moran_i = get_mass_weighted_moran_i_min_max(target_molecule)

# 6. Compute the final product
final_product = target_energy * min_moran_i * max_moran_i

# Print the results as requested
print(f"Deciphered Element Mapping (alphabetical): {x_mapping}")
print(f"\nIdentified molecule with the lowest Mass-Weighted Barysz Graph Energy: {target_molecule_name}")
print(f"Identified Energy = {target_energy}")
print(f"Minimum Mass-Weighted Moran's I = {min_moran_i}")
print(f"Maximum Mass-Weighted Moran's I = {max_moran_i}")
print(f"\nFinal Calculation: {target_energy} * {min_moran_i} * {max_moran_i}")

print(f"\nFinal Product: {final_product}")
print(f"<<<{final_product}>>>")