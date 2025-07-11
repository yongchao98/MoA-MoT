import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen
import numpy as np

def calculate_eccentric_connectivity_index(mol_with_hs):
    """
    Calculates the Eccentric Connectivity Index for a molecule with hydrogens.
    The Eccentric Connectivity Index (ECI) is defined as the sum of products of
    the vertex degree (V_i) and eccentricity (E_i) over all atoms i in the molecular graph.
    ECI = sum(V_i * E_i)
    """
    # Calculate the all-pairs shortest path distance matrix for the graph
    dist_matrix = Chem.GetDistanceMatrix(mol_with_hs)
    
    eci = 0
    num_atoms = mol_with_hs.GetNumAtoms()

    # If the molecule has no atoms, the ECI is 0.
    if num_atoms == 0:
        return 0

    # Iterate over each atom in the molecule
    for i in range(num_atoms):
        # Get the degree of the atom (V_i)
        degree = mol_with_hs.GetAtomWithIdx(i).GetDegree()
        # Get the eccentricity of the atom (E_i), which is the maximum value in its row/column of the distance matrix
        eccentricity = int(np.max(dist_matrix[i]))
        # Add the product to the total ECI sum
        eci += degree * eccentricity
        
    return eci

# Step 1: Define the SMILES strings for the three molecules shown in the image.
smiles_list = [
    'Nc1cc(F)c(OCc2cc(OC)c(OC)c(OC)c2)cn1',  # Reactant A
    'n1cnc2c1CCN(Cl)2',                        # Reactant B
    'Fc1ccc(Nc2ncnc3c2CCN3)nc1'               # Product C
]

# Step 2: Filter molecules with Crippen logP > 1.
eci_components = []
for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    if mol and Crippen.MolLogP(mol) > 1:
        # Step 3: For each filtered molecule, calculate its ECI, including hydrogens.
        # Add explicit hydrogens to the molecular graph as required.
        mol_h = Chem.AddHs(mol)
        eci = calculate_eccentric_connectivity_index(mol_h)
        eci_components.append(eci)

# Step 4: Sum the ECI values of the qualifying molecules.
total_eci = sum(eci_components)

# Step 5: Print the final calculation as an equation.
# The prompt requires showing the numbers that make up the final sum.
equation_parts = [str(x) for x in eci_components]
print(f"{' + '.join(equation_parts)} = {total_eci}")