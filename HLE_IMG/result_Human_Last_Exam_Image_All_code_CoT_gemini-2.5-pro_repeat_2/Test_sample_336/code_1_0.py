import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen

def solve_eccentric_connectivity_index():
    """
    This function calculates the Sum of Eccentric Connectivity Indices for the
    product molecule depicted in the image, as it satisfies the condition
    of having a Crippen logP > 1.
    """
    # SMILES string for the product molecule:
    # N-(6,7-dihydro-5H-pyrrolo[2,3-d]pyrimidin-4-yl)-5-fluoropyridin-2-amine
    smiles = 'Fc1cc(Nc2nccc3c2CCNC3)ccn1'

    # Create a molecule object from SMILES
    mol = Chem.MolFromSmiles(smiles)

    # The problem requires including hydrogens in the calculation
    mol_h = Chem.AddHs(mol)
    num_atoms = mol_h.GetNumAtoms()

    # Calculate the topological distance matrix for the graph with hydrogens
    dist_matrix = Chem.GetDistanceMatrix(mol_h)

    # Calculate eccentricity for each atom
    # Eccentricity of a vertex is its maximum distance to any other vertex
    eccentricities = np.max(dist_matrix, axis=1).astype(int)

    # Get the degree for each atom
    degrees = [atom.GetDegree() for atom in mol_h.GetAtoms()]

    # Calculate the Eccentric Connectivity Index (ECI)
    total_eci = 0
    equation_parts = []
    
    # ECI is the sum of (degree * eccentricity) over all atoms
    for i in range(num_atoms):
        deg = degrees[i]
        ecc = eccentricities[i]
        term = deg * ecc
        total_eci += term
        equation_parts.append(f"({deg} * {ecc})")
        
    print("Calculating the Sum of Eccentric Connectivity Indices for the product molecule.")
    print(f"SMILES: {smiles}")
    print(f"Formula: {Descriptors.MolFormula(mol_h)}")
    print(f"Number of atoms (including H): {num_atoms}")
    print("\nThe final equation is the sum of (degree * eccentricity) for each atom:")
    # We print the equation with all the numbers as requested.
    # To keep the line length reasonable, we'll print it in parts.
    line_break_every = 10
    full_equation_str = "ECI = "
    for i, part in enumerate(equation_parts):
        full_equation_str += part
        if (i + 1) % line_break_every == 0 and (i + 1) < len(equation_parts):
            full_equation_str += " + \n"
        elif (i + 1) < len(equation_parts):
            full_equation_str += " + "
    print(full_equation_str)
    
    print(f"\nSum of Eccentric Connectivity Indices = {total_eci}")

solve_eccentric_connectivity_index()