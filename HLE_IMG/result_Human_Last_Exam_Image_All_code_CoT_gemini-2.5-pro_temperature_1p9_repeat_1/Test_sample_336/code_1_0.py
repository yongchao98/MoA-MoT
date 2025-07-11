import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Crippen
    import numpy as np
except ImportError:
    print("RDKit is not installed. Please install it, for example using 'pip install rdkit-pypi'")
    sys.exit(1)

def solve_chemistry_problem():
    """
    This function calculates the Crippen logP and Eccentric Connectivity Index for a given molecule.
    """
    # The product molecule identified from the image is N-(5-fluoropyridin-2-yl)-6,7-dihydro-5H-pyrrolo[2,3-d]pyrimidin-4-amine.
    smiles = "Fc1cc(Nc2ncnc3c2CCNC3)ncc1"
    
    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Error: Could not create molecule from SMILES string: {smiles}")
        return

    # 1. Calculate and check the Crippen logP value
    logp = Crippen.MolLogP(mol)
    print(f"The product molecule is: N-(5-fluoropyridin-2-yl)-6,7-dihydro-5H-pyrrolo[2,3-d]pyrimidin-4-amine")
    print(f"SMILES: {smiles}")
    print(f"Calculated Crippen logP: {logp:.4f}")

    if logp <= 1:
        print("Skipping ECI calculation as Crippen logP is not > 1.")
        return

    print("Crippen logP > 1. Proceeding with ECI calculation.\n")

    # 2. Add explicit hydrogens for the ECI calculation
    mol_h = Chem.AddHs(mol)
    num_atoms = mol_h.GetNumAtoms()
    print(f"Total number of atoms (including hydrogens): {num_atoms}")

    # 3. Calculate the Eccentric Connectivity Index (ECI)
    # Get the all-pairs shortest path distance matrix (topological distance)
    dist_matrix = Chem.GetDistanceMatrix(mol_h)
    
    # Get the degree for each atom (V_i)
    degrees = [atom.GetDegree() for atom in mol_h.GetAtoms()]
    
    # Get the eccentricity for each atom (E_i), which is the max value in each row of the distance matrix
    eccentricities = [int(np.max(row)) for row in dist_matrix]

    # Calculate the sum of (V_i * E_i)
    total_eci = 0
    equation_parts = []
    for i in range(num_atoms):
        v = degrees[i]
        e = eccentricities[i]
        term = v * e
        total_eci += term
        equation_parts.append(f"{v}*{e}")

    # Print the full equation and the final result
    print("The Eccentric Connectivity Index (ECI) is calculated as Sum(V_i * E_i):")
    print("ECI = " + " + ".join(equation_parts))
    print(f"\nSum of Eccentric Connectivity Indices: {total_eci}")
    print(f"<<<{total_eci}>>>")


# Run the function
solve_chemistry_problem()