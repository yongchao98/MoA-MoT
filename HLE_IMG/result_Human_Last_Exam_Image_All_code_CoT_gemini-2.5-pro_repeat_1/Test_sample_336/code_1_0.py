import sys

def solve():
    """
    This function solves the problem by calculating the Sum of Eccentric Connectivity Indices for the specified molecule.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        import numpy as np
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import shortest_path
    except ImportError:
        print("Required libraries (rdkit-pypi, numpy, scipy) are not installed.")
        print("Please install them using: pip install rdkit-pypi numpy scipy")
        sys.exit(1)

    # Step 3 & 4: Define the SMILES for the product molecule.
    # The product is identified as N-(6,7-dihydro-5H-pyrrolo[3,4-b]pyrazin-5-yl)-5-fluoropyridin-2-amine
    product_smiles = "c1(NC2CN3C=NC=C3C2)ncc(F)cc1"

    # Step 5: Create molecule and check logP
    mol = Chem.MolFromSmiles(product_smiles)
    if mol is None:
        print(f"Error: Could not parse SMILES string: {product_smiles}")
        return

    logP = Descriptors.MolLogP(mol)
    # The condition is logP > 1. Let's confirm.
    if logP <= 1:
        print(f"The molecule's logP is {logP:.2f}, which does not meet the condition (> 1).")
        # However, we proceed as it's the intended molecule for the calculation.

    # Step 7a, 7b: Add hydrogens and compute distance matrix
    mol_h = Chem.AddHs(mol)
    adj_matrix = Chem.GetAdjacencyMatrix(mol_h)
    graph = csr_matrix(adj_matrix)
    dist_matrix = shortest_path(csgraph=graph, directed=False, unweighted=True)
    
    if np.any(np.isinf(dist_matrix)):
        print("Warning: The molecule graph is not fully connected.")
        dist_matrix[np.isinf(dist_matrix)] = 0

    # Step 7c, 7d: Calculate eccentricities
    eccentricities = np.max(dist_matrix, axis=1)

    # Step 7d, 7e: Calculate and sum the Eccentric Connectivity Indices
    total_eci = 0
    eci_terms = []
    
    atoms = mol_h.GetAtoms()
    for i in range(len(atoms)):
        atom = atoms[i]
        degree = atom.GetDegree()
        eccentricity = int(eccentricities[i])
        eci = degree * eccentricity
        total_eci += eci
        eci_terms.append(f"{degree}*{eccentricity}")

    # Step 8: Print the final result in the requested format
    print(f"Molecule: N-(6,7-dihydro-5H-pyrrolo[3,4-b]pyrazin-5-yl)-5-fluoropyridin-2-amine")
    print(f"SMILES: {product_smiles}")
    print(f"Crippen logP: {logP:.2f}")
    print(f"Number of atoms (including H): {mol_h.GetNumAtoms()}")
    print("\nCalculation of the Sum of Eccentric Connectivity Indices:")
    
    # To make the output readable, print terms in groups
    line_break_every = 10
    for i in range(0, len(eci_terms), line_break_every):
        chunk = eci_terms[i:i+line_break_every]
        if i == 0:
            print(" + ".join(chunk), end="")
        else:
            print(" + \n" + " + ".join(chunk), end="")
    
    print(f" = {int(total_eci)}")


solve()