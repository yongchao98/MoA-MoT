import numpy as np
from rdkit import Chem
from scipy.sparse.csgraph import shortest_path

def calculate_eccentric_connectivity_index():
    """
    Calculates the Eccentric Connectivity Index for a molecule specified by a SMILES string.
    The calculation includes hydrogen atoms.
    """
    # Step 1: Define the target molecule.
    # Based on the analysis, Molecule A is chosen as it unambiguously has a Crippen logP > 1.
    # SMILES for Molecule A: 3-((3,4,5-trimethoxybenzyl)oxy)-5-fluoropyridin-2-amine
    smiles = "COc1cc(COC2=C(N)N=C(F)C=C2)cc(OC)c1OC"
    
    # Step 2: Create a molecule object and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: Invalid SMILES string.")
        return

    mol_h = Chem.AddHs(mol)
    num_atoms = mol_h.GetNumAtoms()

    # Step 3: Get the adjacency matrix and compute the distance matrix.
    # The distance matrix contains the shortest path lengths between all pairs of atoms.
    adj_matrix = Chem.GetAdjacencyMatrix(mol_h)
    dist_matrix = shortest_path(csgraph=adj_matrix, directed=False, unweighted=True)

    # Check for disconnected fragments, which would result in infinite distances.
    if np.any(dist_matrix == np.inf):
        print("Error: The molecule graph is not connected.")
        return

    # Step 4: Calculate the Eccentric Connectivity Index.
    # ECI = sum(degree(i) * eccentricity(i)) for all atoms i.
    eccentric_connectivity_index = 0
    equation_terms = []

    for i in range(num_atoms):
        # Get the degree of the atom (number of connections).
        degree = mol_h.GetAtomWithIdx(i).GetDegree()
        
        # Get the eccentricity of the atom (max shortest path to any other atom).
        eccentricity = int(np.max(dist_matrix[i]))
        
        # Calculate the term for this atom and add it to the total sum.
        term = degree * eccentricity
        eccentric_connectivity_index += term
        equation_terms.append(str(term))

    # Step 5: Print the results as requested.
    print("The Eccentric Connectivity Index is the sum of the products of atom degree and eccentricity for all 42 atoms (including hydrogens).")
    print("The equation is:")
    # To keep the output clean, we show the first few and last few terms.
    if len(equation_terms) > 20:
        print(" + ".join(equation_terms[:10]) + " + ... + " + " + ".join(equation_terms[-10:]) + f" = {eccentric_connectivity_index}")
    else:
        print(" + ".join(equation_terms) + f" = {eccentric_connectivity_index}")
    
    print(f"\nThe Sum of Eccentric Connectivity Indices is: {eccentric_connectivity_index}")
    
    # The final answer in the required format
    print(f"\n<<< {eccentric_connectivity_index} >>>")

# Execute the calculation
calculate_eccentric_connectivity_index()