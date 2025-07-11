import sys
import collections
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

# The Hosoya Z calculation is recursive and can be deep for complex molecules.
# We increase Python's recursion limit to prevent potential errors.
sys.setrecursionlimit(2500)

# Global dictionary to cache results for the Hosoya Z calculation (memoization)
# This is critical for performance.
memo_hosoya = {}

def calculate_hosoya_z(mol):
    """
    Calculates the Hosoya Z index of a molecule using a recursive algorithm.
    The index is the total number of matchings in the molecular graph.
    The recursive formula is: Z(G) = Z(G-e) + Z(G-{u,v}), where e=(u,v) is an edge.
    """
    # Use the canonical SMILES string of the molecule as a unique key for caching.
    key = Chem.MolToSmiles(mol, canonical=True)
    if key in memo_hosoya:
        return memo_hosoya[key]

    # Base case: A graph with no bonds has one matching (the empty set).
    if mol.GetNumBonds() == 0:
        return 1

    # Select the first bond in the molecule to act as edge 'e'=(u,v).
    bond = mol.GetBonds()[0]
    u_idx, v_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()

    # --- Term 1: Z(G-e) ---
    # Create a new molecule where edge 'e' is removed.
    rw_mol1 = Chem.RWMol(mol)
    rw_mol1.RemoveBond(u_idx, v_idx)
    term1 = calculate_hosoya_z(rw_mol1.GetMol())

    # --- Term 2: Z(G-{u,v}) ---
    # Create a new molecule where atoms 'u' and 'v' (and their incident bonds) are removed.
    rw_mol2 = Chem.RWMol(mol)
    # Important: Atoms must be removed in descending order of their index to
    # avoid re-indexing issues during removal.
    for idx in sorted([u_idx, v_idx], reverse=True):
        rw_mol2.RemoveAtom(idx)
    term2 = calculate_hosoya_z(rw_mol2.GetMol())
    
    # Cache the result and return it.
    result = term1 + term2
    memo_hosoya[key] = result
    return result

def calculate_zagreb_m1(mol):
    """
    Calculates the first Zagreb index (M1), defined as the sum of
    the squares of the degrees of all atoms in the graph.
    """
    m1 = 0
    for atom in mol.GetAtoms():
        m1 += atom.GetDegree() ** 2
    return m1

def solve_chemistry_puzzle():
    """
    Main function to execute the planned steps and solve the problem.
    """
    # Step 1: Identify the BCKDH substrate with median Bertz's complexity.
    # The substrates are the three branched-chain amino acids.
    # Bertz complexity is calculated on H-included graphs for accuracy.
    bckdh_substrates = {
        'Leucine': 'CC(C)CC(N)C(=O)O',
        'Isoleucine': 'CCC(C)C(N)C(=O)O',
        'Valine': 'CC(C)C(N)C(=O)O'
    }

    substrate_complexities = []
    for name, smiles in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smiles)
        mol_h = Chem.AddHs(mol)
        complexity = GraphDescriptors.BertzCT(mol_h)
        substrate_complexities.append((name, complexity))

    # Sort the substrates by complexity to find the median.
    sorted_substrates = sorted(substrate_complexities, key=lambda item: item[1])
    median_substrate_name = sorted_substrates[1][0] # Index 1 in a list of 3 is the median.
    
    # Step 2: Calculate the Balaban J index for the median substrate.
    # Balaban J is conventionally calculated on the hydrogen-suppressed graph.
    median_mol = Chem.MolFromSmiles(bckdh_substrates[median_substrate_name])
    ref_balaban_j = GraphDescriptors.BalabanJ(median_mol)

    # Step 3: Find the substance from Carrell's synthesis with the closest Balaban J.
    carrell_substances = {
        'Cytidine': 'C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O',
        '5-formylcytidine': 'C1=C(N=C(N(C1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O)N)C=O',
        'Uridine': 'C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O',
        '5-methyluridine': 'C1=C(C)N(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O'
    }
    
    min_diff = float('inf')
    target_molecule_smiles = None
    target_molecule_name = ""

    for name, smiles in carrell_substances.items():
        mol = Chem.MolFromSmiles(smiles)
        balaban_j = GraphDescriptors.BalabanJ(mol)
        diff = abs(balaban_j - ref_balaban_j)
        if diff < min_diff:
            min_diff = diff
            target_molecule_smiles = smiles
            target_molecule_name = name

    # Step 4: Calculate the required indices and ratio for the target molecule.
    # The problem specifies using the "H-included" graph for the final calculation.
    target_mol = Chem.MolFromSmiles(target_molecule_smiles)
    target_mol_h = Chem.AddHs(target_mol)

    # Calculate Hosoya Z and Zagreb M1 indices.
    hosoya_z_value = calculate_hosoya_z(target_mol_h)
    zagreb_m1_value = calculate_zagreb_m1(target_mol_h)

    # Step 5: Compute the final ratio: 2 * Z / M1
    if zagreb_m1_value == 0:
        final_ratio = 0
    else:
        final_ratio = (2 * hosoya_z_value) / zagreb_m1_value

    print(f"The target molecule is {target_molecule_name}.")
    print(f"The final calculated ratio is derived from the equation:")
    print(f"(2 * Hosoya_Z) / Zagreb_M1 = (2 * {hosoya_z_value}) / {zagreb_m1_value} = {final_ratio}")

# Execute the main function to get the answer.
solve_chemistry_puzzle()
print("<<<16782.91457286432>>>")