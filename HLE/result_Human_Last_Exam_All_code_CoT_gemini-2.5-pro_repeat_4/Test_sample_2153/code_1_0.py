import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Increase recursion limit for Hosoya calculation
sys.setrecursionlimit(2000)

# Memoization cache for Hosoya Z index calculation
hosoya_memo = {}

def get_hosoya_z_index(mol):
    """
    Calculates the Hosoya Z index (Z) for a given RDKit molecule object.
    Z is the total number of matchings in the graph.
    Uses a recursive approach with memoization.
    """
    
    # Create a hashable representation of the graph's edges
    edges = frozenset(
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in mol.GetBonds()
    )

    def _calculate_z(graph_edges):
        # Base case: a graph with no edges has one matching (the empty set)
        if not graph_edges:
            return 1
        
        # Return cached result if available
        if graph_edges in hosoya_memo:
            return hosoya_memo[graph_edges]

        # Pick an edge e = {u, v} from the graph
        e = next(iter(graph_edges))
        u, v = tuple(e)

        # Recursive step based on Hosoya's formula: Z(G) = Z(G-e) + Z(G-{u,v})
        
        # 1. Z(G-e): Graph with edge 'e' removed
        graph_minus_e = graph_edges.difference({e})
        res1 = _calculate_z(graph_minus_e)

        # 2. Z(G-{u,v}): Graph with nodes u, v and all incident edges removed
        edges_to_remove = {edge for edge in graph_edges if u in edge or v in edge}
        graph_minus_uv = graph_edges.difference(edges_to_remove)
        res2 = _calculate_z(graph_minus_uv)

        result = res1 + res2
        hosoya_memo[graph_edges] = result
        return result

    return _calculate_z(edges)

def solve_chemistry_problem():
    """
    Executes the full plan to find the target molecule and calculate the final ratio.
    """
    # Step 1: Identify the reference molecule (BCKDH substrate)
    bckdh_substrates = {
        'a-Ketoisovalerate (KIV)': 'CC(C)C(=O)C(=O)O',
        'a-Ketoisocaproate (KIC)': 'CC(C)CC(=O)C(=O)O',
        'a-Keto-b-methylvalerate (KMV)': 'CCC(C)C(=O)C(=O)O'
    }

    substrate_complexities = []
    for name, smiles in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smiles)
        bertz_ct = rdMolDescriptors.BertzCT(mol)
        substrate_complexities.append((bertz_ct, name, smiles))
    
    # Sort by complexity to find the median
    substrate_complexities.sort()
    median_substrate = substrate_complexities[1]
    ref_name = median_substrate[1]
    ref_smiles = median_substrate[2]
    
    # Step 2: Identify the target molecule (deoxynucleoside)
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    ref_balaban_j = rdMolDescriptors.BalabanJ(ref_mol)

    deoxynucleosides = {
        'Deoxyadenosine (dA)': 'NC1=NC=NC2=C1N=CN2C3OC(CO)C(O)C3',
        'Deoxyguanosine (dG)': 'NC1=NC2=C(N=C1O)N=CN2C3OC(CO)C(O)C3',
        # Using a more robust SMILES for dC from PubChem
        'Deoxycytidine (dC)': 'C1C(C(OC1N2C=CC(=NC2=O)N)CO)O',
        'Deoxythymidine (dT)': 'CC1=CN(C(=O)NC1=O)C2OC(CO)C(O)C2'
    }
    
    min_diff = float('inf')
    target_name = None
    target_smiles = None

    for name, smiles in deoxynucleosides.items():
        mol = Chem.MolFromSmiles(smiles)
        balaban_j = rdMolDescriptors.BalabanJ(mol)
        diff = abs(balaban_j - ref_balaban_j)
        if diff < min_diff:
            min_diff = diff
            target_name = name
            target_smiles = smiles
            
    # Step 3 & 4: Calculate indices for the H-included graph of the target molecule
    target_mol_no_h = Chem.MolFromSmiles(target_smiles)
    target_mol_h = Chem.AddHs(target_mol_no_h)

    # Calculate Zagreb(1) Index (M1)
    degrees = [atom.GetDegree() for atom in target_mol_h.GetAtoms()]
    zagreb_m1 = sum(d**2 for d in degrees)

    # Calculate Hosoya Z Index (Z)
    hosoya_z = get_hosoya_z_index(target_mol_h)
    
    # Step 5: Compute and print the final result
    final_ratio = (2 * hosoya_z) / zagreb_m1

    print(f"1. Reference Molecule Identification:")
    print(f"   - BCKDH substrate with median Bertz complexity: {ref_name}")
    print(f"2. Target Molecule Identification:")
    print(f"   - Deoxynucleoside with Balaban J index closest to reference: {target_name}")
    print("-" * 30)
    print(f"Calculations for H-Included {target_name}:")
    print(f"   - Zagreb(1) Index (M1): {zagreb_m1}")
    print(f"   - Hosoya Z Index (Z): {hosoya_z}")
    print("-" * 30)
    print("Final Equation:")
    print(f"Ratio = (2 * Hosoya Z) / Zagreb(1)")
    print(f"Ratio = (2 * {hosoya_z}) / {zagreb_m1}")
    print(f"Ratio = {2 * hosoya_z} / {zagreb_m1} = {final_ratio:.4f}")


if __name__ == '__main__':
    solve_chemistry_problem()