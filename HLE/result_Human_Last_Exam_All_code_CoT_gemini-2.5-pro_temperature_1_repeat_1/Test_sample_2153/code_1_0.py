import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.GraphDescriptors import BalabanJ
import networkx as nx

def solve_chemistry_problem():
    """
    This script solves the entire problem by identifying the molecules and calculating the required indices.
    It requires the rdkit-pypi and networkx libraries.
    """
    # Set a higher recursion limit for the Hosoya Z calculation, which is recursive.
    # The default limit might be too low for the size of the molecule.
    try:
        sys.setrecursionlimit(2000)
    except Exception:
        # Some environments might restrict this. The default might be sufficient.
        pass

    # --- Part 1: Identify the reference molecule and its Balaban J index ---

    # Define the BCKDH substrates with non-linear chains
    bckdh_substrates = {
        'KIV': 'CC(C)C(=O)C(=O)O',   # α-ketoisovalerate
        'KIC': 'CC(C)CC(=O)C(=O)O',  # α-ketoisocaproate
        'KMV': 'CCC(C)C(=O)C(=O)O'   # α-keto-β-methylvalerate
    }

    # Calculate Bertz Complexity for each substrate
    bertz_complexities = {}
    for name, smiles in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smiles)
        bertz_complexities[name] = Descriptors.BertzCT(mol)

    # Find the substrate with the median Bertz complexity
    sorted_substrates = sorted(bertz_complexities.items(), key=lambda item: item[1])
    # The list has 3 items, so the median is the 2nd one (index 1)
    median_substrate_name = sorted_substrates[1][0]
    median_substrate_smiles = bckdh_substrates[median_substrate_name]

    # Calculate the Balaban J index for the reference molecule
    ref_mol = Chem.MolFromSmiles(median_substrate_smiles)
    ref_balaban_j = BalabanJ(ref_mol)

    # --- Part 2: Identify the target molecule from Carell's 2018 work ---

    # Define the DNA nucleosides from the 2018 paper
    dna_nucleosides = {
        'deoxyadenosine': 'NC1=NC=NC2=C1N=CN2C1CC(O)C(CO)O1',
        'deoxyguanosine': 'NC1=NC(=O)C2=C(N1)N=CN2C1CC(O)C(CO)O1',
        'deoxycytidine': 'NC1=NC(=O)C=CN1C1CC(O)C(CO)O1',
        'deoxythymidine': 'CC1=CNC(=O)N=C1C1CC(O)C(CO)O1'
    }

    # Find the nucleoside with the Balaban J index closest to the reference
    min_diff = float('inf')
    target_molecule_smiles = None
    
    for name, smiles in dna_nucleosides.items():
        mol = Chem.MolFromSmiles(smiles)
        j_val = BalabanJ(mol)
        diff = abs(j_val - ref_balaban_j)
        if diff < min_diff:
            min_diff = diff
            target_molecule_smiles = smiles

    # --- Part 3: Calculate the final ratio for the target molecule ---

    # Create the H-included molecular graph of the target molecule
    target_mol_h = Chem.AddHs(Chem.MolFromSmiles(target_molecule_smiles))

    # Calculate the Zagreb(1) index (M1)
    m1_index = sum(atom.GetDegree()**2 for atom in target_mol_h.GetAtoms())

    # Create a NetworkX graph to calculate the Hosoya Z index
    G = nx.Graph()
    for atom in target_mol_h.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in target_mol_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Memoization cache for the recursive Hosoya Z function
    memo = {}

    def hosoya_z(graph):
        """Calculates the Hosoya Z index of a graph using a recursive formula with memoization."""
        # Use a frozenset of sorted edge tuples as a canonical graph representation for the cache key
        graph_key = frozenset(map(lambda e: tuple(sorted(e)), graph.edges()))

        if graph_key in memo:
            return memo[graph_key]
        
        # Base case: A graph with no edges has one matching (the empty set).
        if graph.number_of_edges() == 0:
            return 1

        # Recursive step: Z(G) = Z(G-e) + Z(G-{u,v})
        # Pick any edge e = (u,v)
        u, v = next(iter(graph.edges()))

        # Calculate Z(G-e)
        graph_minus_e = graph.copy()
        graph_minus_e.remove_edge(u, v)
        
        # Calculate Z(G-{u,v})
        graph_minus_uv = graph.copy()
        graph_minus_uv.remove_nodes_from([u, v])

        result = hosoya_z(graph_minus_e) + hosoya_z(graph_minus_uv)
        memo[graph_key] = result
        return result

    # Calculate the Hosoya Z index (Z)
    z_index = hosoya_z(G)

    # Calculate the final ratio
    final_ratio = (2 * z_index) / m1_index
    
    # Print the final result in the specified format
    print("The final calculation is based on the target molecule Deoxycytidine.")
    print("Its H-included Zagreb(1) index is {} and its Hosoya Z index is {}.".format(m1_index, z_index))
    print("The final equation is:")
    print("(2 * {}) / {} = {}".format(z_index, m1_index, final_ratio))
    
    return final_ratio

# Execute the main function
final_answer = solve_chemistry_problem()
print(f"\n<<< {final_answer} >>>")
