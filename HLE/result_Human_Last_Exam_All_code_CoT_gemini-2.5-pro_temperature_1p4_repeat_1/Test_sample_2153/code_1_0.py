import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, GraphDescriptors, AllChem, rdMolDescriptors
import networkx as nx
import numpy as np

def solve_cheminformatics_problem():
    """
    Solves the entire multi-step cheminformatics problem.
    """

    # Memoization cache for the recursive Hosoya Z calculation
    hosoya_memo = {}

    def count_matchings(graph):
        """
        Calculates the Hosoya Z index (number of matchings) of a graph
        using a recursive deletion-contraction method with memoization.
        """
        # A canonical representation of the graph (as a frozenset of its edges) is used as a key for memoization.
        graph_key = frozenset(tuple(sorted(edge)) for edge in graph.edges())

        if graph_key in hosoya_memo:
            return hosoya_memo[graph_key]
        
        # Base case: A graph with no edges has one matching (the empty set).
        if not graph.edges():
            return 1

        # Pick a deterministic edge to start the recursion.
        edge = sorted(list(graph.edges()))[0]
        u, v = edge

        # Z(G) = Z(G-e) + Z(G-{u,v})
        # Case 1: The edge `e` is NOT in the matching.
        # We calculate matchings in the graph G-e (graph with the edge removed).
        graph_minus_e = graph.copy()
        graph_minus_e.remove_edge(u, v)
        res1 = count_matchings(graph_minus_e)

        # Case 2: The edge `e` IS in the matching.
        # We must remove nodes u and v and all their incident edges, then find matchings in the remaining graph.
        graph_minus_uv = graph.copy()
        nodes_to_remove = [u, v]
        # This check is necessary for nodes that might already be removed in a recursive step
        existing_nodes_to_remove = [node for node in nodes_to_remove if node in graph_minus_uv]
        graph_minus_uv.remove_nodes_from(existing_nodes_to_remove)
        res2 = count_matchings(graph_minus_uv)
        
        result = res1 + res2
        hosoya_memo[graph_key] = result
        return result

    # Step 1: Identify BCKDH substrate with median Bertz complexity.
    # Substrates are the branched-chain amino acids.
    bckdh_substrates = {
        'Valine': 'C[C@H](C)[C@H](N)C(=O)O',
        'Leucine': 'C[C@H](C)C[C@H](N)C(=O)O',
        'Isoleucine': 'CC[C@H](C)[C@H](N)C(=O)O'
    }

    substrate_data = []
    for name, smiles in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            bertz_ct = Descriptors.BertzCT(mol)
            substrate_data.append({'name': name, 'smiles': smiles, 'bertz_ct': bertz_ct})

    substrate_data.sort(key=lambda x: x['bertz_ct'])
    median_substrate = substrate_data[len(substrate_data) // 2]
    
    # Step 2: Calculate Balaban J for the reference molecule (Leucine).
    ref_mol = Chem.MolFromSmiles(median_substrate['smiles'])
    target_balaban_j = GraphDescriptors.BalabanJ(ref_mol)

    # Step 3: Identify the nucleobase from Carrell's work with the closest Balaban J.
    # These are the canonical nucleobases.
    carrell_substances = {
        'Adenine': 'Nc1ncnc2n[cH]nc12',
        'Guanine': 'Nc1nc2[nH]c(=O)nc2[nH]1',
        'Cytosine': 'Nc1[nH]c(=O)cc[nH]1',
        'Uracil': 'O=c1[nH]ccc(=O)[nH]1'
    }

    closest_substance = None
    min_diff = float('inf')
    for name, smiles in carrell_substances.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            balaban_j = GraphDescriptors.BalabanJ(mol)
            diff = abs(balaban_j - target_balaban_j)
            if diff < min_diff:
                min_diff = diff
                closest_substance = {'name': name, 'smiles': smiles, 'mol': mol}
    
    final_substance_mol = closest_substance['mol']
    final_substance_name = closest_substance['name']
    
    # Step 4: Perform final calculations for the identified substance (Guanine).
    mol_h = Chem.AddHs(final_substance_mol)

    # Calculate Zagreb(1) index (M1), H-included.
    degrees = [atom.GetDegree() for atom in mol_h.GetAtoms()]
    zagreb_m1 = sum(d**2 for d in degrees)

    # Create a NetworkX graph from the H-included molecule.
    G = nx.Graph()
    for atom in mol_h.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Calculate Hosoya Z index (H-included).
    hosoya_z = count_matchings(G)

    # Calculate the final ratio.
    final_ratio = (2 * hosoya_z) / zagreb_m1
    
    # Print the final result in the requested format.
    print(f"The selected molecule is {final_substance_name}.")
    print(f"Its twice Hosoya Z (H-included) to Zagreb(1) index ratio is calculated as follows:")
    print(f"(2 * Hosoya Z) / Zagreb(1) = (2 * {hosoya_z}) / {zagreb_m1} = {final_ratio:.4f}")

solve_cheminformatics_problem()