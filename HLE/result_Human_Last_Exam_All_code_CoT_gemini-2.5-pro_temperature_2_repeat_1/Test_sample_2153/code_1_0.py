import rdkit
from rdkit import Chem
from rdkit.Chem import GraphDescriptors
import networkx as nx
import numpy as np

# Dictionary to cache Hosoya Z index calculations for performance
hosoya_memo = {}

def get_molecule(smiles):
    """Generate an RDKit molecule object and add hydrogens."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.AddHs(mol)
    return None

def calculate_bertz_complexity(mol):
    """Calculate the Bertz complexity index."""
    return GraphDescriptors.BertzCT(mol)

def calculate_balaban_j(mol):
    """Calculate the Balaban J index."""
    return GraphDescriptors.BalabanJ(mol)

def calculate_zagreb1_index(mol):
    """Calculate the first Zagreb index (M1)."""
    return sum(atom.GetDegree() ** 2 for atom in mol.GetAtoms())

def mol_to_nx(mol):
    """Convert an RDKit molecule to a NetworkX graph."""
    g = nx.Graph()
    for atom in mol.GetAtoms():
        g.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        g.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return g

def calculate_hosoya_z(graph):
    """
    Calculate the Hosoya Z index of a graph using a recursive formula with memoization.
    Z(G) = Z(G-e) + Z(G-{u,v}) for an edge e=(u,v).
    """
    # Use a canonical graph representation (sorted tuple of sorted edge tuples) for memoization key
    graph_key = tuple(sorted(tuple(sorted(e)) for e in graph.edges()))
    
    if graph_key in hosoya_memo:
        return hosoya_memo[graph_key]
    
    if graph.number_of_edges() == 0:
        return 1

    # Pick an edge to start recursion
    edge = list(graph.edges())[0]
    u, v = edge

    # G-e: Graph with the edge removed
    graph_minus_edge = graph.copy()
    graph_minus_edge.remove_edge(u, v)

    # G-{u,v}: Graph with the edge's vertices (and their incident edges) removed
    graph_minus_nodes = graph.copy()
    graph_minus_nodes.remove_nodes_from([u, v])
    
    result = calculate_hosoya_z(graph_minus_edge) + calculate_hosoya_z(graph_minus_nodes)
    hosoya_memo[graph_key] = result
    return result

def solve():
    """Main function to orchestrate the solution."""
    # Step 1: Identify the reference BCKDH substrate
    bckdh_substrates = {
        "Leucine": "CC(C)CC(C(=O)O)N",
        "Isoleucine": "CCC(C)C(C(=O)O)N",
        "Valine": "CC(C)C(C(=O)O)N",
    }

    substrate_complexities = {}
    for name, smiles in bckdh_substrates.items():
        mol = get_molecule(smiles)
        substrate_complexities[name] = calculate_bertz_complexity(mol)

    median_complexity = np.median(list(substrate_complexities.values()))
    
    target_bcaa_name = [name for name, comp in substrate_complexities.items() if np.isclose(comp, median_complexity)][0]
    target_bcaa_mol = get_molecule(bckdh_substrates[target_bcaa_name])

    # Step 2: Calculate Balaban J for the reference molecule
    target_balaban_j = calculate_balaban_j(target_bcaa_mol)
    
    # Step 3: Identify the molecule of interest from Carrell's synthesis
    carrell_products = {
        "Glycine": "C(C(=O)O)N",
        "Alanine": "CC(C(=O)O)N",
        "Uridine": "C1=CN(C(=O)NC1=O)C2C(C(C(O2)CO)O)O",
        "Cytidine": "C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)O",
    }
    
    product_balaban_j_diffs = {}
    for name, smiles in carrell_products.items():
        mol = get_molecule(smiles)
        balaban_j = calculate_balaban_j(mol)
        product_balaban_j_diffs[name] = abs(balaban_j - target_balaban_j)
    
    final_molecule_name = min(product_balaban_j_diffs, key=product_balaban_j_diffs.get)
    
    # Step 4: Perform final calculations on the identified molecule (Alanine)
    final_mol_smiles = carrell_products[final_molecule_name]
    final_mol = get_molecule(final_mol_smiles)

    # Calculate Zagreb(1) Index
    zagreb1_val = calculate_zagreb1_index(final_mol)
    
    # Calculate Hosoya Z Index
    nx_graph = mol_to_nx(final_mol)
    hosoya_z_val = calculate_hosoya_z(nx_graph)
    
    # Calculate the final ratio
    result = (2 * hosoya_z_val) / zagreb1_val
    
    # Print the result in the specified format
    print(f"The molecule of interest is {final_molecule_name}.")
    print(f"The final calculation is (2 * Hosoya Z index) / Zagreb(1) index.")
    print(f"(2 * {hosoya_z_val}) / {zagreb1_val} = {result}")
    
    # Return final answer
    print(f"\n<<<{result:.5f}>>>")

solve()
