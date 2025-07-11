import networkx as nx
from rdkit import Chem
from rdkit.Chem import GraphDescriptors
from functools import lru_cache
import numpy as np

def solve_cheminformatics_problem():
    """
    Solves the multi-step cheminformatics problem to find the specified molecular index ratio.
    """
    # Step 1: Identify the reference BCKDH substrate and its Balaban J index
    # Define BCKDH substrates (branched-chain alpha-keto acids) by their SMILES strings
    bckdh_substrates = {
        'KIV': 'CC(C)C(=O)C(=O)O',  # alpha-ketoisovalerate (from Valine)
        'KIC': 'CC(C)CC(=O)C(=O)O', # alpha-ketoisocaproate (from Leucine)
        'KMV': 'CCC(C)C(=O)C(=O)O'   # alpha-keto-beta-methylvalerate (from Isoleucine)
    }

    # Calculate Bertz complexity for each substrate to find the median one
    bertz_complexities = {}
    for name, smiles in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            bertz_complexities[name] = GraphDescriptors.BertzCT(mol)

    # Sort substrates by complexity to find the median
    sorted_substrates = sorted(bertz_complexities.items(), key=lambda item: item[1])
    median_substrate_name = sorted_substrates[len(sorted_substrates) // 2][0]
    median_substrate_smiles = bckdh_substrates[median_substrate_name]
    
    # Calculate the target Balaban J index from the reference molecule
    ref_mol = Chem.MolFromSmiles(median_substrate_smiles)
    target_balaban_j = GraphDescriptors.BalabanJ(ref_mol)

    # Step 2 & 3: Identify the target molecule from the Carrell 2018 discovery
    # The paper describes simultaneous synthesis of two pairs of ribonucleosides.
    candidate_molecules = {
        'Uridine': 'C1=CN(C(=O)NC1=O)C2C(C(C(O2)CO)O)O',
        'Cytidine': 'C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)O',
        'Adenosine': 'c1nc(c2c(n1)n(cn2)C3C(C(C(O3)CO)O)O)N',
        'Inosine': 'c1nc2c(c(=O)[nH]c1)n(cn2)C3C(C(C(O3)CO)O)O'
    }

    # Find the candidate whose Balaban J index is closest to the target
    best_candidate_name = None
    min_diff = float('inf')

    for name, smiles in candidate_molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            balaban_j = GraphDescriptors.BalabanJ(mol)
            diff = abs(balaban_j - target_balaban_j)
            if diff < min_diff:
                min_diff = diff
                best_candidate_name = name

    target_molecule_name = best_candidate_name
    target_molecule_smiles = candidate_molecules[target_molecule_name]

    # Step 4: Calculate indices for the identified target molecule
    target_mol = Chem.MolFromSmiles(target_molecule_smiles)

    # --- Index Calculation Functions ---

    def calculate_zagreb_m1(mol):
        """Calculates the Zagreb M1 index."""
        m1 = 0
        for atom in mol.GetAtoms():
            m1 += atom.GetDegree() ** 2
        return m1

    def smiles_to_nx_graph(smiles):
        """Converts an RDKit molecule to a NetworkX graph."""
        mol = Chem.MolFromSmiles(smiles)
        G = nx.Graph()
        for atom in mol.GetAtoms():
            G.add_node(atom.GetIdx())
        for bond in mol.GetBonds():
            G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        return G

    def get_edge_frozenset(G):
        """Creates a hashable representation of a graph's edges."""
        return frozenset(map(frozenset, G.edges()))

    @lru_cache(maxsize=None)
    def calculate_hosoya_z(graph_edges_fs):
        """Calculates the Hosoya Z index using a recursive formula with memoization."""
        if not graph_edges_fs:
            return 1  # Z(empty graph) = 1

        G = nx.Graph(list(map(list, graph_edges_fs)))
        u, v = next(iter(G.edges())) # Pick an edge e = (u,v)

        # Recurrence relation: Z(G) = Z(G-e) + Z(G-{u,v})
        G_minus_e = G.copy()
        G_minus_e.remove_edge(u, v)
        z1 = calculate_hosoya_z(get_edge_frozenset(G_minus_e))

        G_minus_uv = G.copy()
        G_minus_uv.remove_node(u)
        G_minus_uv.remove_node(v)
        z2 = calculate_hosoya_z(get_edge_frozenset(G_minus_uv))

        return z1 + z2

    # --- Final Calculations ---
    zagreb_m1_index = calculate_zagreb_m1(target_mol)
    nx_graph = smiles_to_nx_graph(target_molecule_smiles)
    hosoya_z_index = calculate_hosoya_z(get_edge_frozenset(nx_graph))
    
    # Compute the final ratio
    final_ratio = (2 * hosoya_z_index) / zagreb_m1_index

    # --- Print Results ---
    print(f"The identified target molecule is {target_molecule_name}.")
    print(f"Hosoya Z Index (Z): {hosoya_z_index}")
    print(f"Zagreb(1) Index (M1): {zagreb_m1_index}")
    print("\nFinal Equation:")
    print(f"(2 * Z) / M1 = (2 * {hosoya_z_index}) / {zagreb_m1_index} = {final_ratio}")

if __name__ == '__main__':
    solve_cheminformatics_problem()