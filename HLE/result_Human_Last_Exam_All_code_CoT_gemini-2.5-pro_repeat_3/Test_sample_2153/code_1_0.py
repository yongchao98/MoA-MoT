import sys
from rdkit import Chem
from rdkit.Chem import GraphDescriptors
import numpy as np
import networkx as nx

def solve_chemistry_riddle():
    """
    This script solves the entire problem by identifying the target molecule
    and performing the final calculations as per the user's request.
    """
    # To handle potential deep recursion in the Hosoya Z calculation
    sys.setrecursionlimit(2000)

    # --- Part 1: Identify the Reference Molecule ---

    # Define the BCKDH complex substrates
    bckdh_substrates = {
        "KIV": "CC(C)C(=O)C(=O)O",    # from Valine
        "KIC": "CC(C)CC(=O)C(=O)O",   # from Leucine
        "KMV": "CCC(C)C(=O)C(=O)O"    # from Isoleucine
    }

    # Calculate Bertz complexity for each substrate
    complexities = {}
    for name, smiles in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smiles)
        complexities[name] = GraphDescriptors.BertzCT(mol)

    # Find the median complexity value
    median_complexity = np.median(list(complexities.values()))

    # Identify the reference molecule(s) with median complexity. Both KIC and KMV match.
    # We select one, as their relevant properties are identical for this problem.
    reference_name = [name for name, c in complexities.items() if c == median_complexity][0]
    reference_smiles = bckdh_substrates[reference_name]
    reference_mol = Chem.MolFromSmiles(reference_smiles)

    # Calculate the Balaban J index for the reference molecule (H-suppressed graph is standard)
    reference_balaban_j = GraphDescriptors.BalabanJ(reference_mol)

    # --- Part 2: Identify the Target Molecule ---

    # Define the candidate molecules from Carrell's 2018 synthesis
    carrell_molecules = {
        "Cytidine": "C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)O",
        "Uridine": "C1=CN(C(=O)NC1=O)C2C(C(C(O2)CO)O)O",
        "Glycine": "C(C(=O)O)N",
        "Alanine": "CC(C(=O)O)N",
        "Aspartic Acid": "C(C(C(=O)O)N)C(=O)O"
    }

    # Calculate Balaban J for each candidate and find the one closest to the reference
    min_diff = float('inf')
    target_molecule_name = None
    target_molecule_smiles = None

    for name, smiles in carrell_molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        j_index = GraphDescriptors.BalabanJ(mol)
        diff = abs(j_index - reference_balaban_j)
        if diff < min_diff:
            min_diff = diff
            target_molecule_name = name
            target_molecule_smiles = smiles

    # --- Part 3: Perform Final Calculations on the Target Molecule ---

    # Create an RDKit molecule object for the target, including hydrogens
    target_mol = Chem.MolFromSmiles(target_molecule_smiles)
    target_mol_h = Chem.AddHs(target_mol)

    # Calculate the Zagreb M1 index (H-included)
    zagreb_m1 = 0
    for atom in target_mol_h.GetAtoms():
        zagreb_m1 += atom.GetDegree() ** 2

    # --- Hosoya Z Index Calculation ---
    hosoya_cache = {}

    def get_graph_id(graph):
        # A frozenset of frozensets of edge endpoints is a hashable, canonical graph ID
        return frozenset(map(frozenset, graph.edges()))

    def calculate_hosoya_z(graph):
        graph_id = get_graph_id(graph)
        if graph_id in hosoya_cache:
            return hosoya_cache[graph_id]

        if graph.number_of_edges() == 0:
            return 1

        # Pick an edge to recurse on
        u, v = list(graph.edges())[0]
        
        # Z(G) = Z(G-e) + Z(G-{u,v})
        graph_minus_e = graph.copy()
        graph_minus_e.remove_edge(u, v)

        graph_minus_uv = graph.copy()
        graph_minus_uv.remove_nodes_from([u, v])

        result = calculate_hosoya_z(graph_minus_e) + calculate_hosoya_z(graph_minus_uv)
        hosoya_cache[graph_id] = result
        return result

    # Convert the H-included RDKit molecule to a NetworkX graph
    G = nx.Graph()
    for atom in target_mol_h.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in target_mol_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # Calculate the Hosoya Z index
    hosoya_z = calculate_hosoya_z(G)

    # Calculate the final ratio
    final_ratio = (2 * hosoya_z) / zagreb_m1

    # --- Part 4: Print the Final Output ---
    print(f"Target Molecule Identified: {target_molecule_name}")
    print("-" * 40)
    print(f"Hosoya Z Index (H-included): {hosoya_z}")
    print(f"Zagreb M1 Index (H-included): {zagreb_m1}")
    print("-" * 40)
    print("Final Equation:")
    print(f"Ratio = (2 * Hosoya_Z) / Zagreb_M1")
    print(f"Ratio = (2 * {hosoya_z}) / {zagreb_m1}")
    print(f"Ratio = {final_ratio}")


if __name__ == "__main__":
    solve_chemistry_riddle()