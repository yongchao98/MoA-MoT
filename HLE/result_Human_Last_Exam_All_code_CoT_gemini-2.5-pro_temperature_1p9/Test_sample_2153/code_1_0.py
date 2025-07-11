import sys
from rdkit import Chem
from rdkit.Chem import GraphDescriptors
import networkx as nx
from itertools import combinations

def solve_cheminformatics_problem():
    """
    This script solves the user's multi-step cheminformatics problem.
    It identifies a target molecule by comparing topological indices and then
    calculates a specific ratio of two other indices for that molecule.
    """

    # Step 1: Define the molecules involved in the problem.
    # BCKDH substrates (derived from Valine, Leucine, Isoleucine).
    bckdh_substrates = {
        'KIV': 'CC(C)C(=O)C(=O)O',  # alpha-ketoisovalerate
        'KIC': 'CC(C)CC(=O)C(=O)O', # alpha-ketoisocaproate
        'KMV': 'CCC(C)C(=O)C(=O)O', # alpha-keto-beta-methylvalerate
    }

    # Assumed synthesized substances based on the prompt's description.
    synthesized_molecules = {
        'Pyrrole': 'c1cc[nH]c1',
        'Pyridine': 'c1ccncc1',
        'Furan': 'c1ccoc1',
        'Thiophene': 'c1ccsc1',
    }

    # Convert SMILES strings to RDKit molecule objects.
    bckdh_mols = {name: Chem.MolFromSmiles(smi) for name, smi in bckdh_substrates.items()}
    synth_mols = {name: Chem.MolFromSmiles(smi) for name, smi in synthesized_molecules.items()}

    # Step 2: Find the BCKDH substrate with the median Bertz's complexity.
    bertz_values = {name: GraphDescriptors.BertzCT(mol) for name, mol in bckdh_mols.items()}
    sorted_bertz = sorted(bertz_values.items(), key=lambda item: item[1])
    median_bckdh_name = sorted_bertz[1][0]
    median_bckdh_mol = bckdh_mols[median_bckdh_name]

    # Step 3: Get the target Balaban J index from the identified BCKDH substrate.
    target_balaban_j = GraphDescriptors.BalabanJ(median_bckdh_mol)

    # Step 4: Find the synthesized molecule with the Balaban J index closest to the target.
    balaban_values_synth = {name: GraphDescriptors.BalabanJ(mol) for name, mol in synth_mols.items()}
    
    closest_molecule_name = min(
        balaban_values_synth.keys(),
        key=lambda name: abs(balaban_values_synth[name] - target_balaban_j)
    )
    final_mol = synth_mols[closest_molecule_name]

    print(f"Identified reference molecule (median Bertz complexity): {median_bckdh_name}")
    print(f"Identified target substance (closest Balaban J): {closest_molecule_name}\n")

    # Step 5: Calculate Zagreb(1) and Hosoya Z for the final molecule.

    # Manual implementation of Zagreb(1) for standard definition.
    def calculate_zagreb1(mol):
        m1 = 0
        for atom in mol.GetAtoms():
            m1 += atom.GetDegree() ** 2
        return m1

    zagreb_m1 = calculate_zagreb1(final_mol)

    # Hosoya Z index calculation for the H-included graph.
    memo_hosoya = {}
    def calculate_hosoya_z(graph):
        edges = frozenset(graph.edges())
        if not edges:
            return 1
        if edges in memo_hosoya:
            return memo_hosoya[edges]

        e = list(edges)[0]
        u, v = e

        graph_minus_e = graph.copy()
        graph_minus_e.remove_edge(u, v)

        graph_minus_uv = graph.copy()
        graph_minus_uv.remove_nodes_from([u, v])

        result = calculate_hosoya_z(graph_minus_e) + calculate_hosoya_z(graph_minus_uv)
        memo_hosoya[edges] = result
        return result

    # Create the H-included graph using NetworkX.
    final_mol_h = Chem.AddHs(final_mol)
    G = nx.Graph()
    for atom in final_mol_h.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in final_mol_h.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    
    memo_hosoya.clear()
    hosoya_z = calculate_hosoya_z(G)

    # Step 6: Calculate and print the final ratio.
    twice_H = 2 * hosoya_z
    
    if zagreb_m1 == 0:
        ratio = float('inf')
        print("Zagreb(1) index is zero, ratio is undefined.")
    else:
        ratio = twice_H / zagreb_m1

    print("--- Final Calculation ---")
    print(f"Molecule: {closest_molecule_name}")
    print(f"Hosoya Z Index (H-included): {hosoya_z}")
    print(f"Zagreb(1) Index: {zagreb_m1}")
    print("\nFinal Equation:")
    print(f"(2 * {hosoya_z}) / {zagreb_m1}")
    print(f"= {twice_H} / {zagreb_m1}")
    print(f"= {ratio}")
    print(f"\n<<<${ratio:.4f}>>>")

if __name__ == "__main__":
    solve_cheminformatics_problem()