import networkx as nx
from rdkit import Chem

def calculate_indices():
    """
    Identifies the major reduction product of di(perylene-3-yl) disulfide,
    constructs its molecular graph (including H atoms), and calculates the
    Szeged/Wiener index ratio.
    """
    # Note: This script requires the 'rdkit' and 'networkx' libraries.
    # You can install them using pip:
    # pip install rdkit-pypi networkx

    # The reduction of di(perylene-3-yl) disulfide yields perylene-3-thiol.
    # We define this molecule using its SMILES string.
    smiles = "C1=CC=C2C3=C(C=C1)C=C4C5=C(C=C(C=C5)S)C=C3C2=C4"
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("Error: Could not generate molecule from SMILES string.")
        return

    # Add explicit hydrogen atoms to the molecule, as required by the problem.
    mol_with_hs = Chem.AddHs(mol)

    # Convert the molecule to a NetworkX graph for topological analysis.
    G = nx.Graph()
    for atom in mol_with_hs.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol_with_hs.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    print(f"Analyzing the major reduction product: perylene-3-thiol (C20H12S)")
    print(f"Total number of atoms (graph nodes): {G.number_of_nodes()}")
    print("-" * 30)

    # Calculate the Wiener Index (W)
    wiener_index = nx.wiener_index(G)
    print(f"Wiener Index (W) = {wiener_index}")

    # Calculate the Szeged Index (Sz)
    szeged_index = 0
    # Pre-calculating all-pairs shortest paths is more efficient for the next step.
    all_paths = dict(nx.all_pairs_shortest_path_length(G))

    for u, v in G.edges():
        n_u = 0
        n_v = 0
        dist_u = all_paths[u]
        dist_v = all_paths[v]

        for w in G.nodes():
            if dist_u[w] < dist_v[w]:
                n_u += 1
            elif dist_v[w] < dist_u[w]:
                n_v += 1
        szeged_index += n_u * n_v

    print(f"Szeged Index (Sz) = {szeged_index}")
    print("-" * 30)

    # Calculate and print the final ratio as an equation.
    if wiener_index > 0:
        ratio = szeged_index / wiener_index
        print(f"Szeged/Wiener Ratio Equation:")
        print(f"{szeged_index} / {wiener_index} = {ratio}")
    else:
        print("Cannot calculate ratio: Wiener index is zero.")

if __name__ == "__main__":
    calculate_indices()