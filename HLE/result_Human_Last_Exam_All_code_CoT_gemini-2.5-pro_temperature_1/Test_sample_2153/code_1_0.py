import sys
import networkx as nx
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

# Suppress RDKit warnings for cleaner output
from rdkit import rdBase
rdBase.DisableLog('rdApp.warning')

# --- Helper Functions ---

def mol_to_nx(mol):
    """Converts an RDKit molecule to a NetworkX graph."""
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return G

@lru_cache(maxsize=None)
def count_matchings_recursive(edges_tuple):
    """
    Recursively counts the number of matchings in a graph represented by a tuple of edges.
    The Hosoya Z index is the total number of matchings, including the empty one.
    """
    if not edges_tuple:
        return 1

    edges = set(edges_tuple)
    u, v = next(iter(edges))

    # Case 1: Matchings that DO NOT contain the edge (u, v)
    edges1 = edges.copy()
    edges1.remove((u, v))
    # Using tuple(sorted(list(...))) for canonical representation for the cache key
    count1 = count_matchings_recursive(tuple(sorted(list(edges1))))

    # Case 2: Matchings that DO contain the edge (u, v)
    # Remove all edges incident to u or v
    edges2 = edges.copy()
    edges_to_remove = {edge for edge in edges2 if u in edge or v in edge}
    edges2 -= edges_to_remove
    count2 = count_matchings_recursive(tuple(sorted(list(edges2))))

    return count1 + count2

def get_hosoya_z(G):
    """Calculates the Hosoya Z index for a NetworkX graph."""
    # Create a canonical representation of edges (tuple of sorted tuples) for caching
    edges = tuple(sorted([tuple(sorted(edge)) for edge in G.edges()]))
    return count_matchings_recursive(edges)

# --- Main Analysis ---

# Step 1: Identify the BCKDH substrate with median Bertz complexity
bckdh_smiles = {
    "KIV": "CC(C)C(=O)C(=O)O",      # from Valine
    "KIC": "CC(C)CC(=O)C(=O)O",     # from Leucine
    "KMV": "CCC(C)C(=O)C(=O)O"       # from Isoleucine
}

bckdh_data = []
for name, smi in bckdh_smiles.items():
    mol = Chem.MolFromSmiles(smi)
    bertz = GraphDescriptors.BertzCT(mol)
    bckdh_data.append({'name': name, 'mol': mol, 'bertz': bertz})

# Sort by Bertz complexity to find the median
sorted_bckdh = sorted(bckdh_data, key=lambda x: x['bertz'])
median_bckdh_substrate = sorted_bckdh[1]

# Step 2: Calculate Balaban J for the median substrate
balaban_j_bckdh = GraphDescriptors.BalabanJ(median_bckdh_substrate['mol'])

# Step 3: Identify the ribonucleoside with the closest Balaban J index
ribonucleosides_smiles = {
    "Cytidine": "C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O",
    "Uridine": "C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O",
    "Adenosine": "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N=CN=C2N",
    "Guanosine": "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N=C(NC2=O)N"
}

closest_molecule_data = {}
min_diff = float('inf')

for name, smi in ribonucleosides_smiles.items():
    mol = Chem.MolFromSmiles(smi)
    balaban_j = GraphDescriptors.BalabanJ(mol)
    diff = abs(balaban_j - balaban_j_bckdh)
    if diff < min_diff:
        min_diff = diff
        closest_molecule_data = {'name': name, 'mol': mol, 'smi': smi}

target_molecule_name = closest_molecule_data['name']
target_molecule_mol = closest_molecule_data['mol']
target_molecule_smi = closest_molecule_data['smi']

# Step 4: Calculate Zagreb M1 and Hosoya Z for the target molecule

# Calculate Zagreb M1 index (sum of squared degrees of heavy atoms)
zagreb_m1 = sum(atom.GetDegree()**2 for atom in target_molecule_mol.GetAtoms())

# Calculate Hosoya Z index (H-included)
# Add explicit hydrogens to the molecule
mol_h = Chem.AddHs(target_molecule_mol)

# For efficiency, we decompose the graph at the bridge between the base and the sugar
G_h = mol_to_nx(mol_h)
# Find bridge bond (N-C) between Uracil base and Ribose sugar
# In the Uridine SMILES, this connects atom 1 (N) and atom 8 (C)
bridge_bond = (1, 8) 
u, v = bridge_bond

# Decompose the graph
G_decomposed = G_h.copy()
G_decomposed.remove_edge(u, v)
components = list(nx.connected_components(G_decomposed))
c1_nodes, c2_nodes = components[0], components[1]

# Ensure G1 is the smaller component (uracil) and G2 is the larger (ribose)
G1_nodes = c1_nodes if u in c1_nodes else c2_nodes
G2_nodes = c2_nodes if v in c2_nodes else c1_nodes

G1 = G_h.subgraph(G1_nodes)
G2 = G_h.subgraph(G2_nodes)

G1_minus_u = G1.copy()
G1_minus_u.remove_node(u)
G2_minus_v = G2.copy()
G2_minus_v.remove_node(v)

# Hosoya Z is m(G) = m(G-e) + m(G-{u,v}), where m(G-e) = m(G1) * m(G2)
# and m(G-{u,v}) = m(G1-u) * m(G2-v)
m_G1 = get_hosoya_z(G1)
m_G2 = get_hosoya_z(G2)
m_G1_minus_u = get_hosoya_z(G1_minus_u)
m_G2_minus_v = get_hosoya_z(G2_minus_v)

hosoya_z = (m_G1 * m_G2) + (m_G1_minus_u * m_G2_minus_v)

# Step 5: Compute the final ratio
final_ratio = (2 * hosoya_z) / zagreb_m1

# --- Output Results ---
print(f"1. Identified median BCKDH substrate: {median_bckdh_substrate['name']}")
print(f"2. Identified target ribonucleoside (closest Balaban J): {target_molecule_name}")
print("-" * 30)
print("Calculating final ratio for the target molecule...")
print(f"   - Target Molecule: {target_molecule_name}")
print(f"   - Hosoya Z Index (H-included): {hosoya_z}")
print(f"   - Zagreb M1 Index: {zagreb_m1}")
print("-" * 30)
print(f"Final Equation: (2 * {hosoya_z}) / {zagreb_m1}")
print(f"Result: {final_ratio}")

# Final answer in the specified format
sys.stdout.write(f"\n<<<{final_ratio}>>>")