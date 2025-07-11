import sys
import networkx as nx
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

def solve_task():
    """
    This script solves the multi-step chemoinformatics problem by:
    1. Identifying the target BCKDH substrate based on median Bertz's complexity.
    2. Identifying the target ribonucleoside based on the closest Balaban J index.
    3. Calculating the (2 * Hosoya Z) / Zagreb(1) ratio for the H-included graph of the target ribonucleoside.
    Requires rdkit-pypi and networkx. Install with: pip install rdkit-pypi networkx
    """

    # Step 1 & 2: Identify the BCKDH substrate with median Bertz's complexity
    bckdh_substrates = {
        'KIV': 'CC(C)C(=O)C(=O)O',    # a-Ketoisovalerate
        'KMV': 'CCC(C)C(=O)C(=O)O',   # a-Keto-b-methylvalerate
        'KIC': 'CC(C)CC(=O)C(=O)O',  # a-Ketoisocaproate
    }

    substrate_complexities = []
    for name, smiles in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smiles)
        bertz_ct = GraphDescriptors.BertzCT(mol)
        substrate_complexities.append({'name': name, 'smiles': smiles, 'bertz_ct': bertz_ct})

    # Sort by complexity to find the median
    substrate_complexities.sort(key=lambda x: x['bertz_ct'])
    median_substrate = substrate_complexities[1]

    # Step 3: Calculate Balaban J for the median substrate
    # RDKit can return 0 for acyclic molecules unless force_compute is True
    median_mol = Chem.MolFromSmiles(median_substrate['smiles'])
    target_balaban_j = GraphDescriptors.BalabanJ(median_mol, forceCompute=True)

    # Step 4: Find the ribonucleoside with the closest Balaban J index
    nucleosides = {
        'Adenosine': 'C1=NC2=C(N1C3C(C(C(O3)CO)O)O)N=CN=C2N',
        'Guanosine': 'C1=NC2=C(N1C3C(C(C(O3)CO)O)O)NC(=NC2=O)N',
        'Cytidine':  'C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)O',
        'Uridine':   'C1=CN(C(=O)NC1=O)C2C(C(C(O2)CO)O)O',
    }

    closest_nucleoside = None
    min_diff = float('inf')

    for name, smiles in nucleosides.items():
        mol = Chem.MolFromSmiles(smiles)
        balaban_j = GraphDescriptors.BalabanJ(mol)
        diff = abs(balaban_j - target_balaban_j)
        if diff < min_diff:
            min_diff = diff
            closest_nucleoside = {'name': name, 'smiles': smiles}
    
    final_molecule_smiles = closest_nucleoside['smiles']

    # Step 5: Calculate the final ratio for the chosen molecule (Cytidine)
    
    # 5a. Create the H-included graph
    mol = Chem.MolFromSmiles(final_molecule_smiles)
    h_mol = Chem.AddHs(mol)
    
    G = nx.Graph()
    for atom in h_mol.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in h_mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # 5b. Calculate Zagreb(1) Index (M1)
    m1_index = sum(d**2 for n, d in G.degree())

    # 5c. Calculate Hosoya Z Index
    # Increase recursion limit for deep calculations
    sys.setrecursionlimit(G.number_of_nodes() * G.number_of_edges())
    
    memo_z = {}
    def calculate_hosoya_z(graph):
        # For disconnected graphs, Z(G) is the product of Z of its components
        if nx.number_connected_components(graph) > 1:
            component_hashes = []
            for c in nx.connected_components(graph):
                subg = graph.subgraph(c)
                # Sort edges within each component to create a canonical frozenset
                component_hashes.append(frozenset(tuple(sorted(e)) for e in subg.edges()))
            
            # Sort the component hashes to create a canonical key for the set of components
            graph_id = tuple(sorted(component_hashes, key=lambda fs: len(fs)))

            if graph_id in memo_z:
                return memo_z[graph_id]
            
            z = 1
            for c in nx.connected_components(graph):
                subgraph_copy = graph.subgraph(c).copy()
                z *= calculate_hosoya_z(subgraph_copy)
            memo_z[graph_id] = z
            return z
        
        # For a single connected component, use the sorted edge tuple as the key
        graph_id = tuple(sorted(tuple(sorted(e)) for e in graph.edges()))
        if graph_id in memo_z:
            return memo_z[graph_id]

        if graph.number_of_edges() == 0:
            return 1

        # Recursive step: Z(G) = Z(G-e) + Z(G-{u,v})
        # Pick a canonical edge to improve memoization
        u, v = graph_id[0]
        
        G_minus_e = graph.copy()
        G_minus_e.remove_edge(u, v)
        term1 = calculate_hosoya_z(G_minus_e)

        G_minus_uv = graph.copy()
        G_minus_uv.remove_nodes_from([u, v])
        term2 = calculate_hosoya_z(G_minus_uv)
        
        result = term1 + term2
        memo_z[graph_id] = result
        return result

    z_index = calculate_hosoya_z(G)

    # 5d. Compute the final ratio
    final_ratio = (2 * z_index) / m1_index

    print(f"The chosen substance is {closest_nucleoside['name']}.")
    print("The final calculation is (2 * H-included Hosoya Z) / Zagreb(1):")
    print(f"(2 * {z_index}) / {m1_index} = {final_ratio}")

# Run the solver
solve_task()
<<<170.48863636363635>>>