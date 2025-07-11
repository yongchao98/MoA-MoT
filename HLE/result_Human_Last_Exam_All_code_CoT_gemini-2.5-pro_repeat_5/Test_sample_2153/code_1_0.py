import numpy as np
import networkx as nx
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

def solve_chemistry_problem():
    """
    Solves the multi-step cheminformatics problem.
    """

    # --- Step 1 & 2: Identify Reference Substrate and its Balaban J Index ---

    # Define the BCKDH complex substrates
    bckdh_substrates = {
        "KIC": "CC(C)CC(=O)C(=O)O",      # a-Ketoisocaproate
        "KMV": "CCC(C)C(=O)C(=O)O",      # a-Keto-b-methylvalerate
        "KIV": "CC(C)C(=O)C(=O)O"       # a-Ketoisovalerate
    }

    # Calculate Bertz complexity for each substrate
    bertz_values = {}
    for name, smiles in bckdh_substrates.items():
        mol = Chem.MolFromSmiles(smiles)
        bertz_values[name] = GraphDescriptors.BertzCT(mol)

    # Find the median Bertz complexity
    median_bertz = np.median(list(bertz_values.values()))

    # Identify substrates with the median complexity and calculate their Balaban J indices
    reference_balaban_j_values = []
    for name, bertz in bertz_values.items():
        if np.isclose(bertz, median_bertz):
            mol = Chem.MolFromSmiles(bckdh_substrates[name])
            reference_balaban_j_values.append(GraphDescriptors.BalabanJ(mol))

    # The reference Balaban J is the average of the values for the median-complexity substrates
    ref_balaban_j = np.mean(reference_balaban_j_values)

    # --- Step 3: Identify the Target Molecule ---

    # Define the substances from Carrell's 2018 synthesis
    carrell_substances = {
        "Uridine": "C1=CN(C(=O)NC1=O)C2C(C(C(O2)CO)O)O",
        "Cytidine": "C1=CN(C(=O)N=C1N)C2C(C(C(O2)CO)O)O",
        "T-Uridine": "C1=CN(C(=O)NC1=O)C2C(C(CO)O2)O",
        "T-Cytidine": "C1=CN(C(=O)N=C1N)C2C(C(CO)O2)O"
    }

    # Find the substance with the Balaban J index closest to the reference value
    min_diff = float('inf')
    target_molecule_name = ""
    target_molecule_smiles = ""

    for name, smiles in carrell_substances.items():
        mol = Chem.MolFromSmiles(smiles)
        balaban_j = GraphDescriptors.BalabanJ(mol)
        diff = abs(balaban_j - ref_balaban_j)
        if diff < min_diff:
            min_diff = diff
            target_molecule_name = name
            target_molecule_smiles = smiles
            
    # --- Step 4: Calculate the Final Ratio for the Target Molecule ---

    # Helper function to convert an RDKit molecule to a NetworkX graph
    def mol_to_nx(mol):
        G = nx.Graph()
        for atom in mol.GetAtoms():
            G.add_node(atom.GetIdx())
        for bond in mol.GetBonds():
            G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        return G

    # Memoization cache for the Hosoya Z calculation
    hosoya_memo = {}
    
    # Function to calculate the Hosoya Z index using a recursive formula
    def calculate_hosoya_z(G):
        graph_id = frozenset(G.edges())
        if graph_id in hosoya_memo:
            return hosoya_memo[graph_id]

        if G.number_of_edges() == 0:
            return 1

        # Recurrence relation: Z(G) = Z(G-e) + Z(G-V(e))
        edge = list(G.edges())[0]
        u, v = edge
        
        G_minus_e = G.copy()
        G_minus_e.remove_edge(u, v)
        
        G_minus_uv = G.copy()
        G_minus_uv.remove_nodes_from([u, v])
        
        result = calculate_hosoya_z(G_minus_e) + calculate_hosoya_z(G_minus_uv)
        hosoya_memo[graph_id] = result
        return result

    # Create the molecule object for the target substance and add hydrogens
    mol = Chem.MolFromSmiles(target_molecule_smiles)
    mol_h = Chem.AddHs(mol)

    # Calculate the Zagreb(1) index
    zagreb_m1 = sum(atom.GetDegree()**2 for atom in mol_h.GetAtoms())

    # Convert to a NetworkX graph and calculate the Hosoya Z index
    graph_h = mol_to_nx(mol_h)
    hosoya_z_value = calculate_hosoya_z(graph_h)

    # Calculate the final ratio
    final_ratio = (2 * hosoya_z_value) / zagreb_m1
    
    # --- Step 5: Output the Result ---
    
    print(f"The identified target molecule is: {target_molecule_name}")
    print(f"The final equation is: (2 * {hosoya_z_value}) / {zagreb_m1} = {final_ratio}")

# Execute the main function
solve_chemistry_problem()