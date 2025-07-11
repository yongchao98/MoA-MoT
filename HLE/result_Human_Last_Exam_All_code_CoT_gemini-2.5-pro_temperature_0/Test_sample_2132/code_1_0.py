import networkx as nx
from itertools import combinations

def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone (cyclopentanecarboxylic acid).
    """
    # 1. Define the molecular graph (non-hydrogen atoms)
    # The product is cyclopentanecarboxylic acid.
    # Nodes:
    # 1-5: Carbons in the cyclopentane ring
    # 6:   Carbon of the carboxylic acid group
    # 7:   Carbonyl oxygen
    # 8:   Hydroxyl oxygen
    G = nx.Graph()
    G.add_edges_from([
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 1),  # Cyclopentane ring
        (1, 6),                                  # Link to carboxyl group
        (6, 7), (6, 8)                           # Carboxyl group bonds
    ])

    # 2. Determine N_atoms and N_bonds
    n_atoms = G.number_of_nodes()
    n_bonds = G.number_of_edges()

    # 3. Calculate N_paths (total number of simple paths between all pairs)
    n_paths = 0
    # Iterate over all unique pairs of nodes
    for u, v in combinations(G.nodes(), 2):
        # Find all simple paths between the pair and add the count to the total
        num_paths_for_pair = len(list(nx.all_simple_paths(G, source=u, target=v)))
        n_paths += num_paths_for_pair

    # 4. Calculate Böttcher Molecular Complexity (BMC)
    if n_paths == 0:
        bmc = 0
    else:
        bmc = (n_bonds * (n_atoms ** 2)) / n_paths

    # 5. Print the results in the required equation format
    print("The product is cyclopentanecarboxylic acid.")
    print(f"Number of non-hydrogen atoms (N_atoms): {n_atoms}")
    print(f"Number of non-hydrogen bonds (N_bonds): {n_bonds}")
    print(f"Total number of paths (N_paths): {n_paths}")
    print("\nBöttcher Molecular Complexity Calculation:")
    # The final equation as requested, showing each number
    print(f"{n_bonds} * ({n_atoms}^2) / {n_paths} = {bmc}")

if __name__ == '__main__':
    calculate_bottcher_complexity()