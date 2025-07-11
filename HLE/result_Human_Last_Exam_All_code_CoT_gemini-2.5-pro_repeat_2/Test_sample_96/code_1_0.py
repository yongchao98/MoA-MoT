import networkx as nx
from itertools import combinations

def solve_torsion_element_count():
    """
    This function solves the problem by reducing it to a graph-theoretic counting problem.
    It counts the number of induced subgraphs of the E8 Dynkin diagram that are
    isomorphic to the D6 Dynkin diagram.
    """
    
    # Step 1: Explain the theoretical reduction of the problem.
    print("The problem asks for the number of certain torsion elements in A/Z.")
    print("Based on the theory of Artin-Tits groups, this count is equivalent to the number of parabolic subgroups of type D6 within the E8 group.")
    print("This, in turn, is a graph theory problem: counting the number of induced subgraphs of the E8 Dynkin diagram that are isomorphic to the D6 Dynkin diagram.")
    print("We will now solve this counting problem using Python.")
    print("-" * 20)

    # Step 2: Define the E8 and D6 Dynkin diagrams as graphs using networkx.
    
    # The E8 Dynkin diagram:
    #       3
    #       |
    # 1--2--4--5--6--7--8
    E8 = nx.Graph()
    E8_edges = [(1, 2), (2, 3), (2, 4), (4, 5), (5, 6), (6, 7), (7, 8)]
    E8.add_nodes_from(range(1, 9))
    E8.add_edges_from(E8_edges)

    # The D6 Dynkin diagram:
    # 1--2--3--4--5
    #          |
    #          6
    D6 = nx.Graph()
    D6_edges = [(1, 2), (2, 3), (3, 4), (4, 5), (4, 6)]
    D6.add_nodes_from(range(1, 7))
    D6.add_edges_from(D6_edges)

    n_d6 = len(D6.nodes())
    e8_nodes = list(E8.nodes())
    
    # Step 3: Iterate through all subsets of E8's vertices of the correct size.
    count = 0
    found_subgraphs_nodes = []
    
    for nodes_subset in combinations(e8_nodes, n_d6):
        subgraph = E8.subgraph(nodes_subset)
        
        # An irreducible parabolic subgroup must correspond to a connected subgraph.
        if nx.is_connected(subgraph):
            # Check if the subgraph is isomorphic to D6.
            if nx.is_isomorphic(subgraph, D6):
                count += 1
                found_subgraphs_nodes.append(sorted(list(nodes_subset)))

    # Step 4: Output the result of the count.
    print(f"Searching for induced subgraphs of E8 isomorphic to D6...")
    print(f"Number of D6 parabolic subgroups found: {count}")
    
    if found_subgraphs_nodes:
        print(f"The vertex set(s) for the D6 subdiagram(s) are: {found_subgraphs_nodes}")

    print("-" * 20)
    print("Each such subgroup corresponds to a unique torsion element satisfying the conditions.")
    print(f"Therefore, the total number of such torsion elements is {count}.")


solve_torsion_element_count()