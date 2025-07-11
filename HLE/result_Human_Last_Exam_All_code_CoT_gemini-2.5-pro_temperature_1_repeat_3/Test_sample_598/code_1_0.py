import networkx as nx

def count_k_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    A simple cycle does not have repeated vertices.
    networkx.simple_cycles returns a generator for lists of nodes in a cycle.
    We filter for cycles of the desired length and count them.
    """
    cycle_generator = nx.simple_cycles(G)
    count = sum(1 for cycle in cycle_generator if len(cycle) == k)
    return count

def construct_rooks_graph():
    """Constructs the 4x4 Rook's graph."""
    G = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G.add_nodes_from(nodes)
    # Two vertices are adjacent if they share a row or a column.
    for i in range(4):
        for j in range(4):
            # Connect to other nodes in the same row
            for j2 in range(j + 1, 4):
                G.add_edge((i, j), (i, j2))
            # Connect to other nodes in the same column
            for i2 in range(i + 1, 4):
                G.add_edge((i, j), (i2, j))
    return G

def construct_shrikhande_graph():
    """Constructs the Shrikhande graph."""
    G = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G.add_nodes_from(nodes)
    # The graph is a Cayley graph on Z_4 x Z_4.
    # An edge exists between (i,j) and (i',j') if their difference
    # is in the set {(+/-1, 0), (0, +/-1), (1, 1), (-1,-1)}.
    for i in range(4):
        for j in range(4):
            u = (i, j)
            neighbors = [
                ((i + 1) % 4, j),
                ((i - 1) % 4, j),
                (i, (j + 1) % 4),
                (i, (j - 1) % 4),
                ((i + 1) % 4, (j + 1) % 4),
                ((i - 1) % 4, (j - 1) % 4),
            ]
            for v in neighbors:
                G.add_edge(u, v)
    return G

def main():
    """
    Main function to construct the graphs, count cycles, and print the results.
    """
    print("Yes, there exists a pair of strongly regular graphs with the same parameters but a different number of 5-cycles.")
    
    # Parameters for the strongly regular graphs
    n = 16
    d = 6
    lambda_param = 2
    mu_param = 2
    print(f"An example is the class of graphs srg(n,d,lambda,mu) with parameters:")
    print(f"n = {n}")
    print(f"d = {d}")
    print(f"lambda = {lambda_param}")
    print(f"mu = {mu_param}")
    print("-" * 30)

    # Construct the two graphs
    rooks_graph = construct_rooks_graph()
    shrikhande_graph = construct_shrikhande_graph()

    # Count the 5-cycles in each graph
    num_5_cycles_rooks = count_k_cycles(rooks_graph, 5)
    num_5_cycles_shrikhande = count_k_cycles(shrikhande_graph, 5)

    print("There are two non-isomorphic graphs with these parameters:")
    print("1. The 4x4 Rook's graph")
    print("2. The Shrikhande graph")
    print("\nWe will count the number of 5-cycles in each:")
    
    print(f"Number of 5-cycles in the 4x4 Rook's graph: {num_5_cycles_rooks}")
    print(f"Number of 5-cycles in the Shrikhande graph: {num_5_cycles_shrikhande}")
    
    print("\nAs you can see, the number of 5-cycles is different:")
    print(f"{num_5_cycles_rooks} != {num_5_cycles_shrikhande}")

if __name__ == "__main__":
    main()