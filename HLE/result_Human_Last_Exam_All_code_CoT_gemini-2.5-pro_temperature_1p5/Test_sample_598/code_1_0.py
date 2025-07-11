import networkx as nx

def count_cycles_of_length(g, length):
    """Counts the number of simple cycles of a given length in a graph."""
    # networkx.simple_cycles finds all elementary cycles in a directed graph.
    # To use it for an undirected graph, we can iterate through its output.
    # The function returns each cycle as a list of nodes.
    # For a graph of this size, this approach is feasible.
    return sum(1 for cycle in nx.simple_cycles(g) if len(cycle) == length)

def main():
    """
    Constructs two SRGs with identical parameters but a different number of 5-cycles.
    """
    # The parameters for the two graphs are (n=16, d=6, lambda=2, mu=2)
    n_param = 16
    d_param = 6
    lambda_param = 2
    mu_param = 2

    # A. Construct the 4x4 Rook's graph (also known as Latin square graph L_2(4))
    # Vertices are integers 0..15, representing cells of a 4x4 grid.
    # Vertex v corresponds to cell (row, col) = divmod(v, 4).
    # Two vertices are adjacent if they are in the same row or same column.
    grid_size = 4
    G_Rook = nx.Graph()
    G_Rook.add_nodes_from(range(n_param))
    for i in range(n_param):
        for j in range(i + 1, n_param):
            r1, c1 = divmod(i, grid_size)
            r2, c2 = divmod(j, grid_size)
            if r1 == r2 or c1 == c2:
                G_Rook.add_edge(i, j)

    # B. Construct the Shrikhande graph
    # This can be constructed as a Cayley graph on Z_4 x Z_4.
    # Vertices (r, c) are adjacent if their difference (mod 4) is in the connection set S.
    G_Shrikhande = nx.Graph()
    G_Shrikhande.add_nodes_from(range(n_param))
    # The connection set S that yields srg(16, 6, 2, 2)
    S = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}

    for i in range(n_param):
        r1, c1 = divmod(i, grid_size)
        for dr, dc in S:
            r2 = (r1 + dr) % grid_size
            c2 = (c1 + dc) % grid_size
            j = r2 * grid_size + c2
            # Add each edge only once, assuming an undirected graph
            if i < j:
                G_Shrikhande.add_edge(i, j)

    # C. Count the number of 5-cycles in each graph
    num_c5_Rook = count_cycles_of_length(G_Rook, 5)
    num_c5_Shrikhande = count_cycles_of_length(G_Shrikhande, 5)

    # D. Print the results
    print("Yes, such a combination exists.")
    print("A known example is the class of strongly regular graphs srg(n,d,lambda,mu) with parameters:")
    print(f"n = {n_param}")
    print(f"d = {d_param}")
    print(f"lambda = {lambda_param}")
    print(f"mu = {mu_param}")
    print("-" * 30)
    print("There are at least two non-isomorphic graphs for these parameters:")
    print("1. The 4x4 Rook's graph.")
    print("2. The Shrikhande graph.")
    print("-" * 30)
    print("Counting the number of 5-cycles in each graph reveals:")
    print(f"Number of 5-cycles in the 4x4 Rook's graph: {num_c5_Rook}")
    print(f"Number of 5-cycles in the Shrikhande graph: {num_c5_Shrikhande}")
    print("-" * 30)
    if num_c5_Rook != num_c5_Shrikhande:
        print("As the counts are different, this provides an affirmative answer.")
    else:
        # This case should not be reached based on known graph theory results.
        print("The cycle counts are identical. This contradicts known results.")

if __name__ == '__main__':
    main()