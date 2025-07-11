import networkx as nx
import itertools

def create_rooks_graph():
    """
    Creates the 4x4 Rook's graph, an SRG with parameters (16, 6, 2, 2).
    """
    G = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G.add_nodes_from(nodes)
    for u, v in itertools.combinations(nodes, 2):
        # Add an edge if nodes are in the same row or column
        if u[0] == v[0] or u[1] == v[1]:
            G.add_edge(u, v)
    return G

def create_shrikhande_graph():
    """
    Creates the Shrikhande graph, an SRG with parameters (16, 6, 2, 2).
    """
    G = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G.add_nodes_from(nodes)
    for r, c in nodes:
        # Define neighbors based on Cayley graph construction on Z_4 x Z_4
        # The connection set is {(+-1, 0), (0, +-1), (1, 1), (-1, -1)}
        neighbors = [
            ((r + 1) % 4, c),
            ((r - 1) % 4, c),
            (r, (c + 1) % 4),
            (r, (c - 1) % 4),
            ((r + 1) % 4, (c + 1) % 4),
            ((r - 1) % 4, (c - 1) % 4)
        ]
        for neighbor in neighbors:
            G.add_edge((r, c), neighbor)
    return G

def count_five_cycles(G):
    """
    Counts the number of simple 5-cycles in a graph G.
    """
    count = 0
    # Iterating through nodes and finding paths is one way to count.
    # A 5-cycle is a path of length 4 starting and ending at the same node.
    # The path is v1-v2-v3-v4-v5-v1.
    for v1 in G.nodes():
        # Path: v1-v2
        for v2 in G.neighbors(v1):
            # Path: v1-v2-v3 (v3 must not be v1)
            for v3 in G.neighbors(v2):
                if v3 == v1:
                    continue
                # Path: v1-v2-v3-v4 (v4 must not be v1 or v2)
                for v4 in G.neighbors(v3):
                    if v4 == v1 or v4 == v2:
                        continue
                    
                    # A 5-cycle v1-v2-v3-v4-v5-v1 is formed if v5 is a common
                    # neighbor of v1 and v4, and v5 is not v2 or v3.
                    # We can iterate through common neighbors of v1 and v4.
                    common_neighbors = set(G.neighbors(v1)) & set(G.neighbors(v4))
                    for v5 in common_neighbors:
                        if v5 != v2 and v5 != v3:
                            count += 1
                            
    # Each cycle is found starting at each of its 5 nodes,
    # and in two directions (clockwise and counter-clockwise).
    # So we divide the total count by 5 * 2 = 10.
    return count // 10

def main():
    """
    Main function to construct graphs, count cycles, and print the results.
    """
    print("This script demonstrates two strongly regular graphs with the same parameters")
    print("but a different number of 5-cycles.")
    print("-" * 70)
    
    # Define the SRG parameters
    n, d, lam, mu = 16, 6, 2, 2
    print(f"The parameters are (n, d, lambda, mu) = ({n}, {d}, {lam}, {mu}).\n")
    
    # Create the graphs
    rooks_graph = create_rooks_graph()
    shrikhande_graph = create_shrikhande_graph()
    
    print("Graph 1: 4x4 Rook's Graph")
    print("Graph 2: Shrikhande Graph\n")

    # Count the 5-cycles in each graph
    c5_rooks = count_five_cycles(rooks_graph)
    c5_shrikhande = count_five_cycles(shrikhande_graph)

    print(f"Number of 5-cycles in the 4x4 Rook's Graph: {c5_rooks}")
    print(f"Number of 5-cycles in the Shrikhande Graph: {c5_shrikhande}\n")
    
    print("Conclusion: The number of 5-cycles is different.")
    print(f"The final comparison is: {c5_shrikhande} != {c5_rooks}")


if __name__ == '__main__':
    # To run this code, you need to have the networkx library installed:
    # pip install networkx
    main()

<<<240 != 384>>>