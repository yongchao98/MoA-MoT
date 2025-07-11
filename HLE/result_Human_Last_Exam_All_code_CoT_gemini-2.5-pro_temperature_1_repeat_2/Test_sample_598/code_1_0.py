import networkx as nx
from collections import deque

def count_simple_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    A simple cycle is one where vertices are not repeated, except for the start/end vertex.
    The algorithm uses a modified breadth-first search (BFS) to find paths of length k-1,
    and then checks if the last node is connected to the first.
    """
    # This algorithm counts each cycle 2*k times (k starting points, 2 directions).
    # We will divide the final count by 2*k.
    total_cycles = 0
    
    for start_node in G.nodes():
        # queue stores tuples of (current_node, path_list)
        queue = deque([(start_node, [start_node])])
        
        while queue:
            current_node, path = queue.popleft()
            
            if len(path) == k:
                # Path of length k-1 found. Check if it forms a cycle.
                if G.has_edge(path[-1], path[0]):
                    total_cycles += 1
                continue

            for neighbor in G.neighbors(current_node):
                # To ensure a simple path, the neighbor must not be in the current path.
                # Exception: if it's the last step to form a cycle, we allow it.
                # But our path length check handles this already.
                if neighbor not in path:
                    new_path = path + [neighbor]
                    queue.append((neighbor, new_path))
                        
    return total_cycles // (2 * k)

def construct_rooks_graph():
    """Constructs the 4x4 Rook's graph, an srg(16,6,2,2)."""
    G = nx.Graph()
    n = 16
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            # Vertices are cells (row, col) of a 4x4 grid.
            # Map vertex index to grid coordinates.
            r1, c1 = i // 4, i % 4
            r2, c2 = j // 4, j % 4
            # Two vertices are adjacent if they are in the same row or column.
            if r1 == r2 or c1 == c2:
                G.add_edge(i, j)
    return G

def construct_shrikhande_graph():
    """Constructs the Shrikhande graph, an srg(16,6,2,2)."""
    G = nx.Graph()
    n = 16
    G.add_nodes_from(range(n))
    # Connection set for the Cayley graph construction of the Shrikhande graph.
    # The graph is Cay(Z_4 x Z_4, S).
    S = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    for i in range(n):
        for j in range(i + 1, n):
            r1, c1 = i // 4, i % 4
            r2, c2 = j // 4, j % 4
            # Two vertices are adjacent if their difference (mod 4) is in the set S.
            diff = (abs(r1 - r2) % 4, abs(c1 - c2) % 4)
            if diff in S:
                G.add_edge(i, j)
    return G
    
def main():
    """
    Main function to construct graphs, count cycles, and print the result.
    """
    print("Investigating strongly regular graphs with parameters (n,d,lambda,mu) = (16,6,2,2)")
    print("-" * 70)
    
    # Construct the two graphs
    rooks_graph = construct_rooks_graph()
    shrikhande_graph = construct_shrikhande_graph()
    
    # Count 5-cycles in each graph
    print("Counting 5-cycles... (this may take a moment)")
    num_5_cycles_rooks = count_simple_cycles(rooks_graph, 5)
    num_5_cycles_shrikhande = count_simple_cycles(shrikhande_graph, 5)
    
    print("\n--- Results ---")
    
    # Output results for the Rook's graph
    n_r, d_r = rooks_graph.number_of_nodes(), list(rooks_graph.degree())[0][1]
    print(f"Graph 1: 4x4 Rook's Graph")
    print(f"Parameters: n={n_r}, d={d_r}, lambda=2, mu=2")
    print(f"Number of 5-cycles = {num_5_cycles_rooks}")

    print("-" * 20)

    # Output results for the Shrikhande graph
    n_s, d_s = shrikhande_graph.number_of_nodes(), list(shrikhande_graph.degree())[0][1]
    print(f"Graph 2: Shrikhande Graph")
    print(f"Parameters: n={n_s}, d={d_s}, lambda=2, mu=2")
    print(f"Number of 5-cycles = {num_5_cycles_shrikhande}")

    print("\n--- Conclusion ---")
    if num_5_cycles_rooks != num_5_cycles_shrikhande:
        print("Yes, there exist two SRGs with the same parameters but a different number of 5-cycles.")
        print(f"The Rook's graph has {num_5_cycles_rooks} 5-cycles, while the Shrikhande graph has {num_5_cycles_shrikhande} 5-cycles.")
    else:
        print("Based on this example, the number of 5-cycles is the same. My analysis might be flawed.")

if __name__ == "__main__":
    main()
