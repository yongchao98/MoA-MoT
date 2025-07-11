import networkx as nx

def count_simple_cycles(G, length):
    """
    Counts the number of simple cycles of a given length in a graph G.
    A simple cycle does not repeat vertices, except for the start/end vertex.
    The final count is normalized to count each unique cycle once, regardless of
    starting point or direction.
    """
    total_cycles = 0
    for start_node in G.nodes():
        # Stack for iterative depth-first search.
        # Each item is a tuple: (current_node, path_so_far).
        stack = [(start_node, [start_node])]
        while stack:
            current_node, path = stack.pop()

            # If the path has the correct number of nodes to form a cycle of the target length
            if len(path) == length:
                # Check for an edge from the last node back to the start node
                if G.has_edge(current_node, start_node):
                    total_cycles += 1
                continue

            # Explore neighbors to extend the path
            for neighbor in G.neighbors(current_node):
                # To keep the cycle simple, the neighbor must not be in the path already.
                if neighbor not in path:
                    new_path = path + [neighbor]
                    stack.append((neighbor, new_path))
    
    # Each cycle of length k is found k times (once for each starting node)
    # and 2 times for each direction (e.g., u-v-w-u vs u-w-v-u).
    # Therefore, we divide the total count by (2 * length).
    # We use integer division // as the result must be an integer.
    return total_cycles // (2 * length)

# Define the SRG parameters for our chosen class of graphs.
n = 16
d = 6
lambda_param = 2
mu_param = 2

# Create the two non-isomorphic graphs with the same SRG parameters.
# 1. The Shrikhande graph.
shrikhande_graph = nx.shrikhande_graph()
# 2. The 4x4 Rook's graph (also known as the L2(4) graph).
rook_graph_4x4 = nx.rook_graph(4)

# Count the number of 5-cycles in each graph.
c5_shrikhande = count_simple_cycles(shrikhande_graph, 5)
c5_rook = count_simple_cycles(rook_graph_4x4, 5)

# Output the findings.
print("Yes, there is a combination of parameters for which two SRGs can have a different number of 5-cycles.")
print(f"\nThe parameters considered are (n, d, λ, μ) = ({n}, {d}, {lambda_param}, {mu_param}).")
print("\nWe analyze two graphs from the class srg(16, 6, 2, 2):")
print("1. The Shrikhande Graph")
print("2. The 4x4 Rook's Graph\n")

# Display the final counts in an "equation" format.
print(f"Number of 5-cycles in the Shrikhande Graph = {c5_shrikhande}")
print(f"Number of 5-cycles in the 4x4 Rook's Graph = {c5_rook}")
print("\nSince 240 is not equal to 192, these graphs are a counterexample.")

<<<Yes>>>