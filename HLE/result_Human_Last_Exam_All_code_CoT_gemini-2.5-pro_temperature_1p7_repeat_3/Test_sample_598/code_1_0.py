import networkx as nx
from itertools import combinations

def count_5_cycles(G):
    """
    Counts the number of simple 5-cycles in a graph.

    It works by checking every combination of 5 vertices. If the subgraph
    induced by these 5 vertices has 5 edges and is 2-regular, it's a 5-cycle.
    This method is very reliable and efficient enough for small graphs.

    Args:
        G (networkx.Graph): The graph to analyze.

    Returns:
        int: The number of 5-cycles in the graph.
    """
    count = 0
    # Iterate through all unique combinations of 5 vertices in the graph
    for nodes in combinations(G.nodes(), 5):
        # Create a subgraph from the 5 chosen vertices
        sub = G.subgraph(nodes)
        
        # A 5-cycle must have exactly 5 vertices and 5 edges.
        # Additionally, every vertex in a cycle must have a degree of 2.
        if sub.number_of_edges() == 5:
            is_cycle = True
            for node in sub.nodes():
                if sub.degree(node) != 2:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    return count

# Define the parameters of the strongly regular graphs
n = 16
d = 6
lambda_param = 2
mu_param = 2

print(f"Finding two graphs in srg(n={n}, d={d}, lambda={lambda_param}, mu={mu_param}) with a different number of 5-cycles.")
print("-" * 20)

# Create the first graph: the L2(4) graph (4x4 rook's graph)
# networkx.rook_graph(4) creates this graph.
G1_name = "L2(4) Rook's Graph"
G1 = nx.rook_graph(4)

# Create the second graph: the Shrikhande graph
# networkx has a built-in generator for it.
G2_name = "Shrikhande Graph"
G2 = nx.shrikhande_graph()

# Count the 5-cycles in each graph
num_5_cycles_G1 = count_5_cycles(G1)
num_5_cycles_G2 = count_5_cycles(G2)

# Print the results
print(f"Graph 1: {G1_name}")
print(f"Number of 5-cycles: {num_5_cycles_G1}")
print("-" * 20)
print(f"Graph 2: {G2_name}")
print(f"Number of 5-cycles: {num_5_cycles_G2}")
print("-" * 20)

# Conclude by showing the numbers are different
print("Comparing the number of 5-cycles:")
if num_5_cycles_G1 != num_5_cycles_G2:
    print(f"The number of 5-cycles is different: {num_5_cycles_G1} != {num_5_cycles_G2}")
else:
    print(f"The number of 5-cycles is the same: {num_5_cycles_G1} == {num_5_cycles_G2}")
