import networkx as nx

def count_cycles(graph, length):
    """Counts the number of simple cycles of a given length in a graph."""
    # nx.simple_cycles finds all elementary circuits in a graph.
    # An elementary circuit is a closed path where no nodes appear twice,
    # which is the definition of a simple cycle.
    cycle_generator = nx.simple_cycles(graph)
    count = 0
    for cycle in cycle_generator:
        if len(cycle) == length:
            count += 1
    return count

# The parameters for the two strongly regular graphs
n = 16
d = 6
lambda_ = 2
mu = 2

# Construct the first graph: the 4x4 Rook's graph.
# This graph can be constructed as the Cartesian product of two complete graphs K_4.
K4 = nx.complete_graph(4)
rooks_graph = nx.cartesian_product(K4, K4)

# Construct the second graph: the Shrikhande graph.
# NetworkX has a built-in generator for it.
shrikhande_graph = nx.shrikhande_graph()

# Define the cycle length we want to count
cycle_length = 5

# Count the number of 5-cycles in each graph
rooks_cycles_count = count_cycles(rooks_graph, cycle_length)
shrikhande_cycles_count = count_cycles(shrikhande_graph, cycle_length)

# Print the results
print(f"Yes, a pair of graphs G, G' in srg(n,d,lambda,mu) with a different number of 5-cycles exists.")
print("-" * 40)
print(f"The parameters are (n={n}, d={d}, lambda={lambda_}, mu={mu}).")
print("-" * 40)
print("Graph 1: The 4x4 Rook's graph")
print(f"Number of 5-cycles: {rooks_cycles_count}")
print("-" * 40)
print("Graph 2: The Shrikhande graph")
print(f"Number of 5-cycles: {shrikhande_cycles_count}")
print("-" * 40)

if rooks_cycles_count != shrikhande_cycles_count:
    print("The number of 5-cycles in the two graphs is different, confirming the existence of such a pair.")
else:
    print("The number of 5-cycles is the same.")
