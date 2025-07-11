import networkx as nx
import numpy as np

def calculate_total_variation(G, x0):
    """
    Calculates the total variation of a signal x0 on a graph G.
    This is defined as the sum of the absolute differences of the signal
    values across all edges.
    """
    # x1 is the signal on edges, defined as |x0_u - x0_v|
    x1 = {}
    tv_sum = 0.0
    
    print("Calculating Total Variation: sum(|x0_u - x0_v| for each edge {u,v})")
    sum_terms = []

    for u, v in G.edges():
        # Get signal values at the nodes of the edge
        x0_u = x0[u]
        x0_v = x0[v]
        
        # Calculate the value for the edge signal x1
        x1_e = abs(x0_u - x0_v)
        x1[(u, v)] = x1_e
        
        # Add to the total variation sum
        tv_sum += x1_e
        sum_terms.append(f"{x1_e:.1f}")

    # Print the full equation for clarity
    equation_str = " + ".join(sum_terms)
    print(f"Total Variation = {equation_str} = {tv_sum:.1f}")
    
    # Based on the problem's premises, this sum must be 0.
    # This implies that x0 is constant over each connected component.
    # Our example x0 reflects this, so the total variation will be 0.

# 1. Create a graph G with nodes V and edges E.
# Let's create a graph with a cycle (a square) and a bridge (an edge attached to it).
G = nx.Graph()
V = [0, 1, 2, 3, 4]
G.add_nodes_from(V)
E = [(0, 1), (1, 2), (2, 3), (3, 0), (2, 4)] # Edges (0,1,2,3) form a cycle, (2,4) is a bridge
G.add_edges_from(E)

# 2. Define a vertex signal x0.
# According to our derivation, x1 must be the zero vector, which means
# |x0_u - x0_v| = 0 for all edges. This happens if x0 is constant.
# The graph is connected, so x0 must be constant everywhere.
c = 5.0
x0 = {node: c for node in G.nodes()}

print(f"Graph has nodes: {G.nodes()}")
print(f"Graph has edges: {G.edges()}")
print(f"Vertex signal x0: {x0}")
print("-" * 20)

# 3. Calculate and print the Total Variation
calculate_total_variation(G, x0)
