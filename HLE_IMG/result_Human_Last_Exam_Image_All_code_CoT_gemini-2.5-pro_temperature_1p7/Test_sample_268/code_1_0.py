import math
from fractions import Fraction

def calculate_clustering_coefficient():
    """
    Calculates the clustering coefficient of the central node in the 64-node graph.
    """
    # Step 1: Calculate the degree 'k' of the central node.
    # The central node is connected to its 3 original blue neighbors.
    blue_neighbors = 3
    
    # It's also connected to all nodes in the three 4-node green triads.
    green_neighbors = 3 * 4
    
    # It's also connected to all nodes in the three 16-node red structures.
    red_neighbors = 3 * 16
    
    k = blue_neighbors + green_neighbors + red_neighbors
    
    # Step 2: Calculate the number of triangles 'T' involving the central node.
    # T is the number of edges between the neighbors of the central node.
    # These edges are the ones internal to the green and red sub-structures.
    
    # Edges within the 3 green triads. Each is a 4-node star graph with 3 edges.
    edges_in_green_substructures = 3 * 3
    
    # Edges within the 3 red structures. Each red structure is a copy of the
    # 16-node (blue+green) structure.
    # We need to count the edges in one 16-node structure.
    # - Edges in its "blue" core triad: 3
    # - Edges in its three "green" peripheral triads: 3 * 3 = 9
    # - Edges connecting the central node to all 12 "green" nodes: 12
    edges_in_one_16_node_structure = 3 + (3 * 3) + (3 * 4)
    edges_in_red_substructures = 3 * edges_in_one_16_node_structure
    
    T = edges_in_green_substructures + edges_in_red_substructures
    
    # Step 3: Calculate the clustering coefficient using the formula.
    # C = 2 * T / (k * (k - 1))
    
    numerator = 2 * T
    denominator = k * (k - 1)
    
    # Use Fraction for an exact rational result
    c_fraction = Fraction(numerator, denominator)
    c_float = float(numerator) / denominator
    
    # Step 4: Print the results step-by-step.
    print("Step 1: Calculate the degree (k) of the central node.")
    print(f"The central node is connected to {blue_neighbors} blue, {green_neighbors} green, and {red_neighbors} red nodes.")
    print(f"k = {blue_neighbors} + {green_neighbors} + {red_neighbors} = {k}")
    print("\nStep 2: Calculate the number of triangles (T).")
    print("T is the number of edges between the neighbors of the central node.")
    print(f"Edges in green substructures = {edges_in_green_substructures}")
    print(f"Edges in red substructures = 3 * (edges in one 16-node structure) = {edges_in_red_substructures}")
    print(f"T = {edges_in_green_substructures} + {edges_in_red_substructures} = {T}")
    print("\nStep 3: Calculate the clustering coefficient (C).")
    print("The formula is C = (2 * T) / (k * (k - 1))")
    print(f"C = (2 * {T}) / ({k} * ({k} - 1))")
    print(f"C = {numerator} / ({k} * {k-1})")
    print(f"C = {numerator} / {denominator}")
    print(f"C = {c_fraction.numerator} / {c_fraction.denominator} (simplified fraction)")
    print(f"C \u2248 {c_float:.5f}")


calculate_clustering_coefficient()
print("\n<<<9/217>>>")