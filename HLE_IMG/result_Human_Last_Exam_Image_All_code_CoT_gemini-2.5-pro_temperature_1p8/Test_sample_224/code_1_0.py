# This script requires the 'networkx' library.
# You can install it by running: pip install networkx

import networkx as nx
import itertools
from collections import defaultdict

def find_and_verify_petersen_cdc():
    """
    This function demonstrates finding a cycle double cover for the Petersen graph.
    """
    print("Finding a Cycle Double Cover for the Petersen Graph...")

    # We use the standard Petersen graph from the networkx library.
    # It has 10 vertices and 15 edges.
    G = nx.petersen_graph()
    all_graph_edges = list(G.edges())
    num_edges = len(all_graph_edges)
    
    # For a cycle double cover, the sum of cycle lengths must be 2 * num_edges.
    # Equation: 2 * 15 = 30
    print(f"The sum of cycle lengths must be 2 * {num_edges} = {2 * num_edges}")

    # We will search for a cover made of five 6-cycles.
    # Equation: 5 * 6 = 30
    print(f"We will search for a cover made of five 6-cycles: 5 * 6 = {5 * 6}\n")

    # Find all simple cycles of length 6 in the graph.
    all_simple_cycles = nx.simple_cycles(G)
    cycles_len_6 = [cycle for cycle in all_simple_cycles if len(cycle) == 6]
    
    print(f"Found {len(cycles_len_6)} cycles of length 6 in the Petersen graph.")

    # There are C(10, 5) = 252 combinations of five 6-cycles. We test each one.
    found_cdc = None
    for combo in itertools.combinations(cycles_len_6, 5):
        edge_counts = defaultdict(int)
        
        for cycle in combo:
            # Convert the cycle's node path into edges
            for i in range(len(cycle)):
                u, v = cycle[i], cycle[(i + 1) % len(cycle)]
                # Normalize edge representation (e.g., (1, 0) becomes (0, 1))
                edge = tuple(sorted((u, v)))
                edge_counts[edge] += 1

        # Check if this combination forms a valid CDC
        if all(count == 2 for count in edge_counts.values()) and len(edge_counts) == num_edges:
            found_cdc = combo
            break
            
    if found_cdc:
        print("\nSuccessfully found a valid Cycle Double Cover:")
        # Note: Vertex labels (0-9) are from networkx, not necessarily the image.
        for i, cycle in enumerate(found_cdc):
            print(f"  Cycle {i+1}: {cycle}")
    else:
        print("Could not find a CDC of this type with this method.")

find_and_verify_petersen_cdc()

# The main question is: "How many cycle double covers does the Petersen Graph have up to isomorphism?"
# This is a known result from graph theory.
# There are 6 such non-isomorphic cycle double covers.
print("\nThe final answer to the number of non-isomorphic cycle double covers is:")
print(6)
<<<6>>>