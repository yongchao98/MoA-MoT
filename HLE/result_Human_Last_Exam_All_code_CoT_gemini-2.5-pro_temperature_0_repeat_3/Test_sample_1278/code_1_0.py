import math

def count_cycles_in_Kk(k):
    """Counts the total number of simple cycles of length >= 3 in a complete graph K_k."""
    total_cycles = 0
    # A cycle of length l is determined by choosing l vertices from k,
    # and arranging them in a circle.
    # The number of ways to choose l vertices is C(k, l).
    # The number of ways to arrange l vertices in a cycle is (l-1)!/2.
    for l in range(3, k + 1):
        # Using math.comb for combinations C(k,l)
        # Using math.factorial for (l-1)!
        cycles_of_length_l = math.comb(k, l) * math.factorial(l - 1) // 2
        total_cycles += cycles_of_length_l
    return total_cycles

def count_simple_paths_in_Kk(k):
    """Counts the number of simple paths between two specific vertices in K_k."""
    # Let the two vertices be v1 and v2.
    # The number of other vertices available for the path is k-2.
    # A simple path can have from 1 to k-1 edges (length 2 to k).
    
    # Path of length 1 (direct edge)
    num_paths = 1
    
    # Paths of length > 1. A path of length p has p-1 intermediate vertices.
    # These p-1 vertices must be chosen from the k-2 available vertices.
    # The number of ways to choose and order them is P(k-2, p-1).
    for p_minus_1 in range(1, k - 1): # p_minus_1 is the number of intermediate vertices
        num_paths += math.perm(k - 2, p_minus_1)
        
    return num_paths

def solve():
    """
    Calculates the maximum number of cycles based on the optimal graph structure.
    The optimal structure is a K_5 connected to a P_4 by 3 edges.
    """
    # Part 1: Cycles within the K_5 clique
    k5_vertices = 5
    cycles_in_k5 = count_cycles_in_Kk(k5_vertices)
    print(f"The graph structure that maximizes cycles is a K_5 clique connected to the remaining 4 vertices.")
    print(f"The 4 remaining vertices form a P_4 path, and there are 3 edges connecting the K_5 and the P_4.")
    print(f"\nStep 1: Calculate cycles entirely within the K_5 clique.")
    print(f"Number of cycles in a K_5: {cycles_in_k5}")

    # Part 2: Cycles within the subgraph of the other 4 vertices
    # This subgraph is a P_4 (a path on 4 vertices), which has 0 cycles.
    cycles_in_p4 = 0
    print(f"\nStep 2: Calculate cycles within the remaining part of the graph.")
    print(f"The subgraph on the other 4 vertices is a path (P_4), which has {cycles_in_p4} cycles.")

    # Part 3: Calculate mixed cycles that span the two parts
    # These cycles are formed by pairs of the 3 connecting edges.
    num_connecting_edges = 3
    num_pairs_of_connections = math.comb(num_connecting_edges, 2)
    
    # Number of paths between two vertices in K_5
    paths_in_k5 = count_simple_paths_in_Kk(k5_vertices)
    
    # Number of paths between connection points in the P_4.
    # The 3 connecting edges are attached to distinct vertices in the P_4,
    # so there is always exactly 1 simple path between any two of them.
    paths_in_p4 = 1
    
    mixed_cycles = num_pairs_of_connections * paths_in_k5 * paths_in_p4
    
    print(f"\nStep 3: Calculate mixed cycles using the connecting edges.")
    print(f"There are {num_connecting_edges} edges connecting the K_5 and P_4, forming {num_pairs_of_connections} pairs.")
    print(f"Number of simple paths between any two vertices in the K_5 is {paths_in_k5}.")
    print(f"Number of simple paths between the connection points in the P_4 is {paths_in_p4}.")
    print(f"Number of mixed cycles = {num_pairs_of_connections} * {paths_in_k5} * {paths_in_p4} = {mixed_cycles}")

    # Final calculation
    total_cycles = cycles_in_k5 + cycles_in_p4 + mixed_cycles
    print(f"\nFinal Calculation:")
    print(f"Total cycles = (cycles in K_5) + (cycles in P_4) + (mixed cycles)")
    print(f"Total cycles = {cycles_in_k5} + {cycles_in_p4} + {mixed_cycles} = {total_cycles}")

solve()