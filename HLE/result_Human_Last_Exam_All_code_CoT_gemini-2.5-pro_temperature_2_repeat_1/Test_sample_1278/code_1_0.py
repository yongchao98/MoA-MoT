import math

def solve_pirate_standoff():
    """
    Calculates the maximum number of Mexican standoffs (cliques of size >= 3)
    for 9 pirates and 16 pairs at gunpoint.
    """
    
    # --- Problem Setup ---
    num_pirates = 9
    num_gunpoint_pairs = 16

    # --- Optimal Structure: K_5 + 4 satellite vertices ---
    # We found the optimal structure is a K_5 (5 pirates in a full clique)
    # with the remaining 4 pirates attached.
    k5_pirates = 5
    satellite_pirates = 4
    
    k5_edges = math.comb(k5_pirates, 2)
    remaining_edges = num_gunpoint_pairs - k5_edges
    
    # --- Part 1: Standoffs within the K_5 base ---
    # Calculate standoffs (cliques) of size 3, 4, and 5 within the K_5.
    standoffs_k3_base = math.comb(k5_pirates, 3)
    standoffs_k4_base = math.comb(k5_pirates, 4)
    standoffs_k5_base = math.comb(k5_pirates, 5)
    
    total_base_standoffs = standoffs_k3_base + standoffs_k4_base + standoffs_k5_base

    # --- Part 2: New standoffs from satellite pirates ---
    # The remaining 6 edges are distributed among the 4 satellite pirates to maximize
    # new standoffs. The optimal distribution of connections to the K_5 is (3, 1, 1, 1).
    connections = [3, 1, 1, 1]
    
    # The pirate with 3 connections to the K_5 forms new standoffs.
    # A standoff is formed by this pirate plus 2 or more pirates from the K_5 group it's connected to.
    d1 = connections[0]
    new_standoffs_k3 = math.comb(d1, 2) # Forms K3s: (satellite, v_i, v_j)
    new_standoffs_k4 = math.comb(d1, 3) # Forms K4s: (satellite, v_i, v_j, v_k)
    
    # The other pirates with 1 connection each do not form any new standoffs of size >= 3.
    # C(1, 2) = 0, C(1, 3) = 0 etc.
    d2 = connections[1]
    other_new_standoffs = 3 * (math.comb(d2, 2) + math.comb(d2, 3)) # This will be 0

    total_new_standoffs = new_standoffs_k3 + new_standoffs_k4 + other_new_standoffs
    
    # --- Final Calculation ---
    total_standoffs = total_base_standoffs + total_new_standoffs

    # --- Print the explanation and final equation ---
    print("The problem is to find the maximum number of cliques of size >= 3 in a connected graph with 9 vertices and 16 edges.")
    print("The optimal structure is a K_5 clique connected to the 4 other vertices.")
    print("\n1. Standoffs (cliques) within the base K_5:")
    print(f"   - 3-pirate standoffs (K3): C(5,3) = {standoffs_k3_base}")
    print(f"   - 4-pirate standoffs (K4): C(5,4) = {standoffs_k4_base}")
    print(f"   - 5-pirate standoffs (K5): C(5,5) = {standoffs_k5_base}")
    
    print("\n2. New standoffs formed by attaching the remaining 4 pirates:")
    print("   To maximize new standoffs, one pirate connects to 3 pirates in the K_5.")
    print(f"   - New 3-pirate standoffs (K3): C(3,2) = {new_standoffs_k3}")
    print(f"   - New 4-pirate standoffs (K4): C(3,3) = {new_standoffs_k4}")
    
    print("\n--- Total Calculation ---")
    print(f"Total Standoffs = (Standoffs from K_5) + (New Standoffs)")
    print(f"Total = ({standoffs_k3_base} + {standoffs_k4_base} + {standoffs_k5_base}) + ({new_standoffs_k3} + {new_standoffs_k4})")
    print(f"Total = ({total_base_standoffs}) + ({total_new_standoffs}) = {total_standoffs}")

solve_pirate_standoff()