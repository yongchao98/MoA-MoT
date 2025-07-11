import math

def calculate_symmetry_factors():
    """
    Calculates the symmetry factors for the second-order connected vacuum bubble
    diagrams in phi^4 theory and their sum.
    """

    print("Calculating symmetry factors for second-order vacuum bubble diagrams in φ⁴ theory.")
    print("--------------------------------------------------------------------------\n")

    # --- Diagram 1: The "Saturn" Diagram ---
    # This diagram has two vertices, each with a self-loop, connected by two propagators.
    print("1. The 'Saturn' Diagram:")
    print("   This diagram has two vertices connected by two propagators, with a self-loop on each vertex.")
    
    # Symmetry factors for the Saturn diagram:
    # 1. Permuting the two identical vertices: factor of 2
    s1_vtx_perm = 2
    # 2. Permuting the two identical propagators connecting the vertices: factor of 2
    s1_inter_prop_perm = 2
    # 3. Flipping the orientation of the loop on the first vertex (swapping the two half-lines): factor of 2
    s1_loop1_perm = 2
    # 4. Flipping the orientation of the loop on the second vertex: factor of 2
    s1_loop2_perm = 2
    
    s1_total = s1_vtx_perm * s1_inter_prop_perm * s1_loop1_perm * s1_loop2_perm
    
    print(f"   - Symmetry from swapping the two vertices: {s1_vtx_perm}")
    print(f"   - Symmetry from swapping the two propagators between vertices: {s1_inter_prop_perm}")
    print(f"   - Symmetry from flipping the loop on the first vertex: {s1_loop1_perm}")
    print(f"   - Symmetry from flipping the loop on the second vertex: {s1_loop2_perm}")
    print(f"   Total symmetry factor S1 = {s1_vtx_perm} * {s1_inter_prop_perm} * {s1_loop1_perm} * {s1_loop2_perm} = {s1_total}\n")

    # --- Diagram 2: The "Melon" Diagram ---
    # This diagram has two vertices connected by four propagators.
    print("2. The 'Melon' Diagram:")
    print("   This diagram has two vertices connected by four propagators.")
    
    # Symmetry factors for the Melon diagram:
    # 1. Permuting the two identical vertices: factor of 2
    s2_vtx_perm = 2
    # 2. Permuting the four identical propagators connecting the two vertices: factor of 4!
    s2_inter_prop_perm = math.factorial(4)
    
    s2_total = s2_vtx_perm * s2_inter_prop_perm

    print(f"   - Symmetry from swapping the two vertices: {s2_vtx_perm}")
    print(f"   - Symmetry from permuting the four propagators between vertices: {s2_inter_prop_perm}")
    print(f"   Total symmetry factor S2 = {s2_vtx_perm} * {s2_inter_prop_perm} = {s2_total}\n")
    
    # --- Sum of Symmetry Factors ---
    total_sum = s1_total + s2_total
    print("--------------------------------------------------------------------------")
    print("Sum of all symmetry factors:")
    print(f"Total = S1 + S2 = {s1_total} + {s2_total} = {total_sum}")

calculate_symmetry_factors()
<<<64>>>