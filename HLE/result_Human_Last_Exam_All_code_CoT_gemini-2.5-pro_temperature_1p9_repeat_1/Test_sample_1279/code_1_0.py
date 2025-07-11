import math

def calculate_symmetry_factors():
    """
    Calculates the symmetry factors for the two second-order vacuum bubble diagrams
    in phi^4 theory and prints the sum.
    """

    # --- Diagram 1: The "Dumbbell" diagram ---
    # Symmetries:
    # 1. Swapping the two identical vertices (factor of 2)
    # 2. Swapping the two identical propagators connecting the vertices (factor of 2)
    # 3. Symmetry of the loop at the first vertex (factor of 2)
    # 4. Symmetry of the loop at the second vertex (factor of 2)
    s1_vtx_perms = 2
    s1_prop_perms = 2
    s1_loop1_sym = 2
    s1_loop2_sym = 2
    s1 = s1_vtx_perms * s1_prop_perms * s1_loop1_sym * s1_loop2_sym

    print("Calculation for the 'Dumbbell' diagram (S_1):")
    print(f"S_1 = (vertex permutations) * (propagator permutations) * (loop 1 sym) * (loop 2 sym)")
    print(f"S_1 = {s1_vtx_perms} * {s1_prop_perms} * {s1_loop1_sym} * {s1_loop2_sym} = {s1}")
    print("-" * 30)

    # --- Diagram 2: The "Setting Sun" diagram ---
    # Symmetries:
    # 1. Swapping the two identical vertices (factor of 2)
    # 2. Permuting the four identical propagators connecting the vertices (factor of 4!)
    s2_vtx_perms = 2
    s2_prop_perms = math.factorial(4)
    s2 = s2_vtx_perms * s2_prop_perms
    
    print("Calculation for the 'Setting Sun' diagram (S_2):")
    print(f"S_2 = (vertex permutations) * (propagator permutations)")
    print(f"S_2 = {s2_vtx_perms} * {s2_prop_perms} = {s2}")
    print("-" * 30)

    # --- Total Sum ---
    total_s = s1 + s2
    print("Total sum of symmetry factors (S_total):")
    print(f"S_total = S_1 + S_2")
    print(f"S_total = {s1} + {s2} = {total_s}")

calculate_symmetry_factors()
<<<64>>>