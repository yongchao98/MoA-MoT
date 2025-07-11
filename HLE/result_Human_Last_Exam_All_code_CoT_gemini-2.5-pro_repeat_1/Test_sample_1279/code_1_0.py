import math

def calculate_symmetry_factors():
    """
    Calculates the symmetry factors for all second-order vacuum bubble diagrams
    in phi^4 theory and sums them.
    """
    print("Calculating symmetry factors for second-order vacuum bubble diagrams in phi^4 theory.")
    print("-" * 70)

    # --- Diagram A: The "Saturn" Diagram ---
    # This connected diagram has two vertices connected by two propagators,
    # and each vertex also has a self-loop.
    # The symmetries are:
    # 1. Swapping the two propagators connecting the vertices (factor of 2).
    # 2. Flipping the loop at the first vertex (swapping its legs, factor of 2).
    # 3. Flipping the loop at the second vertex (swapping its legs, factor of 2).
    # Note: By standard convention, we do not count the symmetry of swapping the two vertices,
    # as they are indistinguishable and integrated over the whole spacetime volume.
    s_a_swap_lines = 2
    s_a_flip_loop1 = 2
    s_a_flip_loop2 = 2
    s_a = s_a_swap_lines * s_a_flip_loop1 * s_a_flip_loop2
    print("Diagram A (Saturn Diagram):")
    print(f"  - Symmetries: Swap 2 connecting lines ({s_a_swap_lines}) * Flip loop 1 ({s_a_flip_loop1}) * Flip loop 2 ({s_a_flip_loop2})")
    print(f"  - Symmetry Factor S_A = {s_a}\n")

    # --- Diagram B: The "Melon" Diagram ---
    # This connected diagram has two vertices connected by four propagators.
    # The symmetry is the number of ways to permute these four identical propagators.
    s_b_permute_lines = math.factorial(4)
    s_b = s_b_permute_lines
    print("Diagram B (Melon Diagram):")
    print(f"  - Symmetries: Permutations of the 4 connecting lines (4!)")
    print(f"  - Symmetry Factor S_B = {s_b}\n")

    # --- Diagram C: The Disconnected Diagram ---
    # This diagram consists of two identical, disconnected first-order vacuum bubbles (figure-eights).
    # First, we calculate the symmetry factor of a single figure-eight diagram (S_1).
    # A figure-eight has one vertex and two identical loops.
    # Symmetries of a figure-eight:
    # 1. Swapping the two identical loops (factor of 2).
    # 2. Flipping the first loop (factor of 2).
    # 3. Flipping the second loop (factor of 2).
    s1_swap_loops = 2
    s1_flip_loop1 = 2
    s1_flip_loop2 = 2
    s_1 = s1_swap_loops * s1_flip_loop1 * s1_flip_loop2
    print("Component: First-Order Figure-Eight Diagram")
    print(f"  - Symmetries: Swap 2 loops ({s1_swap_loops}) * Flip loop 1 ({s1_flip_loop1}) * Flip loop 2 ({s1_flip_loop2})")
    print(f"  - Symmetry Factor S_1 = {s_1}\n")

    # The symmetry factor for Diagram C is S_1 * S_1 * 2! (for swapping the two identical components).
    s_c_swap_components = math.factorial(2)
    s_c = s_1 * s_1 * s_c_swap_components
    print("Diagram C (Two Disconnected Figure-Eights):")
    print(f"  - Symmetries: (Symm of 1st) * (Symm of 2nd) * (Swap components)")
    print(f"  - Symmetry Factor S_C = {s_1} * {s_1} * {s_c_swap_components} = {s_c}\n")

    # --- Sum of all symmetry factors ---
    total_s = s_a + s_b + s_c
    print("-" * 70)
    print("Sum of all second-order symmetry factors:")
    print(f"Total = S_A + S_B + S_C")
    # The final print statement fulfills the requirement to output each number in the final equation.
    print(f"Total = {s_a} + {s_b} + {s_c} = {total_s}")
    print("-" * 70)

    # Returning the final answer as a string to be captured
    return total_s

# Execute the function and print the final result in the desired format
final_answer = calculate_symmetry_factors()
print(f'<<<{final_answer}>>>')
