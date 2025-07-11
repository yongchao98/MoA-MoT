import math

def calculate_and_print_symmetry_factors():
    """
    Calculates and prints the symmetry factors for second-order vacuum bubble
    diagrams in phi^4 theory and their sum.
    """

    print("Calculating the symmetry factors for the two connected second-order vacuum bubble diagrams in phi^4 theory.")
    print("-" * 80)

    # --- Diagram 1: The "Melon" Diagram ---
    # This diagram consists of two vertices connected by four identical propagators.
    # Its symmetry factor S1 is the product of:
    # 1. Permutations of the 4 identical propagators: 4!
    # 2. Permutations of the 2 identical vertices: 2!
    
    s1_propagator_perms = math.factorial(4)
    s1_vertex_perms = math.factorial(2)
    s1 = s1_propagator_perms * s1_vertex_perms

    print("1. The 'Melon' Diagram:")
    print(f"   - Symmetry from permuting 4 identical propagators connecting the vertices: 4! = {s1_propagator_perms}")
    print(f"   - Symmetry from swapping the 2 identical vertices: 2! = {s1_vertex_perms}")
    print(f"   - Total Symmetry Factor S1 = {s1_propagator_perms} * {s1_vertex_perms} = {s1}")
    print("")

    # --- Diagram 2: The "Figure-8" Diagram ---
    # This diagram has two vertices, each with a self-loop, and two propagators
    # connecting the vertices. Its symmetry factor S2 is the product of:
    # 1. Permutations of the 2 identical connecting propagators: 2!
    # 2. A factor of 2 for each of the 2 self-loops: 2 * 2
    # 3. Permutations of the 2 identical vertices: 2!

    s2_propagator_perms = math.factorial(2)
    s2_loop_factor = 2 * 2
    s2_vertex_perms = math.factorial(2)
    s2 = s2_propagator_perms * s2_loop_factor * s2_vertex_perms

    print("2. The 'Figure-8' Diagram:")
    print(f"   - Symmetry from permuting 2 identical propagators connecting the vertices: 2! = {s2_propagator_perms}")
    print(f"   - Symmetry from 2 self-loops (a factor of 2 for each): 2 * 2 = {s2_loop_factor}")
    print(f"   - Symmetry from swapping the 2 identical vertices: 2! = {s2_vertex_perms}")
    print(f"   - Total Symmetry Factor S2 = {s2_propagator_perms} * {s2_loop_factor} * {s2_vertex_perms} = {s2}")
    print("")

    # --- Summation ---
    total_sum = s1 + s2

    print("-" * 80)
    print("The sum of the symmetry factors is:")
    print(f"Total = S1 + S2 = {s1} + {s2} = {total_sum}")


if __name__ == "__main__":
    calculate_and_print_symmetry_factors()
    total = math.factorial(4) * math.factorial(2) + math.factorial(2) * (2*2) * math.factorial(2)
    # The final answer is wrapped in <<<>>>
    print(f"\n<<< {total} >>>")
