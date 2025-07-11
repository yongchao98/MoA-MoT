import math

def calculate_symmetry_factors():
    """
    Calculates the symmetry factors for the two connected second-order
    vacuum bubble diagrams in phi^4 scalar field theory and sums them.
    """
    print("Calculating symmetry factors for second-order vacuum bubbles in phi^4 theory.\n")

    # --- Diagram 1: The "Sunset" Diagram ---
    # This diagram has two vertices and four propagators connecting them.
    # The symmetry factor S is the product of:
    # 1. Permutations of identical vertices. There are 2, so 2! = 2.
    # 2. Permutations of identical propagators. There are 4 identical propagators
    #    connecting the two vertices, so 4! = 24.
    
    perms_vertices_sunset = math.factorial(2)
    perms_propagators_sunset = math.factorial(4)
    s_sunset = perms_vertices_sunset * perms_propagators_sunset
    
    print("--- Sunset Diagram ---")
    print(f"Permutations of 2 identical vertices: {perms_vertices_sunset}")
    print(f"Permutations of 4 identical propagators: {perms_propagators_sunset}")
    print(f"Symmetry factor for the Sunset diagram = {perms_vertices_sunset} * {perms_propagators_sunset} = {s_sunset}\n")

    # --- Diagram 2: The "Figure-8" Diagram ---
    # This diagram has two vertices, two propagators forming loops (one at each vertex),
    # and two propagators connecting the vertices.
    # The symmetry factor S is the product of:
    # 1. Permutations of identical vertex-loop structures. The two vertex-loop complexes are identical and can be swapped: factor of 2.
    # 2. Permutations of identical propagators connecting the vertices. There are 2, so factor of 2.
    # 3. Permutations of propagator ends within a loop. The two lines forming the loop at the first vertex can be swapped: factor of 2.
    # 4. Permutations of propagator ends within the second loop. Similarly, the two lines for the second loop can be swapped: factor of 2.

    perms_vertices_fig8 = 2
    perms_propagators_fig8 = 2
    perms_loop1_ends_fig8 = 2
    perms_loop2_ends_fig8 = 2
    s_fig8 = perms_vertices_fig8 * perms_propagators_fig8 * perms_loop1_ends_fig8 * perms_loop2_ends_fig8
    
    print("--- Figure-8 Diagram ---")
    print(f"Permutation of the 2 identical vertex-loop structures: {perms_vertices_fig8}")
    print(f"Permutations of the 2 identical propagators between vertices: {perms_propagators_fig8}")
    print(f"Swapping the 2 lines in the first loop: {perms_loop1_ends_fig8}")
    print(f"Swapping the 2 lines in the second loop: {perms_loop2_ends_fig8}")
    print(f"Symmetry factor for the Figure-8 diagram = {perms_vertices_fig8} * {perms_propagators_fig8} * {perms_loop1_ends_fig8} * {perms_loop2_ends_fig8} = {s_fig8}\n")
    
    # --- Summing the factors ---
    total_sum = s_fig8 + s_sunset
    print("--- Total Sum ---")
    print("The sum of the symmetry factors is:")
    print(f"{s_fig8} (Figure-8) + {s_sunset} (Sunset) = {total_sum}")

    return total_sum

if __name__ == '__main__':
    final_answer = calculate_symmetry_factors()
    print(f"\n<<<>>>\n{final_answer}\n<<<>>>")
