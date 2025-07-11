import math

def calculate_phi4_vacuum_bubbles():
    """
    Calculates the symmetry factors for the two second-order vacuum bubble
    diagrams in phi^4 theory and prints the sum.
    """
    print("In phi^4 theory, there are two distinct connected second-order vacuum bubble diagrams.")
    print("Let's calculate the symmetry factor (S) for each.\n")

    # --- Diagram 1: Figure-eight Diagram ---
    print("1. The 'Figure-eight' Diagram:")
    print("   This diagram has two vertices. Each vertex has a self-loop,")
    print("   and the two vertices are connected by two propagators.\n")
    
    # The symmetry factor is the product of the orders of its independent symmetries.
    # Symmetry 1: Swapping the two vertices. They are identical.
    vertex_swap_factor = 2
    # Symmetry 2: Swapping the two propagators connecting the vertices. They are identical.
    propagator_swap_factor = 2
    # Symmetry 3: Flipping the orientation of the first self-loop propagator.
    loop1_flip_factor = 2
    # Symmetry 4: Flipping the orientation of the second self-loop propagator.
    loop2_flip_factor = 2
    
    s1 = vertex_swap_factor * propagator_swap_factor * loop1_flip_factor * loop2_flip_factor
    
    print("   Its symmetry factor S1 is calculated as follows:")
    print(f"   - Factor from swapping the 2 identical vertices: {vertex_swap_factor}")
    print(f"   - Factor from swapping the 2 identical propagators between vertices: {propagator_swap_factor}")
    print(f"   - Factor from reversing the first self-loop: {loop1_flip_factor}")
    print(f"   - Factor from reversing the second self-loop: {loop2_flip_factor}")
    print(f"   S1 = {vertex_swap_factor} * {propagator_swap_factor} * {loop1_flip_factor} * {loop2_flip_factor} = {s1}\n")
    
    # --- Diagram 2: Setting-sun Diagram ---
    print("2. The 'Setting-sun' Diagram:")
    print("   This diagram has two vertices connected by four identical propagators.\n")

    # Symmetry 1: Swapping the two vertices. They are identical.
    vertex_swap_factor_2 = 2
    # Symmetry 2: Permuting the four identical propagators connecting the vertices.
    num_propagators = 4
    propagator_perm_factor = math.factorial(num_propagators)

    s2 = vertex_swap_factor_2 * propagator_perm_factor
    
    print("   Its symmetry factor S2 is calculated as follows:")
    print(f"   - Factor from swapping the 2 identical vertices: {vertex_swap_factor_2}")
    print(f"   - Factor from permuting the {num_propagators} identical propagators: {num_propagators}! = {propagator_perm_factor}")
    print(f"   S2 = {vertex_swap_factor_2} * {propagator_perm_factor} = {s2}\n")

    # --- Summation ---
    total_sum = s1 + s2
    print("--------------------------------------------------")
    print("Finally, we sum the symmetry factors of all diagrams:")
    print(f"Total Sum = S1 + S2")
    print(f"Total Sum = {s1} + {s2} = {total_sum}")

if __name__ == '__main__':
    calculate_phi4_vacuum_bubbles()