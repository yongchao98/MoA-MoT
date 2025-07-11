import math

def calculate_symmetry_factors():
    """
    Calculates the symmetry factors for the two second-order vacuum bubble
    diagrams in phi^4 theory and prints their sum.
    """

    print("Calculating the symmetry factors for second-order vacuum bubble diagrams in phi^4 theory.")
    print("-" * 70)

    # --- Diagram 1: "Doughnut" Diagram ---
    # Two vertices, connected by four propagators.
    print("Diagram 1: 'Doughnut' Diagram")
    
    # Factor for permuting the two identical vertices
    d1_vertex_permutations = math.factorial(2)
    
    # Factor for permuting the four identical propagators between the vertices
    d1_line_permutations = math.factorial(4)
    
    # Calculate the symmetry factor S1
    S1 = d1_vertex_permutations * d1_line_permutations
    
    print(f"Symmetry factor S1 = (vertex permutations) * (line permutations)")
    print(f"S1 = {d1_vertex_permutations} * {d1_line_permutations} = {S1}")
    print("-" * 70)

    # --- Diagram 2: "Double-Scoop" Diagram ---
    # Two vertices, two propagators connecting them, and one loop on each vertex.
    print("Diagram 2: 'Double-Scoop' Diagram")
    
    # Factor for permuting the two identical vertices
    d2_vertex_permutations = math.factorial(2)
    
    # Factor for permuting the two identical propagators between the vertices
    d2_connecting_line_permutations = math.factorial(2)
    
    # Factor for reversing the loop on the first vertex
    d2_loop1_reversal = 2
    
    # Factor for reversing the loop on the second vertex
    d2_loop2_reversal = 2

    # Calculate the symmetry factor S2
    S2 = d2_vertex_permutations * d2_connecting_line_permutations * d2_loop1_reversal * d2_loop2_reversal

    print(f"Symmetry factor S2 = (vertex perms) * (connecting line perms) * (loop 1 reversal) * (loop 2 reversal)")
    print(f"S2 = {d2_vertex_permutations} * {d2_connecting_line_permutations} * {d2_loop1_reversal} * {d2_loop2_reversal} = {S2}")
    print("-" * 70)

    # --- Total Sum ---
    total_sum = S1 + S2
    print("Sum of all symmetry factors:")
    # The user requested that the final equation shows each number.
    print(f"Total Sum = S1 + S2 = {S1} + {S2} = {total_sum}")

calculate_symmetry_factors()
<<<64>>>