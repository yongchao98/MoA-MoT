def calculate_min_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given adjacent and non-adjacent swap rules.
    """
    total_displacement_magnitude = 0
    N = 100

    # The equation for the final answer will be built as a string
    equation_str = []

    for k in range(1, N + 1):
        # Initial position of element k (0-indexed)
        initial_pos = k - 1
        # Final position of element k (0-indexed)
        final_pos = N - k
        
        # Total displacement needed for element k
        displacement = final_pos - initial_pos
        
        # The part of the displacement that must be covered by paid (+/- 1) moves
        # is the remainder when divided by 5. We find the shortest path to a multiple of 5.
        rem = displacement % 5
        
        # For example, a remainder of 4 is 1 move away from a multiple of 5 (4 = 5 - 1).
        # A remainder of 3 is 2 moves away (3 = 5 - 2).
        if rem > 5 / 2:
            paid_displacement = 5 - rem
        else:
            paid_displacement = rem
            
        total_displacement_magnitude += paid_displacement
        if k <= 5: # Show calculation for the first 5 elements as an example
            equation_str.append(f"{paid_displacement}")

    # Each adjacent swap contributes 2 to the total magnitude of displacement (+1 for one element, -1 for the other).
    # So, the total number of swaps is half the total magnitude.
    min_moves = total_displacement_magnitude / 2

    # Print the explanation and the final equation
    print("The final answer is the sum of the minimum paid moves for each element, all divided by 2.")
    print("The individual minimum paid moves for the first 5 elements are:")
    print(" + ".join(equation_str) + " + ...")
    print("Summing these for all 100 elements gives a total displacement magnitude of {}.".format(total_displacement_magnitude))
    print("Since each swap accounts for a displacement magnitude of 2, we divide by 2.")
    print(f"{total_displacement_magnitude} / 2 = {int(min_moves)}")
    print("\nFinal Answer:")
    print(int(min_moves))

calculate_min_moves()