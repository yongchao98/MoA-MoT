def solve_reversal_moves():
    """
    Calculates the minimum number of moves to reverse a sequence of 100 elements
    with the given swap operations.
    """
    N = 100  # Number of elements in the sequence
    D = 5    # The distance for a non-adjacent (free) swap

    # There are D bins of positions, indexed 0 to D-1.
    # The boundaries are between bin k and k+1, for k = 0, 1, 2, 3.
    num_boundaries = D - 1
    moves_per_boundary = []
    total_moves = 0

    # Calculate the number of moves required across each boundary.
    for k in range(num_boundaries):
        crossings = 0
        # Iterate through each element/position from 1 to N.
        for i in range(1, N + 1):
            # Initial position is i.
            # The sequence is reversed, so the final position is N+1-i.
            final_pos = N + 1 - i

            # Determine the initial and final bin for the element at position i.
            # Bins are 0-indexed.
            initial_bin = (i - 1) % D
            final_bin = (final_pos - 1) % D

            # An element must cross the boundary 'k' if it starts in a bin
            # up to k and ends in a bin beyond k.
            if initial_bin <= k and final_bin > k:
                crossings += 1
        
        moves_per_boundary.append(crossings)
        total_moves += crossings

    # Format the output string to show the calculation.
    equation_parts = [str(m) for m in moves_per_boundary]
    equation_str = " + ".join(equation_parts)
    
    print(f"The minimum number of moves is the sum of elements crossing each of the {num_boundaries} boundaries:")
    print(f"{equation_str} = {total_moves}")

solve_reversal_moves()