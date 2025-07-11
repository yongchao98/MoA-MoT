def solve_temporal_fixed_points():
    """
    Solves for the sum of all positive integer temporal fixed points less than 100.
    """
    # The temporal horizon constant, as given in the problem.
    h_T = 48

    # List to store the positive integer fixed points less than 100.
    fixed_points = []

    # A temporal fixed point 'k' is a value where the results of k ⊗ k
    # are the same from both time-flow directions.
    # Forward result: |k + k| = |2k|
    # Backward result: |k - k| = |0| = 0
    # The temporal horizon h(T) implies calculations are modular.
    # So, the condition is: |2k| is congruent to 0 modulo h(T).
    # |2k| % h_T == 0
    
    # We are looking for positive integers k < 100.
    for k in range(1, 100):
        # Since k is positive, |2k| is 2*k.
        # The condition is (2 * k) % h_T == 0.
        if (2 * k) % h_T == 0:
            fixed_points.append(k)

    # Calculate the sum of the found fixed points.
    total_sum = sum(fixed_points)

    # Print the results in the required format.
    print(f"The positive integer temporal fixed points less than 100 are: {fixed_points}")
    
    # Create the equation string.
    equation_str = " + ".join(map(str, fixed_points))
    print(f"The sum is: {equation_str} = {total_sum}")

solve_temporal_fixed_points()
<<<240>>>