def solve_snail_puzzle():
    """
    Calculates the maximal distance a snail can travel in 7 minutes
    based on the given observational constraints.

    The problem can be modeled with a recurrence relation D(t) = D(t-1) + D(t-2),
    where D(t) is the maximum distance in t minutes.
    This is the same recurrence as the Fibonacci sequence.

    Base cases:
    D(1) = 1 meter.
    D(2) = 2 meters.
    """
    
    # Initialize base cases for t=1 and t=2
    # d_n_minus_2 corresponds to D(t-2)
    # d_n_minus_1 corresponds to D(t-1)
    d_n_minus_2 = 1  # Corresponds to D(1)
    d_n_minus_1 = 2  # Corresponds to D(2)

    print("Calculating the maximal distance for t=1 to t=7 minutes:")
    print("D(1) = 1")
    print("D(2) = 2")

    # We already have results for t=1 and t=2, so we loop from t=3 to t=7.
    for t in range(3, 8):
        # Apply the recurrence relation D(t) = D(t-1) + D(t-2)
        d_n = d_n_minus_1 + d_n_minus_2
        
        # Print the equation for the current step
        print(f"D({t}) = D({t-1}) + D({t-2}) = {d_n_minus_1} + {d_n_minus_2} = {d_n}")
        
        # Update the values for the next iteration
        d_n_minus_2 = d_n_minus_1
        d_n_minus_1 = d_n

    final_distance = d_n_minus_1
    print(f"\nThe maximal distance the snail could have advanced in 7 minutes is {final_distance} meters.")

solve_snail_puzzle()