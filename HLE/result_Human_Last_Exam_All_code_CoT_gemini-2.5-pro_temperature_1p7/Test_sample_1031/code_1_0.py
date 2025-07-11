def solve_disease_spread():
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid
    with an infection threshold of 3 neighbors.
    """
    n = 14
    
    # For an n x n grid with a 3-neighbor infection threshold, the minimum number of
    # initial sites is 2n - 2 for even n, and 2n - 1 for odd n.
    # Since n=14 is even, we use the formula 2n - 2.
    
    if n % 2 == 0:
        result = 2 * n - 2
        print(f"For n = {n} (which is even), the minimum number of infected sites is calculated using the formula 2n - 2.")
        print(f"2 * {n} - 2 = {result}")
    else:
        # This part of the conditional is not executed for n=14, but included for completeness.
        result = 2 * n - 1
        print(f"For n = {n} (which is odd), the minimum number of infected sites is calculated using the formula 2n - 1.")
        print(f"2 * {n} - 1 = {result}")

solve_disease_spread()
