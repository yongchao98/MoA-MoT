def solve_grid_infection():
    """
    Calculates the minimum number of sites to infect an n x n grid
    with a 3-neighbor infection rule.
    """
    # The size of the grid is n x n.
    n = 14

    print("The problem is to find the minimum number of initially infected sites to infect an entire n x n grid.")
    print("The rule for infection is: a healthy site becomes infected if it has at least 3 infected neighbors.")
    print("\nThis is a known problem in the mathematical field of bootstrap percolation.")
    print("For an n x n grid and a 3-neighbor rule, the minimum number of sites has been proven to be:")
    print(" - 'n' if n is even")
    print(" - 'n + 1' if n is odd")
    print("\nWe need to solve this for a grid of size n = 14.")

    # Determine if n is even or odd and apply the formula.
    if n % 2 == 0:
        print(f"\nSince n = {n} is an even number, the minimum number of sites is equal to n.")
        result = n
        # The final equation is simply the value of n.
        print(f"Final Answer = {n}")
    else:
        # This case is not applicable for n=14 but included for completeness.
        print(f"\nSince n = {n} is an odd number, the minimum number of sites is n + 1.")
        result = n + 1
        # The final equation shows the addition.
        print(f"Final Answer = {n} + 1 = {result}")

solve_grid_infection()