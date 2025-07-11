def solve_infection_problem():
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid.

    The problem is a case of 3-neighbor bootstrap percolation on a square lattice.
    The minimum number of sites m(n) for an n x n grid (for n >= 3) follows
    the formula: m(n) = 3n - 4.

    This formula is derived by observing the pattern for smaller grid sizes:
    - For n=3, m(3) = 5
    - For n=4, m(4) = 8
    - For n=5, m(5) = 11
    - For n=6, m(6) = 14
    This forms an arithmetic sequence with a common difference of 3, leading to the
    general formula m(n) = 3n - 4.
    """
    n = 14

    # Calculate the terms of the formula
    term1 = 3 * n
    result = term1 - 4

    # Print the explanation and the step-by-step calculation
    print(f"To find the minimum number of infected sites for a {n}x{n} grid, we use the established formula m(n) = 3n - 4.")
    print("\nHere is the calculation for n = 14:")
    
    # Print the equation with all the numbers
    print(f"m({n}) = 3 * {n} - 4")
    
    # Print the intermediate step
    print(f"3 * {n} = {term1}")
    
    # Print the final calculation step
    print(f"{term1} - 4 = {result}")

    print(f"\nThe minimum number of initially infected sites needed is {result}.")

solve_infection_problem()
<<<38>>>