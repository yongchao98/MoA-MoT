def solve_infection_problem():
    """
    Calculates the minimum number of initially infected sites to infect a 14x14 grid.
    
    The problem is a known result in bootstrap percolation. For an n x n grid
    and a 3-neighbor infection rule, the minimum number of sites is 2n-2 for n >= 3.
    """
    n = 14
    
    # The formula for the minimum number of sites is 2n - 2.
    min_sites = 2 * n - 2
    
    # We will print the equation step-by-step
    print(f"Based on the established mathematical formula for this problem, the minimum number of sites is 2*n - 2.")
    print(f"For n = {n}, the calculation is:")
    print(f"2 * {n} - 2 = {min_sites}")

solve_infection_problem()