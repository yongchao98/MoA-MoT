def solve_grid_infection():
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid.

    The problem is a known case of bootstrap percolation. For an n x n grid
    with the rule that a site becomes infected if it has at least 3 infected
    von Neumann neighbors (up, down, left, right), the minimum number of
    sites required to infect the entire grid is n.
    """
    
    # The size of the grid is n x n.
    n = 14
    
    # According to known results in bootstrap percolation, the minimum number
    # of initially infected sites required is n.
    min_sites = n
    
    # The final equation is simply that the result is equal to n.
    # We will print the numbers involved as requested.
    print(f"Grid size (n): {n}")
    print(f"Formula for minimum sites: n")
    print(f"Final Equation: {n} = {min_sites}")
    print(f"The minimum number of initially infected sites is: {min_sites}")

solve_grid_infection()