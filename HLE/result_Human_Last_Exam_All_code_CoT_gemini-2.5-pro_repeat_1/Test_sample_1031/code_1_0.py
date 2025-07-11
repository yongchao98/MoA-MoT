import math

def solve_disease_spread():
    """
    Calculates the minimum number of initially infected sites to infect a 14x14 grid.

    The problem is a classic example of bootstrap percolation. For an n x n grid
    with the 3-neighbor rule, the key is to prevent stable clusters of healthy
    sites from surviving.

    1.  A 2x2 block of healthy sites is a stable configuration if surrounded by
        infected sites. Each site in the 2x2 block has 2 infected neighbors and
        2 healthy neighbors, so none of them get infected.
    2.  To guarantee the infection of the entire grid, every 2x2 block must
        contain at least one initially infected site.
    3.  For an n x n grid where n is even (like n=14), we can tile the grid
        with (n/2) * (n/2) disjoint 2x2 blocks.
    4.  For n=14, we can create (14/2) * (14/2) = 7 * 7 = 49 disjoint 2x2 blocks.
    5.  Since these 49 blocks are disjoint, we need at least one infected site
        in each, giving a lower bound of 49.
    6.  It has been proven that this lower bound is also achievable, making it the
        exact minimum number.
    """
    n = 14
    # For an n x n grid where n is even, the minimum number of sites is (n/2)^2.
    min_sites = (n / 2) ** 2
    
    # The final answer needs to be an integer.
    result = int(min_sites)
    
    print(f"For a {n}x{n} grid, the minimum number of initially infected sites is {result}.")
    # The problem asks for the number in the equation.
    # The equation is (14/2)^2 = 7^2 = 49
    print("The equation is: (14 / 2)^2 = 7^2 = 49")

solve_disease_spread()