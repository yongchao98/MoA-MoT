def solve_disease_spread():
    """
    Calculates the minimum number of initially infected sites to infect a 14x14 grid.
    
    The logic is as follows:
    1.  The grid is n x n, with n=14.
    2.  A healthy site becomes infected if it has at least 3 infected neighbors.
    3.  Corner sites have only 2 neighbors, so they can never become infected by this rule.
        Therefore, all 4 corners must be in the initial set.
    4.  A highly effective strategy is to use a checkerboard pattern. We infect all sites (i, j)
        where i + j is even. On a 14x14 grid, this accounts for half of the sites.
    5.  The corners are (0,0), (0,13), (13,0), (13,13).
        - (0,0): 0+0=0 (even) - included in the checkerboard pattern.
        - (13,13): 13+13=26 (even) - included in the checkerboard pattern.
        - (0,13): 0+13=13 (odd) - not included.
        - (13,0): 13+0=13 (odd) - not included.
    6.  To ensure the "odd" corners are infected, we must add them to the initial set.
    7.  Any "odd" site that is not a corner has at least 3 neighbors, all of which are "even"
        and thus in our initial set. So, all other healthy sites will be infected in the first step.
    8.  The total number is the number of "even" sites plus the number of "odd" corners.
    """
    n = 14
    
    # Total sites on the grid
    total_sites = n * n
    
    # Number of sites where i+j is even (checkerboard "black" squares)
    # For an even n, this is exactly half the total sites.
    even_sites_count = total_sites // 2
    
    # The two corners (0, 13) and (13, 0) have an odd sum of coordinates.
    # They are not included in the 'even_sites' set and must be added.
    extra_corners_needed = 2
    
    # The minimum number of sites is the sum of the "even" sites and the extra corners.
    min_infected_sites = even_sites_count + extra_corners_needed
    
    # Print the explanation and the final equation.
    print(f"For an n x n grid with n = {n}:")
    print(f"The number of sites in a checkerboard pattern (where i+j is even) is n*n / 2.")
    print(f"Number of 'even' sites = {n} * {n} / 2 = {even_sites_count}")
    print(f"Number of corners not covered by this pattern = {extra_corners_needed}")
    print(f"The minimum number of sites is the sum of these two numbers.")
    print(f"Final Equation: {even_sites_count} + {extra_corners_needed} = {min_infected_sites}")
    
solve_disease_spread()
