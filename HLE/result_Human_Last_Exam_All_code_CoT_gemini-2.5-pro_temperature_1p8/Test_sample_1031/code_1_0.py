def solve_disease_spread():
    """
    Calculates the minimum number of initially infected sites to infect a 14x14 grid.
    
    The logic is based on a checkerboard pattern.
    1.  The grid size is n x n, where n = 14.
    2.  For a checkerboard pattern, there are n^2 / 2 squares of each color.
    3.  We infect all squares of one color (the one with fewer or equal corners).
        Let's say we infect the 'white' squares.
    4.  Then, we must also infect any corners of the other color ('black').
    5.  For n=14, there are 2 corners of each color.
    6.  The calculation is (number of white squares) + (number of black corners).
    """
    n = 14
    
    # Total sites
    total_sites = n * n
    
    # Number of sites of one color in a checkerboard pattern
    sites_per_color = total_sites // 2
    
    # For n=14, the corners (1,1) and (14,14) are one color,
    # and (1,14) and (14,1) are the other.
    # So there are 2 corners to add to the initial set of one color.
    additional_corners_needed = 2
    
    # Minimum number of initially infected sites
    min_infected_sites = sites_per_color + additional_corners_needed
    
    print(f"For an n x n grid with n = {n}:")
    print(f"Number of squares of one color (e.g., white) = n^2 / 2 = {sites_per_color}")
    print(f"Number of corners of the other color (black) to add = {additional_corners_needed}")
    print("The final equation for the minimum number of infected sites is:")
    print(f"{sites_per_color} + {additional_corners_needed} = {min_infected_sites}")

solve_disease_spread()