def solve_disease_spread():
    """
    Calculates the minimum number of initially infected sites for an n x n grid.
    This function explains the reasoning step-by-step and prints the final calculation.
    """
    n = 14
    
    print(f"Solving for an {n}x{n} grid.")
    print("Step 1: Analyze corner sites.")
    print("A corner site has only 2 neighbors. To be infected, a site needs at least 3.")
    print("Therefore, all 4 corner sites must be infected initially.")
    corner_sites = 4
    print(f"Corner sites to infect = {corner_sites}")
    
    print("\nStep 2: Analyze edge sites.")
    print("An edge site (not a corner) has 3 neighbors. For it to get infected, all 3 neighbors must be infected.")
    print("If two adjacent edge sites are healthy, neither can infect the other, creating a deadlock.")
    print("So, on each edge, we cannot have two consecutive healthy sites.")

    print("\nStep 3: Calculate the minimum number of infected sites on the edges.")
    sites_between_corners = n - 2
    print(f"Each of the 4 edges has {sites_between_corners} sites between the corners.")
    
    # To avoid two consecutive healthy sites on a path of length L,
    # we need to infect at least L/2 sites.
    min_infected_on_segment = sites_between_corners // 2
    print(f"Minimum sites to infect on each segment of {sites_between_corners} = {min_infected_on_segment}")

    total_edge_segment_sites = 4 * min_infected_on_segment
    print(f"Total sites to infect on all 4 segments = 4 * {min_infected_on_segment} = {total_edge_segment_sites}")

    print("\nStep 4: Calculate the total minimum number of sites.")
    total_minimum_sites = corner_sites + total_edge_segment_sites
    
    print("The total minimum number is the sum of infected corners and infected sites on the edge segments.")
    print(f"Total = (infected corners) + (infected sites on 4 edge segments)")
    print(f"Total = {corner_sites} + {total_edge_segment_sites}")
    print(f"The minimum number of initially infected sites is {total_minimum_sites}.")

solve_disease_spread()