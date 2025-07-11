def solve_infection_problem(n):
    """
    Calculates the minimum number of initially infected sites to infect an n x n grid.
    The rule is a site gets infected if it has at least 3 infected neighbors.
    This function implements the logic for an even n.
    """
    if n % 2 != 0:
        # The logic here is derived for even n, as per the problem (n=14).
        # A different derivation would be needed for odd n.
        print("This solution is derived for even n.")
        return

    # Step 1: Corners must be infected.
    num_corners = 4

    # Step 2: On each edge, every other site must be infected.
    # Number of non-corner sites on one edge is n - 2.
    # We infect half of them.
    num_edge_sites_per_edge = (n - 2) // 2
    num_boundary_sites = num_corners + 4 * num_edge_sites_per_edge
    print(f"The number of initially infected sites on the boundary is {num_boundary_sites}.")
    print(f"This consists of {num_corners} corner sites and {4 * num_edge_sites_per_edge} other sites on the edges.")


    # Step 3: Calculate the number of support sites.
    # These form a frame on rows/cols 1 and n-2.
    # Number of odd indices between 1 and n-2.
    # For n=14, indices are 1..12. Odd numbers are 1,3,5,7,9,11 (6 sites)
    num_odd_indices = (n - 2) // 2
    
    # We have 4 sets of support sites
    # |S1 U S2 U S3 U S4| = |S1|+|S2|+|S3|+|S4| - |overlaps|
    # For even n >= 4, the only overlap is at (1,1).
    num_support_sites = 4 * num_odd_indices - 1
    
    print(f"\nThe number of 'support' sites needed in the next layer is {num_support_sites}.")

    # Step 4: Calculate the total number of sites.
    total_sites = num_boundary_sites + num_support_sites
    
    print(f"\nTotal minimum number of sites = (Boundary Sites) + (Support Sites)")
    print(f"Total minimum number of sites = {num_boundary_sites} + {num_support_sites} = {total_sites}")


# For the given problem, n = 14.
n = 14
solve_infection_problem(n)