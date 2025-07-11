import collections

def solve_disease_spread():
    """
    Calculates the number of initially infected sites for an n x n grid.

    The logic is based on ensuring no stable 2x2 block of healthy cells exists
    and that all four corners are initially infected. A checkerboard pattern with
    a small correction for the corners is used.
    """
    n = 14
    infected_sites = set()

    # A checkerboard pattern where i and j have the same parity ensures
    # no 2x2 block can be entirely healthy.
    # Number of (odd, odd) sites: ceil(n/2) * ceil(n/2) = 7 * 7 = 49
    # Number of (even, even) sites: floor(n/2) * floor(n/2) = 7 * 7 = 49
    # Total for the base checkerboard pattern:
    num_same_parity = (n // 2)**2 + ((n + 1) // 2)**2
    
    # Add sites where coordinates have the same parity.
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if (i % 2) == (j % 2):
                infected_sites.add((i, j))

    # Corner sites must be infected.
    # The same-parity checkerboard pattern infects (1,1) and (14,14).
    # We must explicitly add the other two corners, (1,14) and (14,1),
    # which have coordinates of different parity.
    corners_to_add = set()
    if (1, n) not in infected_sites:
        corners_to_add.add((1,n))
    if (n, 1) not in infected_sites:
        corners_to_add.add((n,1))
    
    infected_sites.update(corners_to_add)

    # Calculate the final number of infected sites.
    # We have 98 sites from the checkerboard pattern.
    # We add 2 corners that were not in the pattern.
    # The final equation is 98 + 2 = 100.
    final_count = len(infected_sites)

    num_added_corners = len(corners_to_add)

    # Print out the calculation step-by-step
    print(f"For a {n}x{n} grid:")
    print(f"1. A checkerboard pattern with same-parity coordinates is chosen. This consists of sites (i,j) where i and j are both even or both odd.")
    print(f"   Number of sites in this base pattern: {num_same_parity}")
    print(f"2. The 4 corner sites must be infected. The base pattern includes (1,1) and (14,14) but not (1,14) and (14,1). We must add the missing corners.")
    print(f"   Number of corners to add: {num_added_corners}")
    print("3. The total minimum number of sites is the sum of sites from the base pattern and the additionally required corners.")
    print(f"Final calculation: {num_same_parity} + {num_added_corners} = {final_count}")
    print("\nFinal Answer:")
    print(final_count)

solve_disease_spread()
