import math

def solve():
    """
    Calculates the minimum number of initially infected sites for an n x n grid.
    The problem is equivalent to finding the minimum size of a set of sites
    that covers all 2x2 subgrids. A known optimal pattern for this is to
    infect all sites (i, j) where (i + 2*j) is divisible by 4.
    This function calculates the number of such sites for n=14.
    """
    n = 14

    # We count the number of sites (i,j) where i,j in {0..n-1}
    # that satisfy (i + 2*j) % 4 == 0.

    # This is equivalent to i % 4 == (-2*j) % 4.
    # We analyze based on the parity of j.

    # Count the number of even and odd indices for j in the range [0, n-1]
    num_j_even = math.ceil(n / 2)
    num_j_odd = math.floor(n / 2)

    # For each even j (j=2k), the condition is i % 4 == (-4k) % 4 == 0.
    # We count how many i in [0, n-1] satisfy i % 4 == 0.
    # i can be 0, 4, 8, ...
    count_i_for_even_j = math.floor((n - 1) / 4) + 1

    # For each odd j (j=2k+1), the condition is i % 4 == (-4k-2) % 4 == 2.
    # We count how many i in [0, n-1] satisfy i % 4 == 2.
    # i can be 2, 6, 10, ...
    count_i_for_odd_j = math.floor((n - 1 - 2) / 4) + 1 if n > 2 else 0

    # Calculate the total number of infected sites
    sites_for_even_j = num_j_even * count_i_for_even_j
    sites_for_odd_j = num_j_odd * count_i_for_odd_j
    total_infected = sites_for_even_j + sites_for_odd_j

    # Print the step-by-step calculation
    print(f"For an n x n grid with n = {n}:")
    print("We use the pattern (i + 2*j) % 4 == 0 to select initial sites.")
    print("\nCalculation steps:")
    print(f"Number of even indices for j in [0, {n-1}]: {num_j_even}")
    print(f"Number of odd indices for j in [0, {n-1}]: {num_j_odd}")
    print(f"For each even j, we need i % 4 == 0. Number of such i's in [0, {n-1}]: {count_i_for_even_j}")
    print(f"For each odd j, we need i % 4 == 2. Number of such i's in [0, {n-1}]: {count_i_for_odd_j}")
    print("\nFinal equation:")
    print(f"Total initially infected sites = ({num_j_even} * {count_i_for_even_j}) + ({num_j_odd} * {count_i_for_odd_j}) = {sites_for_even_j} + {sites_for_odd_j} = {total_infected}")

solve()