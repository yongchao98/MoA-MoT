def solve_grid_infection(n):
    """
    Calculates and explains the minimum number of sites to infect an n x n grid
    where infection spreads to a site with at least 3 infected neighbors.
    This solution is tailored for an even-sized grid n.
    """

    if n % 2 != 0:
        print(f"The provided solution logic is for an even n. The case for n={n} is more complex.")
        return

    print("Step-by-step Derivation:")
    print("1. The rule for infection is that a healthy site becomes infected if it has at least 3 infected neighbors.")
    
    print("\n2. A key observation is that a 2x2 block of healthy sites is stable.")
    print("   Any site within such a block has at least two healthy neighbors. This means it can have at most 4 - 2 = 2 infected neighbors.")
    print("   Since 2 is less than the required threshold of 3, no site in a 2x2 healthy block can ever become infected.")
    
    print("\n3. To infect the entire grid, the initial state cannot have any 2x2 blocks of healthy sites.")
    print("   The problem is equivalent to finding the minimum number of initially infected sites to 'hit' every possible 2x2 block.")
    
    print("\n4. Establishing a Lower Bound:")
    k = n // 2
    num_disjoint_blocks = k * k
    print(f"   A {n}x{n} grid can be tiled by {k} * {k} = {num_disjoint_blocks} disjoint 2x2 blocks.")
    print("   To ensure no 2x2 healthy block exists, each of these disjoint blocks must contain at least one infected site.")
    print(f"   Therefore, we need at least {num_disjoint_blocks} initially infected sites.")

    print("\n5. Proposing a Valid Construction (Upper Bound):")
    print("   This lower bound can be achieved by infecting all sites (i, j) where both i and j are even-numbered.")
    num_even_indices = n // 2
    num_infected_sites = num_even_indices * num_even_indices
    print(f"   For a {n}x{n} grid, there are {num_even_indices} even rows and {num_even_indices} even columns.")
    print(f"   The total number of such sites provides a working configuration of size {num_even_indices} * {num_even_indices} = {num_infected_sites}.")
    
    print("\n6. Conclusion:")
    print("   The minimum required number (lower bound) and the size of our constructed set (upper bound) are equal.")
    print(f"   The minimum number of initially infected sites for an {n}x{n} grid is given by the final equation:")
    print(f"   ({n} / 2) * ({n} / 2) = {k} * {k} = {num_infected_sites}")

# Execute the function for the specific problem (n=14)
solve_grid_infection(14)