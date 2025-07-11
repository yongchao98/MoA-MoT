import math

def count_power_subgroups(group_size):
    """
    Calculates the number of power subgroups in the semidihedral group of a given size.

    Args:
        group_size (int): The order of the semidihedral group, must be a power of 2.
    """
    if group_size <= 0 or (group_size & (group_size - 1)) != 0:
        print("Error: The size of the semidihedral group must be a power of 2.")
        return

    # For a semidihedral group SD_{2^n}, the size is 2^n. We find n.
    n = int(math.log2(group_size))

    print(f"The semidihedral group G of size {group_size} is denoted as SD_{{{group_size}}}, or SD_{{2^{n}}}. So, n = {n}.")
    print("A power subgroup is of the form G^k = {g^k | g in G} for an integer k >= 1.")
    print("-" * 50)

    # Case 1: k is odd
    # For any odd k, G^k = G, since the map g -> g^k is a permutation on the 2-group G.
    # This gives us one unique power subgroup: G itself.
    num_from_odd_k = 1
    print(f"1. When k is odd:")
    print(f"   G^k = G. This accounts for {num_from_odd_k} power subgroup.")

    # Case 2: k is even
    # Any even k can be written as k = 2^j * m, with m odd and j >= 1.
    # G^k = G^(2^j). We need to count the distinct subgroups G^(2^j).
    # For SD_{2^n}, G^(2^j) = <r^(2^j)> for j=1, ..., n-1.
    # where <r> is the maximal cyclic subgroup of order 2^(n-1).
    num_from_even_k = n - 1
    print(f"\n2. When k is even:")
    print(f"   The distinct power subgroups are G^2, G^4, ..., G^(2^{n-1}).")
    print(f"   This gives {num_from_even_k} distinct subgroups, from <r^2> down to the trivial group {{e}}.")

    # Final Calculation
    total = num_from_odd_k + num_from_even_k

    print("-" * 50)
    print("The total number of distinct power subgroups is the sum of these two cases.")
    print(f"\nTotal = {num_from_odd_k} (from odd k's) + {num_from_even_k} (from even k's) = {total}")

# Run the calculation for the semidihedral group of size 512
count_power_subgroups(512)