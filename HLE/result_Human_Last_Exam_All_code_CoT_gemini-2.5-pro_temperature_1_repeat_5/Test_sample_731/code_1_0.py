import math

def solve_power_subgroups():
    """
    Calculates the number of power subgroups in the semidihedral group of size 512.
    """
    group_size = 512

    # For a semidihedral group SD_{2^n}, the order is 2^n.
    # We find n such that 2^n = group_size.
    if group_size <= 1 or (group_size & (group_size - 1)) != 0:
        print("Error: Group size must be a power of 2 greater than 1.")
        return

    n = int(math.log2(group_size))

    # Case 1: The exponent 'm' is odd.
    # For any odd m, the power subgroup G^m is the group G itself.
    # This contributes 1 unique subgroup.
    num_from_odd_m = 1

    # Case 2: The exponent 'm' is even.
    # For any even m, the power subgroup G^m is a cyclic subgroup <r^d>,
    # where d is an even divisor of 2^(n-1).
    # The number of such distinct subgroups is n-1.
    num_from_even_m = n - 1

    # Total number of power subgroups is the sum of the counts from both cases.
    total_subgroups = num_from_odd_m + num_from_even_m

    print("The number of power subgroups is calculated by considering two cases for the exponent 'm'.")
    print(f"1. For odd 'm': This results in {num_from_odd_m} unique subgroup (the group itself).")
    print(f"2. For even 'm': This results in {num_from_even_m} unique cyclic subgroups.")
    print("\nThe total number of power subgroups is the sum of these two cases.")
    print(f"{num_from_odd_m} + {num_from_even_m} = {total_subgroups}")

solve_power_subgroups()