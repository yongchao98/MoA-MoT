import math

def count_power_subgroups_semidihedral():
    """
    Calculates the number of power subgroups in the semidihedral group of size 512.
    """
    group_size = 512

    # A semidihedral group SD_{2^n} has a size of 2^n. We first find n.
    # n = log2(size)
    n = int(math.log2(group_size))

    print(f"The semidihedral group is of size {group_size}, which can be written as 2^n.")
    print(f"First, we find the value of n for the group SD_{{{group_size}}}:")
    print(f"n = log2({group_size}) = {n}")

    # For a semidihedral group SD_{2^n}, the number of power subgroups is n.
    # This count is composed of:
    # 1. The group G itself.
    # 2. n-1 distinct cyclic subgroups.
    num_non_trivial_subgroups = 1
    num_cyclic_subgroups = n - 1
    total_subgroups = num_non_trivial_subgroups + num_cyclic_subgroups

    print("\nThe number of power subgroups for SD_{2^n} is exactly n.")
    print("This total is the sum of:")
    print(f"  {num_non_trivial_subgroups} (for the group itself)")
    print(f"  {num_cyclic_subgroups} (for the distinct cyclic power subgroups)")

    print("\nThe final equation is:")
    print(f"{num_non_trivial_subgroups} + {num_cyclic_subgroups} = {total_subgroups}")

count_power_subgroups_semidihedral()