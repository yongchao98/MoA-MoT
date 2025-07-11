import math

def count_power_subgroups_of_semidihedral():
    """
    Calculates the number of power subgroups in a semidihedral group of a given size.

    The number of power subgroups for the semidihedral group SD_2^n is n.
    This script finds n from the group size and outputs the total count.
    """
    group_size = 512

    # A semidihedral group has size 2^n. We find n for the given size.
    # 2^n = 512
    n = int(math.log2(group_size))

    # The power subgroups of SD_2^n are G itself, and the n-1 cyclic subgroups
    # <r^2>, <r^4>, ..., <r^(2^(n-1))>. The last one is the trivial group {1}.
    # This gives a total of 1 + (n-1) = n subgroups.

    num_proper_subgroups = n - 1
    total_subgroups = 1 + num_proper_subgroups

    print(f"For the semidihedral group of size {group_size}, n is {n}.")
    print("The number of power subgroups is equal to n.")
    print("This is composed of the group G itself, plus n-1 other distinct subgroups.")
    print("\nCalculation:")
    # The final required equation output
    print(f"1 + ({n} - 1) = {total_subgroups}")

count_power_subgroups_of_semidihedral()
