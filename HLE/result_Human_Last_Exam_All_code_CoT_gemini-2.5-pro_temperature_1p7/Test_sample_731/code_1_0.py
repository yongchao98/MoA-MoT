import math

def solve_power_subgroups():
    """
    Calculates the number of power subgroups in the semidihedral group of size 512.
    """
    group_order = 512

    # The semidihedral group is of the form SD_{2^n}.
    # We find n such that 2^n = group_order.
    n = int(math.log2(group_order))

    # For the semidihedral group SD_{2^n}, the number of distinct power subgroups is n.
    # This comes from the set of subgroups {G, G^2, G^4, ..., G^{2^(n-1)}}.
    # The subgroups are:
    # 1. The group G = SD_{2^n} itself.
    # 2. A set of n-1 distinct cyclic subgroups.
    num_non_cyclic = 1
    num_cyclic = n - 1
    total_subgroups = num_non_cyclic + num_cyclic

    print("The semidihedral group is of order 512 = 2^n.")
    print(f"Here, n = log2(512) = {n}.")
    print("\nThe number of power subgroups in SD_{2^n} is n.")
    print("This is composed of the group itself (non-cyclic) and n-1 distinct cyclic subgroups.")
    print("\nCalculation:")
    print(f"Number of non-cyclic power subgroups = {num_non_cyclic}")
    print(f"Number of cyclic power subgroups = n - 1 = {n} - 1 = {num_cyclic}")
    print(f"Total number of power subgroups = {num_non_cyclic} + {num_cyclic} = {total_subgroups}")
    
solve_power_subgroups()
