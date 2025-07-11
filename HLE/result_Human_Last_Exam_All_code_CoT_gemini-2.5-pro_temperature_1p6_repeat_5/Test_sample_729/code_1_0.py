import math

def solve_quaternion_power_subgroups():
    """
    Calculates the number of power subgroups in the generalized quaternion group
    of a given size.
    """
    group_size = 128

    # The size of a generalized quaternion group is always a power of 2, denoted as 2^n.
    # We find n by taking the logarithm base 2 of the group size.
    # The float result is converted to an integer.
    n = int(math.log2(group_size))

    # For a generalized quaternion group Q_{2^n}, the number of distinct power subgroups
    # is exactly n. These subgroups correspond to G^(2^i) for i = 0, 1, ..., n-1.
    num_power_subgroups = n

    print(f"The generalized quaternion group of size {group_size} is Q_{2**n}.")
    print(f"The number of power subgroups in Q_{2**n} is equal to n.")
    print("\nThe final equation is:")
    
    # Output each number in the final equation as requested.
    print(f"Number of power subgroups = log2({group_size}) = {num_power_subgroups}")

solve_quaternion_power_subgroups()