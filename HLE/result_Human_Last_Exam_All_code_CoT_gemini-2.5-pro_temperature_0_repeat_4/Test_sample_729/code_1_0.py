import math

def count_power_subgroups_in_q128():
    """
    Calculates and explains the number of power subgroups in the
    generalized quaternion group of size 128.
    """
    group_size = 128

    # For a generalized quaternion group Q_{2^n}, the size is 2^n.
    # We find n for the given group size.
    # 2^n = 128
    n = int(math.log2(group_size))

    print(f"The group is the generalized quaternion group of size {group_size}, denoted Q_{{{group_size}}}.")
    print(f"This corresponds to the general form Q_{{2^n}} with n = {n}.")
    print("-" * 40)
    print("To find the number of power subgroups, we analyze G^k = {g^k | g in G}:")
    print("")

    # Case 1: k is odd
    num_odd_k_subgroups = 1
    print(f"1. For any odd exponent k, the set of k-th powers G^k is the group G itself.")
    print(f"   This gives {num_odd_k_subgroups} power subgroup: Q_{{{group_size}}}.")
    print("")

    # Case 2: k is even
    # The number of such subgroups is the number of even divisors of 2^(n-1).
    # The even divisors are 2^1, 2^2, ..., 2^(n-1). There are n-1 of them.
    num_even_k_subgroups = n - 1
    order_of_x = 2**(n - 1)
    print(f"2. For any even exponent k, G^k is a cyclic subgroup of the form <x^d>,")
    print(f"   where d is an even divisor of 2^(n-1) = {order_of_x}.")
    print(f"   The even divisors of {order_of_x} are 2, 4, 8, 16, 32, and 64.")
    print(f"   This gives {num_even_k_subgroups} distinct power subgroups.")
    print("-" * 40)

    # Total number of power subgroups
    total_subgroups = num_odd_k_subgroups + num_even_k_subgroups

    print("The total number of distinct power subgroups is the sum from both cases:")
    print(f"   {num_odd_k_subgroups} (from odd k) + {num_even_k_subgroups} (from even k) = {total_subgroups}")
    print("")
    print(f"Final Answer: There are {total_subgroups} power subgroups in the generalized quaternion group of size 128.")

count_power_subgroups_in_q128()