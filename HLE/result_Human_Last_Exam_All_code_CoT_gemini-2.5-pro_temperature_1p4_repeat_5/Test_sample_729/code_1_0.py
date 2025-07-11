def solve_power_subgroups():
    """
    Calculates the number of power subgroups in the generalized quaternion group Q_128.
    """
    group_order = 128
    # The group is Q_4n, so 4n = 128, which means n = 32.
    # The maximal cyclic subgroup <x> has order 2n = 64.
    cyclic_order = 64

    # Case 1: k is odd
    # For any odd exponent k, the power subgroup (Q_128)^k is Q_128 itself.
    # This accounts for 1 unique power subgroup.
    num_from_odd_k = 1
    print(f"For any odd exponent k, the power subgroup is the group Q_{group_order} itself.")
    print(f"Number of subgroups from odd k: {num_from_odd_k}\n")

    # Case 2: k is even
    # For any even exponent k, the power subgroup is <x^d>, where d = gcd(k, 64).
    # The possible values of d are the even divisors of 64.
    # We will now count the number of even divisors of 64.
    print(f"For any even exponent k, the power subgroups are cyclic subgroups of the form <x^d>,")
    print(f"where d is an even divisor of {cyclic_order}.")

    even_divisors = []
    for i in range(1, cyclic_order + 1):
        if cyclic_order % i == 0:  # Check if i is a divisor
            if i % 2 == 0:        # Check if the divisor is even
                even_divisors.append(i)

    num_from_even_k = len(even_divisors)
    print(f"The even divisors of {cyclic_order} are: {even_divisors}")
    print(f"Number of subgroups from even k: {num_from_even_k}\n")

    # Total number of power subgroups is the sum from both cases.
    total_subgroups = num_from_odd_k + num_from_even_k

    # Output the final equation as requested.
    print(f"The total number of power subgroups is the sum from both cases.")
    print(f"Total = {num_from_odd_k} + {num_from_even_k} = {total_subgroups}")

solve_power_subgroups()