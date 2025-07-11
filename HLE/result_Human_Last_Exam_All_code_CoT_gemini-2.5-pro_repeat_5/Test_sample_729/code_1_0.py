import math

def count_power_subgroups_q128():
    """
    This function calculates and explains the number of power subgroups
    in the generalized quaternion group of size 128 (Q_128).
    """

    group_size = 128

    # For Q_{2^n}, the group size is 2^n. For Q_128, n=7.
    n = int(math.log2(group_size))

    # The order of the cyclic subgroup <x> is 2^(n-1).
    x_order = 2**(n-1)

    print(f"The group is the generalized quaternion group Q_{group_size} (or Q_{{2**n}} with n={n}).")
    print(f"It has a presentation <x, y | x^{x_order}=1, x^{x_order//2}=y^2, ...>.")
    print("We want to find the number of its power subgroups, which are subgroups of the form G^m = {g^m | g in G}.\n")

    # Case 1: m is odd
    print("Step 1: Consider the case where the power 'm' is an odd integer.")
    print("For any odd 'm', the set G^m is the entire group G itself.")
    print("This gives us one power subgroup: Q_128.")
    odd_m_subgroups = 1
    print(f"Number of unique subgroups from odd m = {odd_m_subgroups}\n")

    # Case 2: m is even
    print("Step 2: Consider the case where the power 'm' is an even integer.")
    print(f"For any even 'm', the set G^m is a cyclic subgroup of <x> of the form <x^d>, where d = gcd(m, {x_order}).")
    print(f"Since 'm' is even, d must be an even divisor of {x_order}.")
    
    # Find all even divisors of x_order (64)
    divisors = [i for i in range(1, x_order + 1) if x_order % i == 0]
    even_divisors = [d for d in divisors if d % 2 == 0]
    even_m_subgroups = len(even_divisors)

    print(f"The even divisors of {x_order} are: {even_divisors}.")
    print(f"Each of these {even_m_subgroups} divisors corresponds to a distinct power subgroup.")
    print("These subgroups are <x^2>, <x^4>, <x^8>, <x^16>, <x^32>, and <x^64> (which is the trivial subgroup {1}).")
    print(f"Number of unique subgroups from even m = {even_m_subgroups}\n")

    # Final Calculation
    total_subgroups = odd_m_subgroups + even_m_subgroups

    print("Step 3: Calculate the total number of power subgroups.")
    print("The total is the sum of the counts from the odd and even cases.")
    print(f"Total subgroups = (subgroups from odd m) + (subgroups from even m)")
    print(f"Total subgroups = {odd_m_subgroups} + {even_m_subgroups} = {total_subgroups}")

if __name__ == '__main__':
    count_power_subgroups_q128()