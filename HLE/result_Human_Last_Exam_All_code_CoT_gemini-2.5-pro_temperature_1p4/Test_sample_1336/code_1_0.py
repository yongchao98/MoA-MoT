import math

def solve():
    """
    Calculates the total number of smooth coverings for the given group theory problem.

    The number of smooth coverings corresponds to the number of subgroups of the
    center of the covering group G = SL(2, p).
    """

    # For p > 5 prime, the center of G = SL(2, p) is Z(G) = {I, -I},
    # which is isomorphic to the cyclic group Z_n where n=2.
    n = 2

    # The number of subgroups of a cyclic group Z_n is equal to the number of
    # positive divisors of n, which is given by the divisor function tau(n).
    # For a cyclic group, there is exactly one subgroup for each divisor of n.

    # Find the divisors of n.
    divisors = []
    for i in range(1, n + 1):
        if n % i == 0:
            divisors.append(i)

    # For each divisor, there is exactly one subgroup.
    # We create a list of 1s, one for each divisor.
    subgroup_counts_per_divisor = [1] * len(divisors)

    # The total number of subgroups is the sum of these counts.
    total_number_of_coverings = sum(subgroup_counts_per_divisor)

    # To satisfy the instruction "output each number in the final equation",
    # we construct and print the sum explicitly.
    equation_str = " + ".join(map(str, subgroup_counts_per_divisor))

    print(f"The group G = SL(2, p) is a covering group of S = PSL(2, p).")
    print(f"The center Z(G) is isomorphic to the cyclic group Z_n, where n = {n}.")
    print(f"The number of smooth coverings is the number of subgroups of Z_{n}.")
    print(f"The number of subgroups is calculated by summing the count for each divisor of {n}.")
    print(f"The final equation is: {equation_str} = {total_number_of_coverings}")
    print(f"Thus, the total number of such smooth coverings is {total_number_of_coverings}.")

solve()