def solve():
    """
    Calculates the total number of quasi-simple covering groups for PSL(2, p) with p > 5 prime.

    This number corresponds to the number of subgroups of the Schur Multiplier of PSL(2, p),
    which is Z_2 for p > 5.
    """

    # The Schur Multiplier M is the cyclic group of order 2, Z_2.
    schur_multiplier_order = 2

    # The number of non-isomorphic quasi-simple covering groups is equal to the
    # number of subgroups of the Schur Multiplier. For a cyclic group Z_n, this
    # is the number of divisors of n.
    
    # We find the number of divisors of 2.
    num_subgroups = 0
    for i in range(1, schur_multiplier_order + 1):
        if schur_multiplier_order % i == 0:
            num_subgroups += 1

    # The two subgroups are the trivial subgroup {0} and the whole group Z_2.
    # Each subgroup corresponds to one quasi-simple covering group.
    
    # Contribution from the trivial subgroup (leading to SL(2,p))
    count1 = 1
    # Contribution from the full subgroup Z_2 (leading to PSL(2,p))
    count2 = 1
    
    # The total number is the sum
    total_coverings = count1 + count2

    # Print the explanation and the final equation as requested.
    print(f"The Schur Multiplier of PSL(2, p) for p > 5 is Z_2.")
    print(f"The number of quasi-simple covering groups is the number of subgroups of Z_2, which is {num_subgroups}.")
    print(f"The trivial subgroup gives rise to 1 covering group (SL(2,p)).")
    print(f"The full subgroup gives rise to 1 covering group (PSL(2,p)).")
    print("The total number of such smooth coverings is the sum:")
    print(f"{count1} + {count2} = {total_coverings}")

solve()