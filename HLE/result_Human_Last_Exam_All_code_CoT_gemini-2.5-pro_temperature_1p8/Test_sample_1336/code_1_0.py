def count_divisors(n):
    """
    Calculates the number of positive divisors of a given integer n.
    This is also the number of subgroups of the cyclic group Z/nZ.
    """
    if n < 1:
        return 0
    count = 0
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            # If divisors are equal, count only one
            if n / i == i:
                count += 1
            # Otherwise, count both
            else:
                count += 2
    return count

def solve():
    """
    Solves the problem by finding the number of smooth coverings for PSL(2,p).
    """
    # The Schur multiplier of PSL(2, p) for prime p > 5 is Z/2Z.
    # The order of this group is 2.
    schur_multiplier_order = 2

    # The number of non-isomorphic covering groups is the number of subgroups
    # of the Schur multiplier. For a cyclic group of order n, this is the
    # number of divisors of n.
    num_coverings = count_divisors(schur_multiplier_order)

    print(f"The order of the Schur multiplier M(PSL(2,p)) for p > 5 is {schur_multiplier_order}.")
    print(f"The number of subgroups of a cyclic group of order {schur_multiplier_order} is the number of divisors of {schur_multiplier_order}.")
    print(f"Number of divisors of {schur_multiplier_order} is {num_coverings}.")
    print(f"Therefore, the total number of smooth coverings is {num_coverings}.")

solve()
