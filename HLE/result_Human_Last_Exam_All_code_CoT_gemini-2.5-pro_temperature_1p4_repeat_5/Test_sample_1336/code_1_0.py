import math

def get_divisors(n):
    """
    Calculates all positive divisors of an integer n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_covering_groups():
    """
    Solves the problem by finding the number of covering groups for PSL(2,p).
    """
    # For a prime p > 5, the simple group is S = PSL(2, p).
    # G = SL(2, p) is a quasi-simple group that acts as a covering group of S.
    # The number of non-isomorphic covering groups of a perfect group S is
    # determined by the number of subgroups of its Schur multiplier, M(S).

    # For S = PSL(2, p) with p > 5, the Schur multiplier M(S) is the
    # cyclic group of order 2.
    schur_multiplier_order = 2

    # The number of subgroups of a cyclic group of order n is equal to
    # the number of divisors of n.
    divisors_of_order = get_divisors(schur_multiplier_order)
    num_coverings = len(divisors_of_order)

    print("The total number of such smooth coverings corresponds to the number of non-isomorphic covering groups of PSL(2, p).")
    print("This number is equal to the number of subgroups of the Schur multiplier of PSL(2, p).")
    print(f"For p > 5, the Schur multiplier's order is {schur_multiplier_order}.")
    print(f"The number of subgroups is the number of divisors of {schur_multiplier_order}.")
    
    # The final "equation" is finding the size of the set of divisors.
    # The numbers in the equation are the divisors themselves.
    divisors_str = ", ".join(map(str, divisors_of_order))
    print(f"The divisors of {schur_multiplier_order} are: {divisors_str}.")
    print(f"Counting these divisors ({divisors_str}) gives the total number.")
    print(f"The total number of such smooth coverings is {num_coverings}.")

solve_covering_groups()
<<<2>>>