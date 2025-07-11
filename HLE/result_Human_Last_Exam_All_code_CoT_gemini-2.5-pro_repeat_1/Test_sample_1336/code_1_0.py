def solve_covering_groups():
    """
    This function calculates the total number of smooth coverings for the given problem.
    The problem is interpreted as finding the number of non-isomorphic quasi-simple
    covering groups of the simple group S = PSL(2, p) for a prime p > 5.
    """
    
    # The number of such coverings is given by the number of subgroups of the Schur
    # multiplier M(S) of the simple group S = PSL(2, p).
    
    # For a prime p > 5, it is a standard result in group theory that the Schur
    # multiplier of PSL(2, p) is the cyclic group of order 2, Z_2.
    schur_multiplier_order = 2
    
    print(f"The number of coverings is equal to the number of subgroups of the Schur multiplier of PSL(2, p).")
    print(f"For a prime p > 5, the order of the Schur multiplier is {schur_multiplier_order}.")
    
    # The number of subgroups of a cyclic group of order n is given by tau(n),
    # the number of divisors of n. We need to calculate tau(2).
    n = schur_multiplier_order
    
    # Calculate the number of divisors of n
    num_divisors = 0
    # In Python, we can find divisors by iterating up to n
    for i in range(1, n + 1):
        if n % i == 0:
            num_divisors += 1
            
    print(f"The number of subgroups of a cyclic group of order {n} is the number of its divisors, tau({n}).")
    
    # As requested, we output each number in the final equation.
    # The final equation is tau(2) = 2.
    print("\nThe final equation is:")
    print(f"tau({n}) = {num_divisors}")

    total_number = num_divisors
    print(f"\nThus, the total number of such smooth coverings is {total_number}.")

solve_covering_groups()