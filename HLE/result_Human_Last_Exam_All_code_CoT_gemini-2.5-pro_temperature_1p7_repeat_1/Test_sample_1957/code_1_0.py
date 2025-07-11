def solve_min_cost():
    """
    Calculates the minimum value of 1000m + n based on the problem description.
    """
    # The degree 'd' of the irreducible polynomial P defining the set S is 10.
    d = 10

    # The constraint derived from the degrees of the polynomials is n * m = d * k,
    # where k is a positive integer. So, n * m = 10 * k.

    # We want to minimize the cost function C = 1000*m + n.
    # The term '1000*m' dominates the cost. To minimize C, we should choose
    # the smallest possible positive integer for m.
    m_final = 1

    # With m=1, the constraint becomes n = 10 * k.
    # The cost becomes C = 1000 * 1 + 10 * k.
    # To minimize this, we must choose the smallest positive integer for k.
    k_final = 1

    # This determines the value of n.
    n_final = d * k_final

    # The existence of a polynomial map F with n=10, m=1 is guaranteed by
    # the fact that any polynomial of degree d can be written as the
    # determinant of a d x d matrix of linear forms.

    # Now, we calculate the minimum value.
    min_value = 1000 * m_final + n_final
    
    print(f"To minimize the expression 1000*m + n, we determined the optimal integer values for m and n.")
    print(f"The values are m = {m_final} and n = {n_final}.")
    print(f"The final calculation is: 1000 * {m_final} + {n_final} = {min_value}")

solve_min_cost()