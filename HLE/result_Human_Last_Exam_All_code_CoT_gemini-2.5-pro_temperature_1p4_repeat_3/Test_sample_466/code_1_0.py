def solve():
    """
    Calculates the number of internal adjunctions in the given 2-category
    from F_11^3 to itself, which is equivalent to computing the order of
    the general linear group GL(3, 11).
    """
    q = 11
    n = 3

    # Calculate the terms in the formula for the order of GL(n, q)
    # |GL_n(F_q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^{n-1})
    q_n = q**n
    
    term1_val = q**0
    term1 = q_n - term1_val
    
    term2_val = q**1
    term2 = q_n - term2_val
    
    term3_val = q**2
    term3 = q_n - term3_val
    
    # Calculate the final result
    order = term1 * term2 * term3

    # Print the step-by-step calculation
    print(f"The number of internal adjunctions is the order of the group GL(n, F_q) for n={n}, q={q}.")
    print(f"|GL({n}, F_{q})| = (q^n - q^0) * (q^n - q^1) * (q^n - q^2)")
    print(f"             = ({q}^{n} - {q}^0) * ({q}^{n} - {q}^1) * ({q}^{n} - {q}^2)")
    print(f"             = ({q_n} - {term1_val}) * ({q_n} - {term2_val}) * ({q_n} - {term3_val})")
    print(f"             = {term1} * {term2} * {term3}")
    print(f"The total number of internal adjunctions is: {order}")

solve()