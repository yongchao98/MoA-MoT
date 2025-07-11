def solve_inertial_quotient_order():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    # Parameters from the problem
    # n is the dimension of the vector space, derived from D being (C_2)^n
    # The order of D is 16 = 2^4, so n = 4.
    n = 4
    # q is the characteristic of the field k, which is 2.
    q = 2

    # The highest possible order for E is the order of GL(n, F_q).
    # We calculate this using the formula:
    # |GL(n, q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))

    q_n = q**n
    terms = []
    for i in range(n):
        terms.append(q_n - q**i)

    # Calculate the final result
    order = 1
    for term in terms:
        order *= term

    # Print the explanation and the step-by-step calculation
    print("The highest possible order for the inertial quotient E is the order of the general linear group GL(4, F_2).")
    print("The order is calculated as follows:")
    
    # Build the formula string
    formula_str = f"|GL({n}, F_{q})| = "
    for i in range(n):
        formula_str += f"({q}^{n} - {q}^{i})"
        if i < n - 1:
            formula_str += " * "
    print(formula_str)

    # Build the calculation string with intermediate values
    calc_str = "= "
    for i in range(n):
        calc_str += f"({q_n} - {q**i})"
        if i < n - 1:
            calc_str += " * "
    print(calc_str)

    # Build the final equation string
    final_eq_str = "= "
    for i in range(len(terms)):
        final_eq_str += str(terms[i])
        if i < len(terms) - 1:
            final_eq_str += " * "
    print(final_eq_str)
    
    # Print the final answer
    print(f"= {order}")

solve_inertial_quotient_order()