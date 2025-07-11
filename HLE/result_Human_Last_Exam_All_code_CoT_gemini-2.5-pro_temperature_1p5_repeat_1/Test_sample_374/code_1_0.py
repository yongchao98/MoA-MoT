def solve_inertial_quotient_order():
    """
    Calculates the highest possible order for the inertial quotient E.
    """
    n = 4
    q = 2

    print("The defect group D is elementary abelian of order 16, which means it can be viewed as a 4-dimensional vector space over the field with 2 elements, F_2.")
    print("The highest possible order for the inertial quotient E is the order of the automorphism group of D, which is the general linear group GL(4, 2).")
    print("\nThe order of GL(n, q) is calculated by the product: (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))")
    print(f"For n={n} and q={q}, this is:")

    # Build the expression string
    expression_str = ""
    terms = []
    order = 1
    for i in range(n):
        term = (q**n - q**i)
        terms.append(term)
        order *= term
        expression_str += f"({q}^{n} - {q}^{i})"
        if i < n - 1:
            expression_str += " * "
    print(expression_str)

    # Build the values string
    values_str = " = "
    for i, term in enumerate(terms):
        values_str += str(term)
        if i < len(terms) - 1:
            values_str += " * "
    print(values_str)
    
    print(f"\nThe highest order that E can have is: {order}")

solve_inertial_quotient_order()
<<<20160>>>