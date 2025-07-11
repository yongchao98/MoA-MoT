def solve():
    """
    Calculates the number of internal adjunctions in the specified 2-category
    from F_11^3 to itself.
    """
    n = 3
    q = 11

    # The number of adjunctions is the order of the general linear group GL(n, q).
    # The formula is |GL(n, q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1)).

    q_n = q**n

    terms = []
    for i in range(n):
        terms.append(q_n - q**i)

    result = 1
    for term in terms:
        result *= term
        
    print(f"The number of internal adjunctions is the order of the general linear group GL({n}, {q}).")
    print("The calculation is as follows:")

    # Build the equation string step-by-step
    # First part: Formula with powers
    q_n_expr = f"{q}^{n}"
    power_exprs = [f"({q_n_expr} - {q}^{i})" for i in range(n)]
    print(" * ".join(power_exprs))
    
    # Second part: Evaluated powers
    q_n_val = q**n
    q_pow_vals = [q**i for i in range(n)]
    val_exprs = [f"({q_n_val} - {val})" for val in q_pow_vals]
    print("= " + " * ".join(val_exprs))
    
    # Third part: Evaluated terms
    term_vals = [(q_n_val - val) for val in q_pow_vals]
    term_exprs = [str(t) for t in term_vals]
    print("= " + " * ".join(term_exprs))

    # Final result
    final_answer = 1
    for t in term_vals:
        final_answer *= t
    
    print(f"= {final_answer}")

solve()