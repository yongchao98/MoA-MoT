def solve_bch_order_10():
    """
    Calculates the number of nonzero coefficients of order 10 in the
    Baker-Campbell-Hausdorff (BCH) expansion using Witt's formula.
    """
    n = 10  # The order
    k = 2   # The number of generators (X and Y)

    # For n=10, the divisors 'd' are 1, 2, 5, 10.
    divisors = [1, 2, 5, 10]

    # Pre-calculated values of the Mobius function for the divisors of 10.
    mu_values = {1: 1, 2: -1, 5: -1, 10: 1}

    equation_terms_str = []
    value_terms = []
    total_sum = 0

    for d in divisors:
        mu_d = mu_values[d]
        power = n // d
        term_value = mu_d * (k ** power)

        total_sum += term_value
        value_terms.append(str(term_value))
        
        # Build the string for each part of the equation
        if mu_d == 1:
            # hide the '1*' for cleaner output
            term_str = f"{k}^{power}"
        else:
            # show '(-1)*...'
            term_str = f"({mu_d}) * {k}^{power}"
        equation_terms_str.append(term_str)
        
    final_result = total_sum // n

    # Join terms for printing the equation. Replace '+ -' with '-'
    final_equation = f" ( {' + '.join(equation_terms_str)} )".replace('+ (-', '- (')
    final_values_eq = f" ( {' + '.join(value_terms)} )".replace('+ -', '- ')

    print("The number of nonzero coefficients is calculated using Witt's formula:")
    print("L_n(k) = (1/n) * sum_{d|n} [mu(d) * k^(n/d)]")
    print("\nFor order n=10 and k=2 generators:")
    # The user wanted to see each number in the final equation.
    print(f"L_10(2) = (1/{n}) *{final_equation}")
    print(f"         = (1/{n}) *{final_values_eq}")
    print(f"         = (1/{n}) * ( {total_sum} )")
    print(f"         = {final_result}")

solve_bch_order_10()