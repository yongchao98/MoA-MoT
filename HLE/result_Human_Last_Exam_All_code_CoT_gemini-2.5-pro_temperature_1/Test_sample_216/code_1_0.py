def solve():
    """
    Calculates and prints the tightest upper bound for the performance difference
    J(pi^*) - J(hat_pi) based on the provided information.
    """
    # The variables are symbolic: H, |A|, lambda.
    # We will represent the formula as a string.
    
    # The tightest upper bound on the performance difference J(pi^*) - J(hat_pi)
    # is derived from the standard analysis of imitation learning, which results
    # in a quadratic dependency on the horizon H.
    
    # The bound is H^2 * R_max * T(hat_pi, pi^*).
    # Assuming rewards are in [0, 1], R_max = 1.
    
    # The given bound for the total variation (TV) risk T(hat_pi, pi^*) is:
    # |A| * (1 - e^(-lambda))
    
    # Therefore, the final bound is H^2 * |A| * (1 - e^(-lambda)).
    
    # We will print the components of this expression.
    term_H_sq = "H^2"
    term_A = "|A|"
    term_lambda = "(1 - exp(-lambda))"
    
    final_equation = f"{term_H_sq} * {term_A} * {term_lambda}"
    
    print("The tightest upper bound on J(pi^*) - J(hat_pi) is given by the expression:")
    print(final_equation)

solve()