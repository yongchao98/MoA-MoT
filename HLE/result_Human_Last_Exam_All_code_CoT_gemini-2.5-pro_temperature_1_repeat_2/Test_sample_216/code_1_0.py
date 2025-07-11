import sympy

def solve_imitation_learning_bound():
    """
    This function formulates and prints the tightest upper bound for the performance
    difference between an expert and a learned policy in imitation learning.
    """
    # Define the symbols used in the equation
    H = sympy.Symbol('H')  # Horizon
    A = sympy.Symbol('|A|')  # Size of the action space
    lmbda = sympy.Symbol('lambda')  # Hyperparameter lambda

    # The performance difference J(pi*) - J(pi_hat) is bounded by a factor
    # that depends quadratically on the horizon H, due to compounding error.
    # This factor is H*(H+1)/2.
    compounding_error_factor = H * (H + 1) / 2

    # The problem provides an upper bound on the population total variation (TV) risk.
    tv_risk_bound = A * (1 - sympy.exp(-lmbda))

    # The final upper bound on the performance difference is the product of these two terms.
    upper_bound = compounding_error_factor * tv_risk_bound

    # Print the final expression for the upper bound.
    # We use sympy's pretty print for a clear mathematical representation.
    print("The tightest upper bound of J(pi^*) - J(pi_hat) is given by the expression:")
    sympy.pprint(upper_bound, use_unicode=False)
    
    # For a machine-readable format as requested
    print("\nIn simple text format, the equation is:")
    print(f"J(pi^*) - J(pi_hat) <= ({H}*({H}+1)/2) * {A} * (1 - exp(-{lmbda}))")
    
    # As per instructions to output each 'number' (symbol in this case) in the final equation.
    # We construct the string representation of the final formula piece by piece.
    final_equation_str = f"H * (H + 1) / 2 * |A| * (1 - exp(-lambda))"
    
    # Let's consider the final answer format requested
    # The final answer is the formula itself.
    print("\nFinal Answer Formula:")
    print(final_equation_str)


solve_imitation_learning_bound()