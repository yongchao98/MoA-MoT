import sympy

def solve_imitation_learning_bound():
    """
    This function calculates and prints the tightest upper bound for the performance
    difference J(pi*) - J(pi_hat) based on the provided problem description.
    """
    # Define the symbolic variables from the problem statement
    H = sympy.Symbol('H')  # Horizon
    A_size = sympy.Symbol('|A|')  # Size of the action space
    lam = sympy.Symbol('lambda')  # Hyperparameter lambda

    # The standard compounding error bound from imitation learning theory is:
    # J(pi*) - J(pi_hat) <= H * (H - 1) * R_max * epsilon_TV
    # where epsilon_TV is the average Total Variation error. We assume R_max = 1.
    # The population TV risk given is interpreted as the expected L1 distance,
    # which is equal to 2 * epsilon_TV.
    # T_risk = 2 * epsilon_TV  =>  epsilon_TV = T_risk / 2.
    # So, J(pi*) - J(pi_hat) <= H * (H - 1) * (T_risk / 2).
    # We are given T_risk <= |A| * (1 - exp(-lambda)).
    # Substituting this gives the final bound.

    # Calculate the coefficient for the expression
    coefficient = sympy.Rational(1, 2)

    # The rest of the expression
    expression = H * (H - 1) * A_size * (1 - sympy.exp(-lam))

    # The final upper bound formula
    upper_bound = coefficient * expression

    # We use pretty print to display the mathematical formula clearly
    print("The tightest upper bound for J(pi*) - J(pi_hat) is:")
    sympy.pprint(upper_bound, use_unicode=True)

    # The instruction asks to output each number in the final equation.
    # In the expression 0.5 * H * (H-1) * |A| * (1 - exp(-lambda)), the numbers are:
    print("\nThe numbers in the final equation are:")
    print(f"Coefficient: {float(coefficient)}")
    print("From (H-1): 1")
    print("From (1 - ...): 1")
    print("In the exponent exp(-lambda): -1")

if __name__ == "__main__":
    solve_imitation_learning_bound()