import sympy

def solve_optimal_feedback():
    """
    Calculates the optimal feedback adjustment factor 'b' for a given
    first transmission power 'P1' and noise correlation 'alpha'.

    The optimal value is derived by maximizing the channel's mutual information,
    which simplifies to maximizing the determinant of the received signal's
    covariance matrix. The result of this optimization is b = -alpha * P1.
    """
    # Example values for the parameters
    # alpha is the weather-induced correlation
    alpha = 0.6
    # P1 is the power used in the first transmission
    P1 = 10.0

    # Calculate the optimal feedback factor b using the derived formula
    b = -alpha * P1

    # Print the parameters and the result
    print(f"Given parameters:")
    print(f"  Noise correlation (alpha): {alpha}")
    print(f"  Power in first transmission (P1): {P1}\n")
    print(f"The optimal feedback adjustment factor (b) is: {b:.4f}\n")

    # Display the final equation with the numbers
    print("The final equation relating these values is:")
    # Using sympy to pretty print the equation
    b_sym, alpha_sym, p1_sym = sympy.symbols('b alpha P1')
    equation = sympy.Eq(b_sym, -alpha_sym * p1_sym)
    # Substitute numerical values
    numerical_equation = equation.subs({b_sym: b, alpha_sym: alpha, p1_sym: P1})
    
    # We are asked to output each number in the final equation.
    # The equation is b = -alpha * P1
    # We will print it in a formatted string.
    print(f"{b:.4f} = -{alpha} * {P1}")

if __name__ == "__main__":
    solve_optimal_feedback()