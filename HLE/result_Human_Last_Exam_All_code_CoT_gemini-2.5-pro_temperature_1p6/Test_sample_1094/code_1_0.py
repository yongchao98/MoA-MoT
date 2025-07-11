def solve_ac_loss():
    """
    This function provides the formula for the normalized AC loss in a superconducting
    elliptic bar with a transport current less than the critical current (i < 1).
    """

    # The normalized loss is represented by the left-hand side (LHS) of the equation.
    lhs = "2 * pi * Q / (mu_0 * Ic**2)"

    # The formula is a function of i = Im/Ic. We construct the right-hand side (RHS).
    # The numbers in the equation are 2, 1, 1, and the exponent 2.
    
    # Term 1 of the formula: (1-i)*ln(1-i)
    term1 = "(1 - i) * ln(1 - i)"
    
    # Term 2 of the formula: (1+i)*ln(1+i)
    term2 = "(1 + i) * ln(1 + i)"

    # Term 3 of the formula: i^2
    term3 = "i**2"

    # The full right-hand side of the equation
    rhs = f"2 * [{term1} + {term2} - {term3}]"

    # Print the final result in a clear format
    print("The normalized AC loss per cycle per unit length is given by the Norris equation for an elliptical wire.")
    print("The equation is independent of the bar's aspect ratio (a/b) for i < 1.")
    print("-" * 50)
    print(f"{lhs} = {rhs}")
    print("-" * 50)
    print("\nWhere:")
    print("  Q: Loss per cycle per unit length.")
    print("  mu_0: Permeability of free space.")
    print("  pi: The mathematical constant pi (approx. 3.14159).")
    print("  Ic: The critical current of the bar (Ic = Jc * pi * a * b).")
    print("  i: The normalized current amplitude (i = Im / Ic).")
    print("  ln: The natural logarithm.")

# Execute the function to display the result
solve_ac_loss()