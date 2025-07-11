def solve_superconductor_field():
    """
    This function generates and prints the mathematical expression for the magnetic field
    outside a stack of superconducting strips.
    """
    # Define the components of the formula as strings for readability.

    # The cosine term relates the flux front 'a' to the applied field 'H_a'.
    cos_term = "cos^2(H_a / H_0)"

    # The hyperbolic sine terms arise from the conformal mapping for the stack.
    sinh_x_term = "sinh^2(pi * x / D)"
    sinh_w_term = "sinh^2(pi * w / D)"

    # The numerator of the argument of the logarithm.
    numerator = f"({sinh_x_term} - {sinh_w_term})"

    # The denominator of the argument of the logarithm.
    denominator = f"({sinh_x_term} - {sinh_w_term} * {cos_term})"

    # The complete induced field part of the equation.
    # Note: Jc*d/(2*pi) is equal to H_0/2
    induced_field = f"(H_0 / 2) * ln[{numerator} / {denominator}]"

    # The total field is the sum of the applied field and the induced field.
    total_field_expression = f"H_z(x) = H_a + {induced_field}"

    print("The expression for the magnetic field H_z(x) for |x| > w is:")
    print(total_field_expression)

solve_superconductor_field()