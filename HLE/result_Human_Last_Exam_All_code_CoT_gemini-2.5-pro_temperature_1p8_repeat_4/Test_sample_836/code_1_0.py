def print_dynamic_structure_factor_equation():
    """
    Prints the name and defining equation for the space-time, double Fourier
    transform of the generalized pair correlation function.
    """
    print("The space-time, double Fourier transform of the generalized pair correlation function G(r, t) is called the Dynamic Structure Factor, denoted as S(k, ω).")
    print("\nThe defining equation is presented below. Note that '∫d³r' denotes a volume integral over all space and '∫dt' denotes an integral over all time.")
    print("-" * 50)

    # Define the components of the equation as strings for clear printing.
    result_function = "S(k, ω)"
    equals_sign = "="
    time_integral = "∫dt"
    space_integral = "∫d³r"
    correlation_function = "G(r, t)"
    operator = "*"
    exponential_function = "exp"
    left_paren = "("
    exponent = "-i(k·r - ωt)" # This contains the imaginary unit 'i' and the implied coefficient '-1'
    right_paren = ")"

    # Print the full equation by assembling the components.
    # S(k, ω) = ∫dt ∫d³r G(r, t) * exp(-i(k·r - ωt))
    print(f"Final Equation: {result_function} {equals_sign} {time_integral} {space_integral} {correlation_function} {operator} {exponential_function}{left_paren}{exponent}{right_paren}")

    # To satisfy the instruction to output each number/symbol in the equation:
    print("\n--- Key Symbols in the Equation ---")
    print(f"S, k, ω : The symbols defining the final function, the Dynamic Structure Factor S, as a function of wavevector k and angular frequency ω.")
    print(f"G, r, t : The symbols defining the original function, the generalized pair correlation function G, as a function of spatial separation r and time difference t.")
    print(f"exp: The exponential function.")
    print(f"-1: The implicit numerical coefficient in the exponent.")
    print(f"i: The imaginary unit, where i² = -1.")
    print("-" * 50)

# Execute the function to display the information.
print_dynamic_structure_factor_equation()