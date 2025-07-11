def generate_it3_vertical_slice_formula():
    """
    This function generates and prints the mathematical formulation for the
    upper bound of a vertical cross-section of an Interval Type-3
    Membership Function (IT3 MF) using a Gaussian-based paradigm.

    The formula describes the membership grade as a function of the
    secondary variable 'u' for a fixed primary variable 'x'.
    """

    # Define the components of the formula as strings for display purposes
    function_name = "f_upper(u; x)"
    primary_variable = "x"
    secondary_variable = "u"
    upper_mf = f"μ_U({primary_variable})"
    scale_parameter = "α_U"

    # The mean of the Gaussian function is the upper membership function of the underlying IT2 set
    mean = upper_mf
    # The standard deviation is proportional to the mean, scaled by a parameter alpha
    std_dev = f"({scale_parameter} * {upper_mf})"

    # Construct the final formula string using the components
    # The formula is for a Gaussian function: exp(-0.5 * ((var - mean) / std_dev)^2)
    exponent_numerator = f"({secondary_variable} - {mean})"
    exponent_term = f"({exponent_numerator} / {std_dev})^2"
    # The constant -0.5 is explicitly included as requested
    equation = f"{function_name} = exp(-0.5 * {exponent_term})"

    print("The mathematical formulation for the upper bound of the vertical cross-section of an IT3 MF is:")
    print(equation)


if __name__ == "__main__":
    generate_it3_vertical_slice_formula()