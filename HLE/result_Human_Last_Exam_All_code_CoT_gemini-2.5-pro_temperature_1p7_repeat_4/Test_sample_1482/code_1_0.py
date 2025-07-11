def generate_it3_mf_formulation():
    """
    This function generates and prints the mathematical formulation
    for the upper bound of a vertical cross-section in a
    Gaussian-based Interval Type-3 Membership Function (IT3 MF).

    The formulation describes how the membership grade (μ) depends on the
    primary input variable 'x' and the secondary variable 'u'.

    Symbols used:
    μ_upper(x, u): The upper bound of the IT3 MF.
    exp(...): The exponential function, e^...
    u: The secondary variable, controlling the level of uncertainty.
    x: The primary input variable.
    c: The center (mean) of the Gaussian function.
    σ_upper: The standard deviation defining the upper bound of the FOU.
    """

    # The formula is constructed as a string. Unicode characters are used for clarity.
    # The upper bound of the Footprint of Uncertainty (FOU) is a Type-1 Gaussian:
    # FOU_upper(x) = exp( -0.5 * ( (x - c) / σ_upper )^2 )
    # The IT3 MF upper bound is derived from the FOU upper bound:
    # μ_upper(x, u) = [ FOU_upper(x) ]^u
    # This simplifies to the expression below.
    equation = "μ_upper(x, u) = exp( -0.5 * u * ( (x - c) / σ_upper )^2 )"

    print("The mathematical formulation for the upper bound of a vertical cross-section of a Gaussian IT3 MF is:")
    print(equation)

generate_it3_mf_formulation()
