def formulate_it3_mf_vertical_cross_section():
    """
    This function formulates and prints the mathematical equation for the upper bound
    of a vertical cross-section in a Gaussian Interval Type-3 Membership Function.
    """
    
    # --- Symbolic Representation of Variables and Parameters ---
    # These strings represent the components of our mathematical equation.
    
    # The secondary input variable, defined on the vertical axis of the primary membership interval.
    u = "u"
    
    # The primary input variable.
    x = "x"
    
    # The center (or mean) of the primary Gaussian membership function for x.
    c_x = "c_x"
    
    # The standard deviation of the primary Upper Membership Function (UMF).
    sigma_bar_x = "σ_bar_x"
    
    # A positive constant parameter that controls the amount of "blurriness" in the third dimension.
    k = "k"

    # --- Constructing the Equation ---
    
    # 1. Formulate the primary Upper Membership Function (UMF), μ_bar_A(x), as a Gaussian.
    # This defines the upper bound of the membership interval for the primary input x.
    primary_umf_expr = f"exp(-0.5 * (({x} - {c_x}) / {sigma_bar_x})**2)"
    
    # 2. Formulate the complete equation for the upper bound of the vertical cross-section.
    # This is a Gaussian function of 'u', centered at the value of the primary UMF,
    # with a standard deviation proportional to the primary UMF.
    final_equation = (
        f"μ_bar_A(x)({u}) = exp(-0.5 * ((({u}) - ({primary_umf_expr})) / (({k}) * ({primary_umf_expr})))**2)"
    )

    print("The mathematical formulation for the upper bound of a Gaussian IT3 MF vertical cross-section is:")
    print("-" * 80)
    print(final_equation)
    print("-" * 80)
    print("Where:")
    print(f"  • μ_bar_A(x)({u}) is the upper membership grade for the secondary variable '{u}' at a given primary input '{x}'.")
    print(f"  • '{u}' is the secondary input variable.")
    print(f"  • '{x}' is the primary input variable.")
    print(f"  • '{c_x}' is the center of the primary Gaussian function.")
    print(f"  • '{sigma_bar_x}' is the standard deviation for the primary upper membership function.")
    print(f"  • '{k}' is a parameter controlling the blurriness of the vertical slice.")

# Execute the function to print the result.
formulate_it3_mf_vertical_cross_section()