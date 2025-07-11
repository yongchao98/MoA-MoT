def generate_it3_mf_formulation():
    """
    Generates and explains the mathematical formulation for the vertical
    cross-sections of an Interval Type-3 Membership Function (IT3 MF)
    using a Gaussian-based paradigm.
    """

    print("Interval Type-3 Fuzzy Logic: Vertical Cross-Section Formulation")
    print("=" * 65)
    print("\nThis script presents the general mathematical formula for a bound (upper or lower)")
    print("of a vertical cross-section in an Interval Type-3 Membership Function (IT3 MF).")
    print("The formulation uses a Gaussian function, as specified.\n")

    # Define the components of the equation for clarity
    bound_symbol = "μ_bound(u; x)"
    scale_func = "scale(x)"
    primary_var = "x"
    secondary_var = "u"
    center_func = "c(x)"
    std_dev_func = "σ(x)"
    exponent_term = f"(({secondary_var} - {center_func}) / {std_dev_func})"

    # Construct the final equation string
    # Using 'exp' for the exponential function
    equation = f"{bound_symbol} = {scale_func} * exp( -0.5 * {exponent_term}² )"
    
    # Using Unicode characters for better display
    equation_pretty = f"μ_bound(u; x) = scale(x) · exp( -0.5 · ( (u - c(x)) / σ(x) )² )"


    print("The general formula is:\n")
    # Print each component of the equation clearly
    print(f"  μ_bound(u; x) = scale(x) * exp( -0.5 * ( (u - c(x)) / σ(x) )^2 )\n")
    
    print("Where each component represents:")
    print("-" * 35)
    print(f"  {bound_symbol:<15}: The value of the membership bound (either upper or lower) for the vertical cross-section.")
    print(f"  {secondary_var:<15}: The secondary input variable, which is the domain of the vertical cross-section.")
    print(f"  {primary_var:<15}: The primary input variable, which is fixed for the cross-section and acts as a parameter.")
    print(f"  {scale_func:<15}: A scaling function that determines the height (peak) of the Gaussian. Its value depends on '{primary_var}'.")
    print(f"  {center_func:<15}: The center (mean) of the Gaussian function. Its position depends on '{primary_var}'.")
    print(f"  {std_dev_func:<15}: The standard deviation (width) of the Gaussian. Its value, representing uncertainty, depends on '{primary_var}'.")
    print(f"  exp(...):<14} The standard exponential function.")
    print(f"  -0.5:<14} A constant factor characteristic of the Gaussian formulation.")

if __name__ == "__main__":
    generate_it3_mf_formulation()