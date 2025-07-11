def formulate_it3_vertical_section():
    """
    This function formulates and prints the mathematical equation for the
    upper bound of a vertical cross-section of a Gaussian Interval Type-3
    Membership Function (IT3 MF).
    """

    print("Formulation for the Vertical Cross-Section of a Gaussian IT3 MF")
    print("="*65)
    print("This formulation describes the upper bound of the membership grade for a secondary")
    print("variable 'u', given a fixed primary input 'x'.\n")

    print("--- Symbol Definitions ---")
    print(" Ã: The Interval Type-3 fuzzy set.")
    print(" x: The primary input variable.")
    print(" u: The secondary variable, representing the membership grade of x in the primary set.")
    print(" μ_P(x): The principal (Type-1) Gaussian membership function of x.")
    print(" c: The center of the principal Gaussian function.")
    print(" σ_x: The standard deviation (spread) of the principal Gaussian function.")
    print(" overline{μ}_{Ã(x)}(u): The Upper Membership Function (UMF) of the vertical cross-section at x.")
    print("                     This function defines the upper bound of the uncertainty.")
    print(" overline{β}: A constant parameter that scales the uncertainty for the UMF. A larger")
    print("             value indicates greater uncertainty.")
    print("\n--- Component Formulas ---")

    # Formula for the Principal Membership Function
    principal_mf = "μ_P(x) = exp(-0.5 * ((x - c) / σ_x)²) "
    print("1. Principal Membership Function (T1 MF):")
    print(f"   {principal_mf}\n")

    # Assembling the final formula for the vertical cross-section's upper bound
    # Note: Unicode characters are used for clarity (μ, σ, β, Ã, ², overline).
    final_equation = "\u0305μ_{\u00c3(x)}(u) = exp(-0.5 * ((u - μ_P(x)) / (\u0305β * (1 - μ_P(x))))²)"

    print("2. Final Formulation for the Upper Bound of the Vertical Cross-Section:")
    print("   The vertical cross-section at 'x' is an Interval Type-2 MF whose upper bound is")
    print("   defined by the following Gaussian function of 'u':")
    print(f"\n   {final_equation}")
    print("\n" + "="*65)
    print("In this equation, the term '-0.5' is the standard coefficient for a Gaussian function.")
    print("The center of this Gaussian is 'μ_P(x)' and its standard deviation is 'β_upper * (1 - μ_P(x))'.")

# Execute the function to print the formulation
formulate_it3_vertical_section()