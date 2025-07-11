import math

def generate_it3_mf_formulation():
    """
    Generates and prints the mathematical formulation for the vertical cross-section
    of a Gaussian-based Interval Type-3 Membership Function (IT3 MF).
    """

    print("This script outlines the mathematical formulation that characterizes the vertical")
    print("cross-section of an upper bound for a Gaussian Interval Type-3 Membership Function (IT3 MF).")
    print("The formulation defines the membership grade as a function of the secondary variable 'u' for a fixed primary input 'x'.")
    print("-" * 70)

    # --- Main Equation ---
    main_equation = "μ_upper_vs(u | x) = exp(-0.5 * ( (u - C(x)) / S(x) )^2)"
    print("\nThe general mathematical formulation for the vertical slice is:")
    print(f"  {main_equation}")
    
    # --- Component Definitions ---
    center_equation = "C(x) = exp(-0.5 * ( (x - c_p) / σ_p )^2)"
    spread_equation = "S(x) = α * C(x)"

    print("\nWhere the components are defined as:")
    print(f"\n1.  Center of the vertical slice, C(x):")
    print(f"    This represents the principal membership grade at input x.")
    print(f"    {center_equation}")

    print(f"\n2.  Spread of the vertical slice, S(x):")
    print(f"    This represents the uncertainty or 'blur' in the membership grade at input x.")
    print(f"    {spread_equation}")

    print("-" * 70)

    # --- Variable and Parameter Explanations ---
    print("\nExplanation of Variables and Parameters:")
    print("  μ_upper_vs(u | x) : The membership grade of the upper bound's vertical slice.")
    print("  u                 : The secondary input variable, where u is in the range [0, 1].")
    print("  x                 : The fixed primary input variable.")
    print("  exp()             : The exponential function, e^().")
    print("  c_p               : A symbolic parameter representing the center of the primary Gaussian function.")
    print("  σ_p               : A symbolic parameter representing the standard deviation of the primary Gaussian function.")
    print("  α                 : A symbolic parameter (alpha) that scales the spread of uncertainty.")

    # --- Final Answer Formatting ---
    # The user wants the final equation as the answer. We provide the main formula.
    final_answer = main_equation
    print(f"\n<<<μ_upper_vs(u | x) = exp(-0.5 * ( (u - C(x)) / S(x) )^2)>>>")

# Execute the function to print the formulation
generate_it3_mf_formulation()