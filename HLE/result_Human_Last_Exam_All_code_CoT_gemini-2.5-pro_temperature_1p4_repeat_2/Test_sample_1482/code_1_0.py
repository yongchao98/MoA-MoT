import math

def generate_it3_formulation():
    """
    This function generates and prints the mathematical formulation
    for the vertical cross-section of a Gaussian-based IT3 MF.
    """

    # --- Step 1: Define Parameters for the Gaussian models ---
    # These are example parameters for a concrete formulation.

    # Parameters for the primary membership function which defines the center of the vertical slice
    # h_p: Height of the primary Gaussian (typically 1.0)
    # c_p: Center (mean) of the primary Gaussian on the 'x' axis
    # sigma_p: Standard deviation of the primary Gaussian on the 'x' axis
    h_p = 1.0
    c_p = 5.0
    sigma_p = 2.0

    # Parameter for the uncertainty model
    # k_sigma: A scaling factor that determines the width of the uncertainty (std. dev. of the vertical slice)
    # relative to the primary membership value.
    k_sigma = 0.25

    # --- Step 2: Construct the string representations of the component functions ---

    # The center of the Gaussian vertical slice, mu_c(x), is defined by a primary Gaussian MF.
    # This function determines the "ridge" or peak location of the membership over 'u' for a given 'x'.
    mu_c_x_str = f"{h_p} * exp(-0.5 * ((x - {c_p}) / {sigma_p})^2)"

    # The standard deviation of the Gaussian vertical slice, sigma_u(x), models the uncertainty.
    # Here, we model it as being proportional to the center value mu_c(x).
    sigma_u_x_str = f"{k_sigma} * mu_c(x)"
    sigma_u_x_full_str = f"{k_sigma} * ({mu_c_x_str})"

    # --- Step 3: Construct and print the final formulation ---

    print("In the context of an Interval Type-3 Membership Function (IT3 MF), a vertical cross-section at a fixed primary input 'x' is an Interval Type-2 MF over the secondary variable 'u'.")
    print("The upper bound of this vertical cross-section's Footprint of Uncertainty, overline{mu}_x(u), can be modeled with a Gaussian function as follows:")
    print("\n[--- General Formulation ---]\n")
    print("  overline{mu}_x(u) = exp(-0.5 * ((u - mu_c(x)) / sigma_u(x))^2)")
    print("\nWhere the parameters are functions of 'x':")

    print("\n  1. The center, mu_c(x), is defined by a primary Gaussian function:")
    print(f"     mu_c(x) = {mu_c_x_str}")

    print("\n  2. The standard deviation, sigma_u(x), models the uncertainty:")
    print(f"     sigma_u(x) = {sigma_u_x_str}")

    # --- Step 4: Assemble the complete equation with all numbers ---
    
    final_equation = f"overline{{mu}}_x(u) = exp(-0.5 * ((u - ({mu_c_x_str})) / ({sigma_u_x_full_str}))^2)"

    print("\n\n[--- Complete Formulation with Numerical Parameters ---]\n")
    print("Substituting the component functions yields the final detailed equation:")
    print(f"  {final_equation}")


# Execute the function to print the formulation
generate_it3_formulation()

# The final answer as a single string
final_answer_string = f"overline{{mu}}_x(u) = exp(-0.5 * ((u - (1.0 * exp(-0.5 * ((x - 5.0) / 2.0)^2))) / (0.25 * (1.0 * exp(-0.5 * ((x - 5.0) / 2.0)^2))))^2)"
print(f"\n\n<<<" + final_answer_string + ">>>")