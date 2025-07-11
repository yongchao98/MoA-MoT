import math

def generate_it3_formulation():
    """
    Determines and prints the mathematical formulation for the vertical cross-section
    of a sample Interval Type-3 Membership Function (IT3 MF).
    """

    # --- 1. Define Model Parameters ---

    # Fixed primary input variable
    x_prime = 6.0

    # Secondary input variable in [0, 1], representing the level of uncertainty
    u = 0.5

    # Parameters for the base Gaussian function G(x) = exp(-0.5 * ((x - c) / s)^2)
    # This function provides the fundamental shape.
    gaussian_mean_c = 5.0
    gaussian_std_dev_sigma = 3.0

    # Scaling factors for the four bounding membership functions.
    # These factors must be ordered s_ll <= s_ul <= s_lu <= s_uu to ensure
    # the membership functions are correctly nested.
    s_uu = 1.0  # Scale for Upper bound of the Upper IT2 MF
    s_lu = 0.9  # Scale for Lower bound of the Upper IT2 MF
    s_ul = 0.7  # Scale for Upper bound of the Lower IT2 MF
    s_ll = 0.5  # Scale for Lower bound of the Lower IT2 MF

    # --- 2. Define Core Functions ---

    def gaussian(x, mean, std_dev):
        """Calculates the value of a Gaussian function."""
        return math.exp(-0.5 * ((x - mean) / std_dev) ** 2)

    # --- 3. Perform Calculations ---

    # Calculate the value of the base Gaussian shape function at x_prime
    G_at_x_prime = gaussian(x_prime, gaussian_mean_c, gaussian_std_dev_sigma)

    # Calculate the values of the four bounding functions at x_prime by scaling G
    mu_upper_U = s_uu * G_at_x_prime
    mu_lower_U = s_lu * G_at_x_prime
    mu_upper_L = s_ul * G_at_x_prime
    mu_lower_L = s_ll * G_at_x_prime

    # Use the linear model to calculate the bounds of the vertical cross-section
    U_bound = u * mu_upper_U + (1 - u) * mu_lower_U
    L_bound = u * mu_upper_L + (1 - u) * mu_lower_L

    # --- 4. Print the Formulation and Results ---

    print("Mathematical Formulation for an IT3 MF Vertical Cross-Section\n")
    print("The vertical cross-section for a fixed primary input x' and a secondary input u is an interval [L(u), U(u)].")
    print("This formulation uses a base Gaussian G(x) and four scaling factors to define bounding IT2 MFs.\n")

    print(f"Parameters:")
    print(f"  Fixed primary input, x' = {x_prime}")
    print(f"  Secondary input, u = {u}")
    print(f"  Base Gaussian G(x) = exp(-0.5 * ((x - {gaussian_mean_c}) / {gaussian_std_dev_sigma})^2)")
    print("-" * 60)

    print("Step 1: Calculate the value of the four bounding functions at x'.")
    print(f"  μ_upper_U(x') = s_uu * G(x') = {s_uu:.2f} * {G_at_x_prime:.4f} = {mu_upper_U:.4f}")
    print(f"  μ_lower_U(x') = s_lu * G(x') = {s_lu:.2f} * {G_at_x_prime:.4f} = {mu_lower_U:.4f}")
    print(f"  μ_upper_L(x') = s_ul * G(x') = {s_ul:.2f} * {G_at_x_prime:.4f} = {mu_upper_L:.4f}")
    print(f"  μ_lower_L(x') = s_ll * G(x') = {s_ll:.2f} * {G_at_x_prime:.4f} = {mu_lower_L:.4f}")
    print("-" * 60)

    print("Step 2: Formulate and compute the final Upper Bound U(u).")
    print(f"  U(u) = u * μ_upper_U(x') + (1 - u) * μ_lower_U(x')")
    print(f"  U({u}) = {u} * {mu_upper_U:.4f} + (1 - {u}) * {mu_lower_U:.4f}")
    print(f"       = {u * mu_upper_U:.4f} + {(1 - u) * mu_lower_U:.4f}")
    print(f"       = {U_bound:.4f}\n")

    print("Step 3: Formulate and compute the final Lower Bound L(u).")
    print(f"  L(u) = u * μ_upper_L(x') + (1 - u) * μ_lower_L(x')")
    print(f"  L({u}) = {u} * {mu_upper_L:.4f} + (1 - {u}) * {mu_lower_L:.4f}")
    print(f"       = {u * mu_upper_L:.4f} + {(1 - u) * mu_lower_L:.4f}")
    print(f"       = {L_bound:.4f}")
    print("-" * 60)

    print(f"Resulting Interval of the Vertical Cross-Section at x'={x_prime}, u={u}:")
    print(f"[{L_bound:.4f}, {U_bound:.4f}]")

# Execute the function to print the solution
generate_it3_formulation()
