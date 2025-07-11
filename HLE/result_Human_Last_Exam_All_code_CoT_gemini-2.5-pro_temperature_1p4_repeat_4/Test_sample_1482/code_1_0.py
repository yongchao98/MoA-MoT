import math

def display_it3_vertical_cross_section_formulation():
    """
    Calculates and displays the mathematical formulation for the vertical
    cross-section of an Interval Type-3 Fuzzy Membership Function (IT3 MF)
    using Gaussian-based paradigms.
    """
    # --- 1. Define Input Parameters ---
    # Primary input variable
    x_val = 5.0
    # Secondary input variable (within the primary membership domain)
    u_val = 0.6

    # Parameters for the primary upper bound Gaussian: μ'_upper(x)
    c_upper = 5.0  # Center
    sigma_upper = 2.0  # Standard deviation

    # Parameters for the primary lower bound Gaussian: μ'_lower(x)
    c_lower = 5.0  # Center
    sigma_lower = 4.0  # Standard deviation
    scale_lower = 0.8  # Scaling factor (height)

    # Parameters for the vertical slice formulation
    # Controls the standard deviation of the vertical slice
    k = 2.0
    # Scales the upper bound to get the lower bound of the vertical slice
    alpha = 0.5

    # --- 2. Print General Formulation and Parameters ---
    print("### Formulation for Vertical Cross-Section of an IT3 MF ###\n")
    print("The vertical cross-section is an interval [μ_lower(x, u), μ_upper(x, u)].")
    print("The general formulation for its bounds is:\n")
    print("  μ_upper(x, u) = exp(-0.5 * ((u - C(x)) / Σ(x))²)")
    print("  μ_lower(x, u) = α * μ_upper(x, u)\n")
    print("Where the parameters C(x) and Σ(x) are derived from the primary membership bounds (the FOU):\n")
    print("  C(x) = (μ'_upper(x) + μ'_lower(x)) / 2")
    print("  Σ(x) = (μ'_upper(x) - μ'_lower(x)) / k\n")
    print("----------------------------------------------------------")
    print(f"Executing calculation for x = {x_val} and u = {u_val}...\n")

    # --- 3. Calculate Primary Membership Bounds (FOU) ---
    mu_prime_upper = math.exp(-0.5 * ((x_val - c_upper) / sigma_upper)**2)
    mu_prime_lower = scale_lower * math.exp(-0.5 * ((x_val - c_lower) / sigma_lower)**2)

    print("Step 1: Calculate the primary membership bounds at x =", x_val)
    print(f"  μ'_upper({x_val}) = exp(-0.5 * (({x_val} - {c_upper}) / {sigma_upper})²) = {mu_prime_upper:.4f}")
    print(f"  μ'_lower({x_val}) = {scale_lower} * exp(-0.5 * (({x_val} - {c_lower}) / {sigma_lower})²) = {mu_prime_lower:.4f}")
    print(f"  Primary FOU at x={x_val} is [{mu_prime_lower:.4f}, {mu_prime_upper:.4f}]\n")

    # Note on 'u' value
    if not (mu_prime_lower <= u_val <= mu_prime_upper):
        print(f"Note: The specified u={u_val} is outside the primary FOU. In a practical system, membership would be 0.")
        print("However, the formulation is calculated here for demonstration.\n")


    # --- 4. Calculate Vertical Slice Parameters ---
    C_x = (mu_prime_upper + mu_prime_lower) / 2
    # Avoid division by zero if bounds are identical
    Sigma_x = (mu_prime_upper - mu_prime_lower) / k if (mu_prime_upper - mu_prime_lower) != 0 else 1e-9


    print("Step 2: Calculate the vertical slice parameters C(x) and Σ(x)")
    print(f"  C({x_val}) = ({mu_prime_upper:.4f} + {mu_prime_lower:.4f}) / 2 = {C_x:.4f}")
    print(f"  Σ({x_val}) = ({mu_prime_upper:.4f} - {mu_prime_lower:.4f}) / {k} = {Sigma_x:.4f}\n")


    # --- 5. Determine and Print the Final Formulation ---
    print("Step 3: Substitute calculated values into the final formulation")

    # Upper Bound
    upper_bound_val = math.exp(-0.5 * ((u_val - C_x) / Sigma_x)**2)
    print("\n--- Upper Bound Formulation ---")
    print(f"μ_upper({x_val}, {u_val}) = exp(-0.5 * (({u_val} - {C_x:.4f}) / {Sigma_x:.4f})²)")
    print(f"                 = {upper_bound_val:.5f}")


    # Lower Bound
    lower_bound_val = alpha * upper_bound_val
    print("\n--- Lower Bound Formulation ---")
    print(f"μ_lower({x_val}, {u_val}) = {alpha} * μ_upper({x_val}, {u_val})")
    print(f"                 = {alpha} * {upper_bound_val:.5f}")
    print(f"                 = {lower_bound_val:.5f}")

    # --- 6. Final Result ---
    print("\n----------------------------------------------------------")
    print("Final Determined Interval for the Vertical Cross-Section:")
    print(f"[μ_lower, μ_upper] = [{lower_bound_val:.5f}, {upper_bound_val:.5f}]")


if __name__ == '__main__':
    display_it3_vertical_cross_section_formulation()
