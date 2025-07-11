import math

def gaussian(x, c, sigma):
    """Calculates the value of a Gaussian membership function."""
    if sigma == 0:
        return 1.0 if x == c else 0.0
    exponent = -0.5 * ((x - c) / sigma)**2
    return math.exp(exponent)

def calculate_vertical_cross_section():
    """
    Calculates and prints the formulation for the vertical cross-section of an IT3 MF.
    """
    # Parameters for the primary Gaussian bounds
    c_primary = 5.0      # Center of the primary Gaussian functions
    sigma_upper = 1.5    # Std deviation for the upper bound (¯μ_U)
    sigma_lower = 1.0    # Std deviation for the lower bound (μ_U)

    # Scaling constant for the secondary MF's standard deviation
    k = 0.5

    # Fixed input values for demonstration
    x_input = 4.5  # Primary input variable
    u_input = 0.9  # Secondary input variable

    # --- Step 1: Calculate the primary bounds at x_input ---
    # The upper bound mu_upper must be >= the lower bound mu_lower.
    # For a Gaussian MF, a larger sigma gives a wider curve, so for x != c,
    # the function with the larger sigma will have a higher value.
    mu_upper_val = gaussian(x_input, c_primary, sigma_upper)
    mu_lower_val = gaussian(x_input, c_primary, sigma_lower)

    # --- Step 2: Calculate parameters for the secondary (vertical) MF ---
    # Check if bounds are distinct to avoid division by zero
    if mu_upper_val == mu_lower_val:
        print("Primary bounds are identical; vertical cross-section is undefined.")
        return

    c_secondary = (mu_upper_val + mu_lower_val) / 2
    sigma_secondary = k * (mu_upper_val - mu_lower_val)

    # --- Step 3: Calculate the final vertical cross-section membership value ---
    result = gaussian(u_input, c_secondary, sigma_secondary)
    
    # --- Step 4: Print the detailed formulation with the calculated numbers ---
    print("--- Mathematical Formulation for the Vertical Cross-Section ---")
    print(f"This formulation models the secondary membership function μ_VCS(u; x) as a Gaussian function.")
    print("The primary bounds are themselves Gaussians defined by a center (c) and standard deviations (σ_U, σ_L).")
    print("\n--- Input Parameters ---")
    print(f"Primary MF Center (c): {c_primary}")
    print(f"Primary Upper Bound Sigma (σ_U): {sigma_upper}")
    print(f"Primary Lower Bound Sigma (σ_L): {sigma_lower}")
    print(f"Secondary MF Width Constant (k): {k}")
    print(f"Primary Input (x): {x_input}")
    print(f"Secondary Input (u): {u_input}")

    print("\n--- Calculated Primary Bounds at x = {x_input} ---")
    print(f"Upper Primary Bound ¯μ_U(x): exp(-0.5 * (({x_input} - {c_primary}) / {sigma_upper})²) = {mu_upper_val:.4f}")
    print(f"Lower Primary Bound μ_U(x):  exp(-0.5 * (({x_input} - {c_primary}) / {sigma_lower})²) = {mu_lower_val:.4f}")

    print("\n--- Calculated Secondary MF Parameters at x = {x_input} ---")
    print(f"Secondary Center c_s(x): ({mu_upper_val:.4f} + {mu_lower_val:.4f}) / 2 = {c_secondary:.4f}")
    print(f"Secondary Sigma σ_s(x): {k} * ({mu_upper_val:.4f} - {mu_lower_val:.4f}) = {sigma_secondary:.4f}")

    print("\n--- Final Equation for μ_VCS(u={u_input}; x={x_input}) ---")
    final_equation = f"μ_VCS = exp(-0.5 * (({u_input} - {c_secondary:.4f}) / {sigma_secondary:.4f})²)"
    print(final_equation)
    print(f"\nResult: {result:.4f}")
    
    # The final answer to be returned in the specified format
    global final_result_value
    final_result_value = result

# Execute the calculation and store the final answer
final_result_value = None
calculate_vertical_cross_section()