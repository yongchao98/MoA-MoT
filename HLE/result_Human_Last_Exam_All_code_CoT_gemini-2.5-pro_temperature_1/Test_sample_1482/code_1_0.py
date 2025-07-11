import math

def calculate_it3_vertical_slice_membership():
    """
    Calculates and explains the formulation for the upper bound of a vertical
    cross-section of an Interval Type-3 Membership Function (IT3 MF).

    The formulation is based on a reference IT2 MF with Gaussian bounds.
    """

    # --- Parameters Definition ---
    # These parameters define the shape of the reference IT2 MF and the IT3 MF.
    
    # Center of the primary Gaussian functions
    c = 5.0
    
    # Standard deviation for the reference Upper Membership Function (UMF).
    # A larger std dev creates a wider curve, forming the upper bound.
    sigma_upper = 2.0
    
    # Standard deviation for the reference Lower Membership Function (LMF).
    # A smaller std dev creates a narrower curve, forming the lower bound.
    sigma_lower = 1.0
    
    # Alpha (α) is a parameter in [0,1] that controls the amount of "blurring" or
    # the uncertainty in the vertical dimension (the "type-3-ness").
    alpha = 0.5
    
    # --- Input Variables ---
    # Example values for the primary and secondary inputs.
    
    # Primary input variable
    x = 4.0
    
    # Secondary input variable (membership grade)
    u = 0.8
    
    # --- Mathematical Formulation and Calculation ---

    print("The mathematical formulation for the upper bound of a vertical cross-section of an IT3 MF, μ_upper(x, u), is:")
    print("μ_upper(x, u) = exp(-0.5 * ((u - μ_bar_A(x)) / σ_vertical(x))^2)")
    print("\nThis formulation depends on a reference Interval Type-2 (IT2) set, A, defined by Gaussian functions:")
    print(f"  - Upper Bound (UMF): μ_bar_A(x) = exp(-0.5 * ((x - c) / σ_upper)^2)")
    print(f"  - Lower Bound (LMF): μ_underline_A(x) = exp(-0.5 * ((x - c) / σ_lower)^2)")
    print("\nThe standard deviation of the vertical slice, σ_vertical(x), is proportional to the uncertainty in the reference set:")
    print("  - σ_vertical(x) = α * (μ_bar_A(x) - μ_underline_A(x))")

    print("\n--- Calculation with Specific Values ---")
    print(f"Parameters: c = {c}, σ_upper = {sigma_upper}, σ_lower = {sigma_lower}, α = {alpha}")
    print(f"Inputs: x = {x}, u = {u}")

    # Step 1: Calculate the bounds of the reference IT2 MF at x.
    # UMF of the reference IT2 MF at x
    mu_bar_A_x_exponent = -0.5 * ((x - c) / sigma_upper)**2
    mu_bar_A_x = math.exp(mu_bar_A_x_exponent)
    
    # LMF of the reference IT2 MF at x
    mu_underline_A_x_exponent = -0.5 * ((x - c) / sigma_lower)**2
    mu_underline_A_x = math.exp(mu_underline_A_x_exponent)

    print("\nStep 1: Calculate reference IT2 MF bounds at x =", x)
    print(f"μ_bar_A({x}) = exp(-0.5 * (({x} - {c}) / {sigma_upper})^2) = exp({mu_bar_A_x_exponent:.4f}) = {mu_bar_A_x:.4f}")
    print(f"μ_underline_A({x}) = exp(-0.5 * (({x} - {c}) / {sigma_lower})^2) = exp({mu_underline_A_x_exponent:.4f}) = {mu_underline_A_x:.4f}")

    # Step 2: Calculate the standard deviation for the vertical slice.
    fou_height = mu_bar_A_x - mu_underline_A_x
    sigma_vertical_x = 0
    if fou_height > 1e-9:  # Check for non-zero height to avoid division by zero
        sigma_vertical_x = alpha * fou_height
    
    print("\nStep 2: Calculate the standard deviation of the vertical slice at x =", x)
    print(f"σ_vertical({x}) = {alpha} * ({mu_bar_A_x:.4f} - {mu_underline_A_x:.4f}) = {sigma_vertical_x:.4f}")

    # Step 3: Calculate the final upper membership grade of the vertical slice.
    mu_bar_IT3_xu = 0.0
    final_exponent = float('-inf')
    if sigma_vertical_x > 1e-9:
        final_exponent = -0.5 * ((u - mu_bar_A_x) / sigma_vertical_x)**2
        mu_bar_IT3_xu = math.exp(final_exponent)

    print("\nStep 3: Calculate the final membership grade μ_upper(x, u)")
    print(f"μ_upper({x}, {u}) = exp(-0.5 * (({u} - {mu_bar_A_x:.4f}) / {sigma_vertical_x:.4f})^2)")
    # The following line shows the equation with all numbers substituted before the final calculation.
    print(f"             = exp({final_exponent:.4f})")
    print(f"             = {mu_bar_IT3_xu:.4f}")

    # Final answer for extraction
    global final_answer
    final_answer = mu_bar_IT3_xu

# This global variable will hold the final numerical result.
final_answer = 0.0
calculate_it3_vertical_slice_membership()
# The final numerical result required by the problem format.
# The explanation and step-by-step breakdown are printed above.
print(f"\n<<<>>>\n{final_answer}")
