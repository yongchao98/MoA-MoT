import math

def calculate_it3_vertical_slice():
    """
    Calculates and explains the vertical cross-section of a Gaussian IT3 MF.
    """
    # --- Model Parameters ---
    # Input values
    x = 6.0
    u = 0.878

    # Parameters for the primary IT2 MF bounds
    c = 5.0          # Center of the primary Gaussians
    h_UMF = 1.0      # Height of the primary Upper Membership Function (UMF)
    sigma_UMF = 2.0  # Std. dev. of the primary UMF (smaller sigma -> narrower curve)
    h_LMF = 0.9      # Height of the primary Lower Membership Function (LMF)
    sigma_LMF = 4.0  # Std. dev. of the primary LMF (larger sigma -> wider curve)

    # Parameters for the vertical slice (secondary IT2 MF)
    alpha_bar = 1.0      # Upper scaling factor for the vertical slice FOU
    alpha_underline = 0.5# Lower scaling factor for the vertical slice FOU
    k = 0.4              # Proportionality constant for the vertical slice std. dev.

    # --- Formulation and Calculation ---

    print("Mathematical Formulation for Vertical Cross-Sections of an IT3 MF")
    print("------------------------------------------------------------------")
    print("The vertical cross-section at a fixed 'x' is an Interval Type-2 Fuzzy Set over 'u'.")
    print("Its membership is an interval [μ_lower(x, u), μ_upper(x, u)], formulated as:")
    print("\n  μ_upper(x, u) = ᾱ * exp( -0.5 * ( (u - cᵤ(x)) / σᵤ(x) )² )")
    print("  μ_lower(x, u) = α̲ * exp( -0.5 * ( (u - cᵤ(x)) / σᵤ(x) )² )")
    print("\nWhere the parameters are derived from primary bounds μ̄(x) and μ̲(x):")
    print("  Primary UMF:       μ̄(x) = h_UMF * exp(-0.5 * ((x - c) / σ_UMF)²) ")
    print("  Primary LMF:       μ̲(x) = h_LMF * exp(-0.5 * ((x - c) / σ_LMF)²) ")
    print("  Vertical Center:   cᵤ(x) = (μ̄(x) + μ̲(x)) / 2")
    print("  Vertical Std. Dev.: σᵤ(x) = k * (μ̄(x) - μ̲(x))")
    
    print("\n--- Calculation with Specific Values ---")
    print(f"Given: x={x}, u={u}, c={c}, h_UMF={h_UMF}, σ_UMF={sigma_UMF}, h_LMF={h_LMF}, σ_LMF={sigma_LMF}, ᾱ={alpha_bar}, α̲={alpha_underline}, k={k}")

    # Step 1: Calculate primary bounds μ̄(x) and μ̲(x)
    mu_bar_x = h_UMF * math.exp(-0.5 * ((x - c) / sigma_UMF)**2)
    mu_underline_x = h_LMF * math.exp(-0.5 * ((x - c) / sigma_LMF)**2)
    print("\nStep 1: Calculate Primary Bounds μ̄(x) and μ̲(x)")
    print(f"μ̄({x}) = {h_UMF} * exp(-0.5 * (({x} - {c}) / {sigma_UMF})²) = {mu_bar_x:.4f}")
    print(f"μ̲({x}) = {h_LMF} * exp(-0.5 * (({x} - {c}) / {sigma_LMF})²) = {mu_underline_x:.4f}")

    # Initialize results
    mu_upper_x_u = 0.0
    mu_lower_x_u = 0.0

    # Step 2: Check domain and calculate vertical slice parameters
    print("\nStep 2: Calculate Vertical Slice Parameters cᵤ(x) and σᵤ(x)")
    if not (mu_underline_x <= u <= mu_bar_x):
        print(f"The value u = {u} is outside the valid domain [{mu_underline_x:.4f}, {mu_bar_x:.4f}]. Membership is [0, 0].")
    else:
        c_u_x = (mu_bar_x + mu_underline_x) / 2.0
        sigma_u_x_width = (mu_bar_x - mu_underline_x)
        sigma_u_x = k * sigma_u_x_width
        print(f"cᵤ({x}) = ({mu_bar_x:.4f} + {mu_underline_x:.4f}) / 2 = {c_u_x:.4f}")
        print(f"σᵤ({x}) = {k} * ({mu_bar_x:.4f} - {mu_underline_x:.4f}) = {sigma_u_x:.4f}")
        
        # Step 3: Calculate the final membership bounds
        print("\nStep 3: Calculate Final Membership Bounds with Substituted Values")
        # Handle case where primary bounds are equal, so sigma_u_x is 0
        if sigma_u_x < 1e-9:
            # Membership is non-zero only if u is exactly at the center
            if abs(u - c_u_x) < 1e-9:
                exp_term = 1.0
            else:
                exp_term = 0.0
        else:
            exp_term = math.exp(-0.5 * ((u - c_u_x) / sigma_u_x)**2)
        
        mu_upper_x_u = alpha_bar * exp_term
        mu_lower_x_u = alpha_underline * exp_term

        # Print the final equation with all numbers
        print(f"μ_upper({x}, {u}) = {alpha_bar} * exp(-0.5 * (({u} - {c_u_x:.4f}) / {sigma_u_x:.4f})²) = {mu_upper_x_u:.4f}")
        print(f"μ_lower({x}, {u}) = {alpha_underline} * exp(-0.5 * (({u} - {c_u_x:.4f}) / {sigma_u_x:.4f})²) = {mu_lower_x_u:.4f}")
    
    # Final Answer
    final_answer = f"[{mu_lower_x_u:.4f}, {mu_upper_x_u:.4f}]"
    print("\n<<<" + final_answer + ">>>")

if __name__ == '__main__':
    calculate_it3_vertical_slice()