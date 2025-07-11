import math

def calculate_it3_vertical_slice(x_val):
    """
    Calculates and prints the UMF and LMF equations for the vertical
    slice of an IT3 MF at a specific primary input value x_val.
    """
    print(f"Determining the formulation for the vertical cross-section at primary input x = {x_val}\n")

    # --- Define Parameter Functions ---
    # These functions determine how the shape of the vertical slice (an IT2 FS)
    # changes as the primary input x changes.
    
    # c(x): Center of the Gaussian for the secondary variable u.
    # For this example, let's assume it's constant.
    center_func = lambda x: 0.5
    
    # sigma_L(x): Standard deviation for the UMF.
    # It's a function of x, causing the uncertainty to vary.
    sigma_lower_func = lambda x: 0.12 + 0.08 * math.exp(-((x - 5.0)**2) / 8.0)
    
    # sigma_U(x): Standard deviation for the LMF.
    # It must be >= sigma_L(x).
    sigma_upper_func = lambda x: 0.25 + 0.10 * math.exp(-((x - 5.0)**2) / 8.0)

    # --- Calculate Parameters for the chosen x_val ---
    c_val = center_func(x_val)
    sigma_L_val = sigma_lower_func(x_val)
    sigma_U_val = sigma_upper_func(x_val)

    # --- Denominators for the exponent ---
    denom_L = 2 * (sigma_L_val ** 2)
    denom_U = 2 * (sigma_U_val ** 2)

    print("--- Calculated Parameters ---")
    print(f"Center c({x_val}) = {c_val:.4f}")
    print(f"UMF Standard Deviation σ_L({x_val}) = {sigma_L_val:.4f}")
    print(f"LMF Standard Deviation σ_U({x_val}) = {sigma_U_val:.4f}\n")
    
    print("--- Mathematical Formulation of the Vertical Cross-Section ---")
    
    # --- Print the Upper Membership Function (UMF) Equation ---
    print("1. Upper Membership Function (UMF):")
    # Show the full equation with each number
    print(f"   overline_μ(u | {x_val}) = exp( -(u - {c_val:.4f})**2 / (2 * {sigma_L_val:.4f}**2) )")
    # Show the equation with the denominator calculated
    print(f"   overline_μ(u | {x_val}) = exp( -(u - {c_val:.4f})**2 / {denom_L:.4f} )\n")

    # --- Print the Lower Membership Function (LMF) Equation ---
    print("2. Lower Membership Function (LMF):")
    # Show the full equation with each number
    print(f"   underline_μ(u | {x_val}) = exp( -(u - {c_val:.4f})**2 / (2 * {sigma_U_val:.4f}**2) )")
    # Show the equation with the denominator calculated
    print(f"   underline_μ(u | {x_val}) = exp( -(u - {c_val:.4f})**2 / {denom_U:.4f} )")


if __name__ == '__main__':
    # Define a fixed value for the primary input variable x
    primary_input_x = 5.0
    
    calculate_it3_vertical_slice(primary_input_x)
