import numpy as np

def calculate_chi(k, h):
    """
    Calculates the stable amplitude factor chi based on the derived formula.

    The formula for chi is:
    chi = - (k*h * cosh(k*h)) / ((k*h)**2 + (cosh(k*h))**2)

    Args:
        k (float): The wavenumber, defined as 2*pi/lambda.
        h (float): The thickness of the tectonic plate.

    Returns:
        float: The value of the dimensionless amplitude factor chi.
    """
    # Dimensionless parameter
    kh = k * h

    # Numerator of the formula for chi
    # The coefficient is -1
    numerator_coeff = -1.0
    # The full numerator is -1 * kh * cosh(kh)
    numerator_value = numerator_coeff * kh * np.cosh(kh)

    # Denominator of the formula for chi
    # The formula is (kh)^2 + (cosh(kh))^2
    denominator_value = (kh)**2 + (np.cosh(kh))**2

    # The stable amplitude factor
    chi = numerator_value / denominator_value

    # --- Outputting the numbers in the final equation ---
    print("For the given k and h, the parameter k*h is calculated.")
    print(f"k*h = {kh:.4f}")
    
    print("\nThe final equation is chi = Numerator / Denominator")
    
    print("\nNumerator components:")
    print(f"Coefficient: {numerator_coeff}")
    print(f"kh * cosh(kh): {kh * np.cosh(kh):.4f}")
    print(f"Final Numerator Value: {numerator_value:.4f}")

    print("\nDenominator components:")
    print(f"(kh)^2: {(kh)**2:.4f}")
    print(f"(cosh(kh))^2: {(np.cosh(kh))**2:.4f}")
    print(f"Final Denominator Value: {denominator_value:.4f}")

    return chi

# --- Example Calculation ---
# Use plausible values for a tectonic setting.
# Let's assume a wavelength (lambda) of 1000 km and a plate thickness (h) of 80 km.
lambda_val = 1000e3  # meters
h_val = 80e3       # meters

# Calculate the wavenumber k
k_val = 2 * np.pi / lambda_val

print(f"Calculating chi for wavelength lambda = {lambda_val/1e3} km and thickness h = {h_val/1e3} km.")
print("-" * 60)

# Calculate the amplitude factor
chi_value = calculate_chi(k_val, h_val)

print("-" * 60)
print("The final result for the stable amplitude factor is:")
print(f"chi = {chi_value:.4f}")

# The final answer is the formula for chi. We extract the formulaic result.
H = 'kh'
chi_formula = f"-({H}*cosh({H}))/(({H})^2 + cosh({H})^2)"
print(f"\nThe symbolic formula for chi is: {chi_formula}")

# Return the symbolic answer as requested by the format
final_answer = f"- (kh * cosh(kh)) / ((kh)^2 + cosh(kh)^2)"
print(f"\n<<<{- (k_val*h_val * np.cosh(k_val*h_val)) / ((k_val*h_val)**2 + np.cosh(k_val*h_val)**2)}>>>")