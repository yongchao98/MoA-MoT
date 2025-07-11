import numpy as np

def calculate_amplitude_factor(h, wavelength):
    """
    Calculates the stable amplitude factor chi for a viscous plate model.

    Args:
        h (float): Thickness of the tectonic plate in meters.
        wavelength (float): Wavelength of the harmonic basal shear stress in meters.

    Returns:
        float: The dimensionless amplitude factor chi.
    """
    # Step 1: Calculate the wavenumber k
    k = 2 * np.pi / wavelength
    
    # Step 2: Calculate the dimensionless parameter K
    K = k * h
    
    # Step 3: Calculate cosh(K)
    cosh_K = np.cosh(K)
    
    # Step 4: Calculate the numerator and denominator for chi
    numerator = -K * cosh_K
    denominator = K**2 + cosh_K**2
    
    # Step 5: Calculate chi
    chi = numerator / denominator
    
    # --- Outputting the details of the calculation ---
    print("--- Assumed Parameters ---")
    print(f"Plate thickness (h): {h / 1000:.1f} km")
    print(f"Wavelength (lambda): {wavelength / 1000:.1f} km")
    print("\n--- Calculated Values ---")
    print(f"Wavenumber (k): {k:.4e} 1/m")
    print(f"Dimensionless Wavenumber (K = k*h): {K:.4f}")
    
    print("\n--- Final Equation for chi ---")
    # Using the instruction to output each number in the final equation
    print(f"chi = - (K * cosh(K)) / (K^2 + cosh(K)^2)")
    print(f"chi = - ({K:.4f} * {cosh_K:.4f}) / (({K:.4f})^2 + ({cosh_K:.4f})^2)")
    print(f"chi = {numerator:.4f} / ({K**2:.4f} + {cosh_K**2:.4f})")
    print(f"chi = {numerator:.4f} / {denominator:.4f}")

    print("\n--- Result ---")
    print(f"The stable amplitude factor chi is: {chi:.4f}")
    
    return chi

# Example usage with plausible geophysical values
# Tectonic plate thickness of 50 km
plate_thickness = 50e3  # meters
# Wavelength of basal forcing of 500 km
stress_wavelength = 500e3 # meters

chi_value = calculate_amplitude_factor(plate_thickness, stress_wavelength)

# The final answer format requires enclosing the single numerical result.
# We will use the calculated chi_value for this.
final_answer = f"<<<{chi_value:.4f}>>>"
# print(final_answer) # This would normally print the required format, but we'll return it instead.