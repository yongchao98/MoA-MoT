import numpy as np

def calculate_chi(k, h):
    """
    Calculates the stable amplitude factor chi.

    Args:
        k (float): The wavenumber (2*pi/lambda).
        h (float): The thickness of the tectonic plate.

    Returns:
        float: The non-dimensional amplitude factor chi.
    """
    # The amplitude factor depends on the non-dimensional parameter kh
    kh = k * h
    
    # The formula derived from solving the Stokes equations with the given boundary conditions
    # chi = -(kh * cosh(kh)) / ( (kh)^2 + (cosh(kh))^2 )
    numerator = kh * np.cosh(kh)
    denominator = kh**2 + np.cosh(kh)**2
    chi = -numerator / denominator
    
    return chi, kh

if __name__ == '__main__':
    # --- Example values for a tectonic plate ---
    # Plate thickness in meters
    h_plate = 100e3  # 100 km
    # Wavelength of the basal shear stress in meters
    lambda_stress = 200e3  # 200 km
    
    # Calculate the wavenumber k
    wavenumber_k = 2 * np.pi / lambda_stress
    
    # Calculate the amplitude factor chi and the intermediate parameter kh
    amplitude_factor_chi, kh_value = calculate_chi(wavenumber_k, h_plate)

    # Output the parameters and the final result
    print("--- Input Parameters ---")
    print(f"Plate thickness (h): {h_plate / 1000} km")
    print(f"Wavelength (lambda): {lambda_stress / 1000} km")
    
    print("\n--- Calculated Values ---")
    print(f"Wavenumber (k = 2*pi/lambda): {wavenumber_k:.6f} 1/m")
    print(f"Non-dimensional parameter (kh): {kh_value:.6f}")

    print("\n--- Final Equation ---")
    print(f"The equation for the amplitude factor is: chi = - (kh * cosh(kh)) / (kh^2 + cosh(kh)^2)")
    print(f"Plugging in the numbers:")
    print(f"chi = - ({kh_value:.4f} * {np.cosh(kh_value):.4f}) / (({kh_value:.4f})^2 + ({np.cosh(kh_value):.4f})^2)")
    print(f"chi = - ({kh_value * np.cosh(kh_value):.4f}) / ({kh_value**2:.4f} + {np.cosh(kh_value)**2:.4f})")
    print(f"chi = - ({kh_value * np.cosh(kh_value):.4f}) / ({kh_value**2 + np.cosh(kh_value)**2:.4f})")
    
    print("\n--- Result ---")
    print(f"The stable amplitude factor chi is: {amplitude_factor_chi:.6f}")