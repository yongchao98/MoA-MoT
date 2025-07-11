import numpy as np

def calculate_amplitude_factor(wavelength_m, thickness_m):
    """
    Calculates the stable amplitude factor chi based on the derived formula.

    Args:
        wavelength_m (float): The wavelength lambda of the harmonic stress in meters.
        thickness_m (float): The thickness h of the tectonic plate in meters.
    """
    # Calculate the wavenumber k
    k = 2 * np.pi / wavelength_m

    # Calculate the dimensionless parameter kh
    kh = k * thickness_m

    # Calculate hyperbolic functions
    cosh_kh = np.cosh(kh)
    sinh_kh = np.sinh(kh)

    # Calculate the numerator of the chi formula
    numerator = kh

    # Calculate the denominator of the chi formula
    denominator = (1 - kh + kh**2) * cosh_kh + kh * sinh_kh

    # Calculate chi
    chi = numerator / denominator

    # Output the parameters and the final equation
    print("--- Input Parameters ---")
    print(f"Wavelength (lambda): {wavelength_m / 1000:.0f} km")
    print(f"Plate Thickness (h): {thickness_m / 1000:.0f} km")
    print("\n--- Calculated Values ---")
    print(f"Wavenumber (k): {k:.4e} 1/m")
    print(f"Dimensionless parameter (kh): {kh:.4f}")
    print(f"cosh(kh): {cosh_kh:.4f}")
    print(f"sinh(kh): {sinh_kh:.4f}")
    print("\n--- Final Equation ---")
    # The format string is complex, so we build it step-by-step
    equation_str = f"chi = {kh:.4f} / ((1 - {kh:.4f} + {kh:.4f}**2) * {cosh_kh:.4f} + {kh:.4f} * {sinh_kh:.4f})"
    print(equation_str)
    equation_eval_str = f"chi = {numerator:.4f} / ({denominator:.4f})"
    print(equation_eval_str)

    print("\n--- Result ---")
    print(f"The stable amplitude factor chi is: {chi:.4f}")
    return chi

if __name__ == '__main__':
    # Example values for a tectonic plate
    # Wavelength of basal shear stress, e.g., 1000 km
    lambda_val = 1000e3
    # Thickness of the plate, e.g., 100 km
    h_val = 100e3

    final_chi = calculate_amplitude_factor(lambda_val, h_val)
    # The final answer is returned in the specified format below
    # print(f"\n<<<{final_chi:.4f}>>>")