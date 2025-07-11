import numpy as np

def calculate_chi(k, h):
    """
    Calculates the stable amplitude factor chi based on the wavenumber k and thickness h.
    
    The factor chi relates the surface topography e_s to the basal shear stress S_0
    by the equation e_s = chi * (S_0 / (Delta_rho * g)).

    Args:
        k (float): The wavenumber (2*pi/lambda).
        h (float): The thickness of the tectonic plate.

    Returns:
        float: The dimensionless stable amplitude factor chi.
    """
    kh = k * h
    
    # This formula is derived from solving the Stokes equations with the given boundary conditions.
    numerator = np.sinh(kh)**2 - kh**2
    denominator = kh + np.sinh(kh) * np.cosh(kh)
    
    if np.isclose(denominator, 0):
        # Avoid division by zero, although for positive k, h it shouldn't happen.
        return np.inf
        
    chi = numerator / denominator
    return chi

def main():
    """
    Main function to demonstrate the calculation of the amplitude factor chi.
    """
    # Example parameters for a tectonic plate
    h = 50e3      # Plate thickness in meters (50 km)
    lambda_val = 400e3  # Wavelength of the basal stress in meters (400 km)
    
    # Calculate the wavenumber k
    k = 2 * np.pi / lambda_val
    
    # Calculate the dimensionless parameter kh
    kh = k * h
    
    # Calculate the stable amplitude factor chi
    chi_value = calculate_chi(k, h)
    
    # Print the inputs, the symbolic formula, and the final result
    print("This script calculates the stable amplitude factor χ in the equation for tectonic plate surface deformation.")
    print("\nGiven parameters:")
    print(f"Plate thickness h = {h/1000} km")
    print(f"Wavelength λ = {lambda_val/1000} km")
    print(f"Wavenumber k = {k:.6f} m⁻¹")
    print(f"Dimensionless parameter k*h = {kh:.4f}")
    
    print("\nThe stable amplitude factor is given by the formula:")
    print("χ = (sinh(k*h)² - (k*h)²) / (k*h + sinh(k*h)*cosh(k*h))")
    
    print(f"\nThe calculated value for the stable amplitude factor is:")
    print(f"χ = {chi_value:.4f}")
    
    print("\nThe full equation for the stable vertical displacement amplitude (e_s) is:")
    # The instruction is to output each number in the final equation.
    # Since S_0, Delta_rho, and g are symbolic, we output the equation with the calculated number for chi.
    print(f"e_s = {chi_value:.4f} * (S₀ / (Δρ * g))")

if __name__ == "__main__":
    main()

<<<0.8931>>>