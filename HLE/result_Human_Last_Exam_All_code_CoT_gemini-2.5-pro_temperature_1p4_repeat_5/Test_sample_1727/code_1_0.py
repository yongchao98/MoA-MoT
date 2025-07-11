import numpy as np

def calculate_chi(k, h):
    """
    Calculates the stable amplitude factor chi for a given wavenumber k and plate thickness h.

    Args:
        k (float): The wavenumber of the basal shear stress (in 1/meters).
        h (float): The thickness of the tectonic plate (in meters).

    Returns:
        float: The non-dimensional stable amplitude factor chi.
    """
    H = k * h
    
    # Numerator of the expression for chi: 2*sinh(H) - H*cosh(H)
    numerator = 2 * np.sinh(H) - H * np.cosh(H)
    
    # Denominator of the expression for chi: H^2 - 2*sinh(H)^2
    # This is equivalent to H^2 + 1 - cosh(2*H)
    denominator = H**2 - 2 * np.sinh(H)**2
    
    if denominator == 0:
        return np.inf # Avoid division by zero
        
    chi = numerator / denominator
    return chi

def main():
    """
    Main function to demonstrate the calculation of chi.
    """
    # Example parameters
    # Plate thickness h = 50 km
    h_val = 50e3  # meters
    # Wavelength lambda = 400 km
    lambda_val = 400e3  # meters
    # Wavenumber k = 2*pi/lambda
    k_val = 2 * np.pi / lambda_val
    
    # Dimensionless parameter H = k*h
    H_val = k_val * h_val
    
    # Calculate chi using the derived formula
    chi_value = calculate_chi(k_val, h_val)
    
    # The final equation for the amplitude factor chi as a function of H = kh is:
    # chi = (2*sinh(H) - 1*H*cosh(H)) / (H**2 - 2*sinh(H)**2)
    
    print("This script calculates the dimensionless amplitude factor 'chi'.")
    print("\n--- Final Equation ---")
    print("The derived equation for chi is:")
    # The prompt requests to output each number in the final equation.
    # Here are the numbers: 2, 1, 2.
    print("chi = (2 * sinh(H) - 1 * H * cosh(H)) / (H**2 - 2 * (sinh(H))**2)")
    print("where H = k * h is the dimensionless wavenumber.")
    
    print("\n--- Example Calculation ---")
    print(f"Given plate thickness h = {h_val/1000} km")
    print(f"Given wavelength lambda = {lambda_val/1000} km")
    print(f"Calculated wavenumber k = {k_val:.5f} 1/m")
    print(f"Dimensionless parameter H = k*h = {H_val:.4f}")
    
    print("\n--- Result ---")
    print(f"The stable amplitude factor chi is: {chi_value:.4f}")

if __name__ == "__main__":
    main()