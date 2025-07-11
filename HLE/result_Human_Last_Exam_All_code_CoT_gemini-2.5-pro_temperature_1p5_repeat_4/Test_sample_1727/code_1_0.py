import numpy as np

def solve_for_chi():
    """
    This function calculates the stable amplitude factor chi for the tectonic plate deformation.
    It first presents the symbolic formula derived from the fluid dynamics analysis.
    Then, it calculates a numerical value for chi using typical geophysical parameters
    for a tectonic plate, showing the intermediate steps of the calculation.
    """
    
    # Define the physical parameters for an example case
    # Plate thickness, h (meters)
    h = 50e3  # 50 km
    # Wavelength of the basal shear stress, lambda (meters)
    lambda_val = 500e3  # 500 km

    # Calculate the dimensionless parameter H = kh
    # Wavenumber, k (1/m)
    k = 2 * np.pi / lambda_val
    # Dimensionless thickness, H
    H = k * h

    def calculate_chi(H):
        """
        Calculates the dimensionless amplitude factor chi based on the derived formula.
        H is the dimensionless parameter kh.
        """
        # Numerator of the expression for chi
        numerator = (1 - H) * np.sinh(H)
        
        # Denominator of the expression for chi
        denominator = np.cosh(H)**2 - np.sinh(2 * H) - H**2
        
        # The factor chi
        chi = numerator / denominator
        return chi, numerator, denominator

    # Calculate the numerical value of chi for the example
    chi_value, num_val, den_val = calculate_chi(H)

    # Print the results
    print("The stable amplitude factor chi is given by the symbolic equation:")
    print("chi = [(1 - kh) * sinh(kh)] / [cosh(kh)^2 - sinh(2*kh) - (kh)^2]")
    print("\nFor a tectonic plate with thickness h = 50 km and a forcing wavelength lambda = 500 km:")

    print(f"\n1. The wavenumber k = 2*pi/lambda is: {k:.4e} 1/m")
    print(f"2. The dimensionless parameter H = kh is: {H:.4f}")

    print("\n3. We evaluate the components of the formula for H = kh:")
    print(f"   sinh(H) = sinh({H:.4f}) = {np.sinh(H):.4f}")
    print(f"   cosh(H) = cosh({H:.4f}) = {np.cosh(H):.4f}")
    print(f"   sinh(2*H) = sinh({2*H:.4f}) = {np.sinh(2*H):.4f}")
    print(f"   H^2 = ({H:.4f})^2 = {H**2:.4f}")

    print("\n4. We compute the numerator and denominator of the equation for chi:")
    print(f"   Numerator = (1 - {H:.4f}) * {np.sinh(H):.4f} = {num_val:.4f}")
    print(f"   Denominator = {np.cosh(H):.4f}^2 - {np.sinh(2*H):.4f} - {H**2:.4f} = {den_val:.4f}")

    print("\n5. Finally, the amplitude factor chi is the ratio of the numerator to the denominator:")
    print(f"   chi = {num_val:.4f} / ({den_val:.4f}) = {chi_value:.4f}")

solve_for_chi()