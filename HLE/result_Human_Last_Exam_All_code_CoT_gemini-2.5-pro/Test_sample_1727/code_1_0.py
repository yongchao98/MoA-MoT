import numpy as np

def calculate_amplitude_factor(wavelength_m, thickness_m):
    """
    Calculates the stable amplitude factor chi for tectonic plate deformation.

    This function implements the solution to the Stokes flow problem for a viscous
    layer subjected to a basal shear stress, assuming a free-slip top surface.

    Args:
        wavelength_m (float): The wavelength (lambda) of the basal shear stress in meters.
        thickness_m (float): The thickness (h) of the viscous layer in meters.

    Returns:
        float: The dimensionless amplitude factor chi.
    """
    # 1. Calculate the wavenumber k from the wavelength lambda
    k = 2 * np.pi / wavelength_m

    # 2. Calculate the dimensionless parameter H = k*h
    H = k * thickness_m

    # 3. Calculate the required hyperbolic functions
    sinh_H = np.sinh(H)
    cosh_H = np.cosh(H)

    # 4. Calculate the numerator and denominator of the formula for chi
    # The formula is derived as: χ = - (H * sinh(H)) / (H + sinh(H) * cosh(H))
    numerator = -H * sinh_H
    denominator = H + sinh_H * cosh_H
    
    # An equivalent form of the denominator using sinh(2H) is H + 0.5 * np.sinh(2*H)

    # 5. Calculate chi
    chi = numerator / denominator

    # --- Output the step-by-step calculation ---
    print("Step 1: The derived formula for the stable amplitude factor χ is:")
    print("χ = - (k*h) * sinh(k*h) / (k*h + sinh(k*h) * cosh(k*h))")
    print("-" * 50)

    print("Step 2: Using the provided input values:")
    print(f"Wavelength λ = {wavelength_m / 1e3:.0f} km")
    print(f"Plate thickness h = {thickness_m / 1e3:.0f} km")
    print(f"Wavenumber k = 2*pi/λ = {k:.4g} m⁻¹")
    print("-" * 50)

    print("Step 3: Calculating the dimensionless parameter H and hyperbolic functions:")
    print(f"H = k*h = {H:.4f}")
    print(f"sinh(H) = sinh({H:.4f}) = {sinh_H:.4f}")
    print(f"cosh(H) = cosh({H:.4f}) = {cosh_H:.4f}")
    print("-" * 50)

    print("Step 4: Substituting these values into the formula:")
    print(f"χ = - ({H:.4f}) * ({sinh_H:.4f}) / (({H:.4f}) + ({sinh_H:.4f}) * ({cosh_H:.4f}))")
    print(f"χ = {numerator:.4f} / {denominator:.4f}")
    print("-" * 50)

    print("Step 5: The final calculated stable amplitude factor is:")
    print(f"χ = {chi:.4f}")
    
    return chi

if __name__ == '__main__':
    # Example values for a typical lithospheric plate
    lambda_val = 1000e3  # meters (1000 km)
    h_val = 100e3       # meters (100 km)
    
    # Calculate and print the result
    final_chi = calculate_amplitude_factor(lambda_val, h_val)
    # The final answer is submitted in the format below
    # print(f"\n<<<{final_chi:.4f}>>>")