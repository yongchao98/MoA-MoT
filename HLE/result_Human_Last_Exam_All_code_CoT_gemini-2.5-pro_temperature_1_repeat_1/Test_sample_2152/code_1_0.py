import math

def calculate_ratio():
    """
    Calculates the ratio of differential cross-sections for scattering off a 
    magnetic monopole vs. a magnetic dipole.
    """
    
    # The problem asks for the ratio R = (dσ/dΩ)_monopole / (dσ/dΩ)_dipole.
    # After derivation using the first Born approximation, the ratio simplifies to:
    # R = (3 * h_bar^2 * e_m^2) / (m^2 * v^2 * mu^2 * sin^2(theta))
    
    # Given parameters:
    # e_m = e / 16
    # v = c / 100
    # mu = 25 * mu_B, where mu_B = e*h_bar / (2*m*c) (in Gaussian units)
    # The particle is a charged particle with charge 'e', we assume it's an electron, so m = m_e.
    
    # Substituting these into the ratio formula:
    # R = (3 * h_bar^2 * (e/16)^2) / (m_e^2 * (c/100)^2 * (25 * e*h_bar/(2*m_e*c))^2 * sin^2(theta))
    # R = (3 * h_bar^2 * e^2 / 256) / (m_e^2 * c^2/10000 * 625 * e^2*h_bar^2/(4*m_e^2*c^2) * sin^2(theta))
    
    # All physical constants (e, h_bar, m_e, c) cancel out:
    # R = (3 / 256) / ( (1/10000) * (625/4) * sin^2(theta) )
    # R = (3 / 256) * (40000 / 625) / sin^2(theta)
    # R = (3 / 256) * 64 / sin^2(theta)
    # R = (3 * 64) / 256 / sin^2(theta)
    # R = 192 / 256 / sin^2(theta)
    # R = 3 / (4 * sin^2(theta))
    
    # The angle is given as theta = pi/30
    theta_rad = math.pi / 30
    
    # Perform the final calculation
    sin_theta = math.sin(theta_rad)
    sin_theta_sq = sin_theta**2
    denominator = 4 * sin_theta_sq
    ratio = 3 / denominator

    print("The final expression for the ratio R is:")
    print("R = 3 / (4 * sin^2(pi/30))")
    print("\nLet's plug in the numbers for the final equation:")
    print(f"The equation is: R = 3 / (4 * sin^2({theta_rad:.5f}))")
    print(f"First, calculate sin(pi/30):")
    print(f"sin({theta_rad:.5f}) = {sin_theta:.5f}")
    print(f"Next, square this value:")
    print(f"sin^2({theta_rad:.5f}) = {sin_theta_sq:.5f}")
    print(f"Then, multiply by 4:")
    print(f"4 * {sin_theta_sq:.5f} = {denominator:.5f}")
    print(f"Finally, divide 3 by this result:")
    print(f"R = 3 / {denominator:.5f}")
    print(f"\nThe final calculated ratio is: {ratio}")

calculate_ratio()