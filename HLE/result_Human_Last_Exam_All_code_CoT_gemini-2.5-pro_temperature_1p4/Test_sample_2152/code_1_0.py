import numpy as np

def solve_scattering_ratio():
    """
    Calculates the ratio of differential cross-sections for scattering by a magnetic monopole
    and a magnetic dipole.
    """
    # The formula for the ratio R, after substituting the physical quantities and simplifying, is:
    # R = (3 * g^2 * hbar^2) / (4 * mu^2 * p^2 * cos^2(theta/2))
    # where:
    # g = e_m = e/16 (monopole charge)
    # mu = 25 * mu_B = 25 * (e*hbar / (2*m_e*c)) (dipole moment in Gaussian units)
    # p = m*v, where v = c/100 (particle momentum)
    # We assume the particle is an electron, so m = m_e.

    # After substituting and cancelling terms (e, hbar, m_e, c), the formula simplifies to:
    # R = (3/256) / (625 / 10000 * cos^2(theta/2))
    # R = 3 / (16 * cos^2(theta/2))
    
    # Given scattering angle
    theta_rad = np.pi / 30

    # Numerator of the simplified ratio expression
    numerator = 3.0
    
    # Denominator components
    denominator_factor1 = 16.0
    cos_theta_half = np.cos(theta_rad / 2)
    denominator_factor2 = cos_theta_half**2

    # Calculate the final ratio
    ratio = numerator / (denominator_factor1 * denominator_factor2)

    print("The simplified formula for the ratio of differential cross-sections is:")
    print("R = 3 / (16 * cos^2(theta/2))")
    print("\nFor the given angle theta = pi/30:")
    
    print("\nThe final equation is:")
    print(f"{numerator} / ({denominator_factor1} * {denominator_factor2:.6f})")

    print(f"\nThe calculated ratio is: {ratio}")

solve_scattering_ratio()
<<<0.1879995163355203>>>