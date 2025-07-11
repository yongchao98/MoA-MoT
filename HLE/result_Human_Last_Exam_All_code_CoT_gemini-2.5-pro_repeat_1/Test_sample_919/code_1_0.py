import numpy as np

def calculate_force_per_area(K0, a, y, d, mu, mu0):
    """
    Calculates the force per unit y-z area on the x=d interface.

    Args:
        K0 (float): Current sheet amplitude (A/m).
        a (float): Spatial frequency of the current (1/m).
        y (float): Position along the y-axis (m).
        d (float): Thickness of the air gap (m).
        mu (float): Permeability of the magnetic material (H/m).
        mu0 (float): Permeability of free space (H/m).

    Returns:
        float: The x-component of the force per unit area (N/m^2).
    """
    # Calculate intermediate terms for the equation
    # Denominator term D = [cosh(ad) + (mu0/mu) * sinh(ad)]
    ad = a * d
    cosh_ad = np.cosh(ad)
    sinh_ad = np.sinh(ad)
    mu_ratio = mu0 / mu
    denominator = cosh_ad + mu_ratio * sinh_ad

    # Numerator term N = - (mu0 / 2) * K0^2 * sin^2(ay)
    ay = a * y
    sin_ay_sq = np.sin(ay)**2
    numerator = - (mu0 / 2) * K0**2 * sin_ay_sq

    # Force per unit area f_x = N / D^2
    force_x = numerator / (denominator**2)

    # Print the equation with numerical values for clarity
    print("Calculating the force per unit area (x-component):")
    print(f"f_x = - (mu0 * K0^2 * sin(a*y)^2) / (2 * [cosh(a*d) + (mu0/mu) * sinh(a*d)]^2)")
    print("Substituting the given values:")
    print(f"Numerator = - ({mu0:.2e} / 2) * {K0}^2 * sin({a}*{y})^2 = {numerator:.4e}")
    print(f"Denominator^2 = [cosh({a}*{d}) + ({mu0:.2e}/{mu:.2e}) * sinh({a}*{d})]^2 = {denominator**2:.4e}")
    
    return force_x

if __name__ == '__main__':
    # Define physical constants and example parameters
    mu0 = 4 * np.pi * 1e-7  # Permeability of free space (H/m)
    
    # Example parameters for the system
    K0 = 100.0   # Current sheet amplitude in A/m
    a = 10.0     # Spatial variation constant in 1/m
    d = 0.01     # Air gap thickness in m (1 cm)
    y = 0.05     # Position on the y-axis in m (5 cm)
    # Relative permeability of the magnetic material (e.g., for a soft iron)
    mu_r = 4000.0
    mu = mu_r * mu0  # Permeability of the magnetic material in H/m

    # Calculate the force
    force_x_component = calculate_force_per_area(K0, a, y, d, mu, mu0)

    print("\nFinal Result:")
    print(f"The force per unit area on the conductor at x={d}m is {force_x_component:.4e} N/m^2 in the x-direction.")