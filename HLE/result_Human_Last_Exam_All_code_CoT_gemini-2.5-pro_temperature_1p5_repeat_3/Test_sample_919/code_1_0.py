import numpy as np

def calculate_force_on_conductor(K0, a, d, y, mu_0, mu):
    """
    Calculates the force per unit y-z area on the conducting plane at x=d.

    Args:
        K0 (float): Current sheet amplitude (in A/m).
        a (float): Spatial frequency of the current (in rad/m).
        d (float): Thickness of the air gap (in m).
        y (float): Position on the y-axis (in m).
        mu_0 (float): Permeability of free space.
        mu (float): Permeability of the magnetic material.
        
    Returns:
        float: The magnitude of the force per unit area in the x-direction (in N/m^2).
    """
    # Arguments for the hyperbolic and trig functions
    ad = a * d
    ay = a * y

    # Numerator of the force equation
    numerator = mu_0 * K0**2 * (np.sin(ay))**2

    # Denominator of the force equation
    cosh_ad = np.cosh(ad)
    sinh_ad = np.sinh(ad)
    denominator_base = cosh_ad + (mu_0 / mu) * sinh_ad
    denominator = 2 * denominator_base**2
    
    force_magnitude = numerator / denominator
    
    # Print the equation with values
    print("Equation for force per unit area (magnitude):")
    print(f"f/A = (μ₀ * K₀² * sin²(a*y)) / (2 * [cosh(a*d) + (μ₀/μ) * sinh(a*d)]²)")
    print("\nSubstituting the given values:")
    print(f"f/A = ({mu_0:.2e} * {K0}² * sin²({a:.2f}*{y:.3f})) / (2 * [cosh({a:.2f}*{d:.3f}) + ({mu_0:.2e}/{mu:.2e}) * sinh({a:.2f}*{d:.3f})]²)")
    
    print(f"\nCalculated terms:")
    print(f"  a*y = {ay:.4f} rad")
    print(f"  sin(ay) = {np.sin(ay):.4f}")
    print(f"  sin²(ay) = {np.sin(ay)**2:.4f}")
    print(f"  a*d = {ad:.4f} rad")
    print(f"  cosh(ad) = {cosh_ad:.4f}")
    print(f"  sinh(ad) = {sinh_ad:.4f}")
    print(f"  μ₀/μ = {mu_0/mu:.4f}")
    print(f"  Denominator Base = {denominator_base:.4f}")

    return force_magnitude

# --- Example Usage ---
# Physical constants and problem parameters
mu_0 = 4 * np.pi * 1e-7  # Permeability of free space (H/m)
mu_r = 150.0             # Relative permeability of the magnetic material
mu = mu_r * mu_0         # Permeability of the magnetic material

K0 = 100.0               # Current sheet amplitude (A/m)
wavelength = 0.1         # Wavelength of the current variation (m)
a = 2 * np.pi / wavelength # Spatial frequency (rad/m)
d = 0.01                 # Air gap thickness (m)
y = 0.0125               # Position on y-axis for evaluation (m)

# Calculate and print the result
force_per_area = calculate_force_on_conductor(K0, a, d, y, mu_0, mu)
print(f"\nFinal Result:")
print(f"The magnitude of the force per unit area at y = {y} m is: {force_per_area:.4f} N/m^2")
