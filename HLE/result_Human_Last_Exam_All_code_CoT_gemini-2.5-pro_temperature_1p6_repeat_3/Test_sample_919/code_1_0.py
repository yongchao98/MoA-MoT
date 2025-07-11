import numpy as np

def calculate_force_on_conductor(K0, a, y, d, mu, mu0):
    """
    Calculates the force per unit area on the x=d interface due to the current sheet.

    Args:
        K0 (float): Current sheet amplitude (A/m).
        a (float): Spatial frequency of the current in y-direction (rad/m).
        y (float): Position on the y-axis (m).
        d (float): Thickness of the air gap (m).
        mu (float): Permeability of the magnetic material (H/m).
        mu0 (float): Permeability of free space (H/m).

    Returns:
        float: The magnitude of the force per unit area (N/m^2).
    """

    # The equation for the force per unit area is:
    # f = (μ₀/2) * (K₀² * sin²(ay)) / [cosh(ad) + (μ₀/μ)sinh(ad)]² in the x-direction.
    # We will calculate each part of this equation.

    ad = a * d
    ay = a * y

    # Calculate numerator terms
    num_term1_mu = mu0 / 2.0
    num_term2_K = K0**2
    num_term3_sin = np.sin(ay)**2
    numerator = num_term1_mu * num_term2_K * num_term3_sin

    # Calculate denominator terms
    den_term1_cosh = np.cosh(ad)
    den_term2_sinh = (mu0 / mu) * np.sinh(ad)
    denominator_base = den_term1_cosh + den_term2_sinh
    denominator = denominator_base**2

    # Calculate the force magnitude
    force_magnitude = numerator / denominator

    # Print the equation and the value of each term
    print("Equation for force magnitude F/A:")
    print("F/A = (μ₀/2) * (K₀²) * (sin²(ay)) / [cosh(ad) + (μ₀/μ)sinh(ad)]²\n")

    print("--- Input Parameters ---")
    print(f"K₀ (Current Amplitude) = {K0} A/m")
    print(f"a (Spatial Frequency) = {a} rad/m")
    print(f"y (Position) = {y} m")
    print(f"d (Gap Thickness) = {d} m")
    print(f"μ (Material Permeability) = {mu:.4e} H/m")
    print(f"μ₀ (Free Space Permeability) = {mu0:.4e} H/m")
    print("-" * 25)

    print("--- Calculated Terms ---")
    print(f"μ₀/2 = {num_term1_mu:.4e}")
    print(f"K₀² = {num_term2_K:.4e}")
    print(f"sin²(ay) = {num_term3_sin:.4f}")
    print(f"cosh(ad) = {den_term1_cosh:.4f}")
    print(f"(μ₀/μ)sinh(ad) = {den_term2_sinh:.4f}")
    print(f"[denominator base]² = {denominator:.4f}")
    print("-" * 25)

    print("--- Final Result ---")
    print(f"The magnitude of the force per unit area is: {force_magnitude:.4e} N/m^2")
    print("The direction of the force is along the positive x-axis (+î_x).")
    
    return force_magnitude

# --- Example Usage ---
# Define physical constants and parameters
mu0 = 4 * np.pi * 1e-7  # Permeability of free space in H/m
mu_r = 1000.0          # Relative permeability of the magnetic material
mu = mu_r * mu0

K0_val = 100.0  # A/m
a_val = 50.0    # rad/m
d_val = 0.01    # 1 cm
y_val = np.pi / (2 * a_val) # y-position where sin(ay) is maximum (1)

if __name__ == "__main__":
    calculate_force_on_conductor(K0=K0_val, a=a_val, y=y_val, d=d_val, mu=mu, mu0=mu0)