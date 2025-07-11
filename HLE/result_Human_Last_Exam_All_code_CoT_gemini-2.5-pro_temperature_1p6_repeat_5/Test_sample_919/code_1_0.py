import numpy as np

def calculate_force_per_area(K0, a, d, y, mu, mu0):
    """
    Calculates the force per unit y-z area on the x=d interface.

    Args:
        K0 (float): Current sheet amplitude (in A/m).
        a (float): Spatial constant (in rad/m).
        d (float): Air gap thickness (in m).
        y (float): Position on the y-axis (in m).
        mu (float): Permeability of the magnetic material (in H/m).
        mu0 (float): Permeability of free space (in H/m).

    Returns:
        tuple: A tuple (fx, fy, fz) representing the force vector per area (in N/m^2).
    """
    # Calculate terms in the formula
    ay = a * y
    ad = a * d

    sin_ay_sq = np.sin(ay)**2
    cosh_ad = np.cosh(ad)
    sinh_ad = np.sinh(ad)
    mu_ratio = mu0 / mu

    # Calculate numerator and denominator of the main expression
    numerator__ = K0**2 * sin_ay_sq
    denominator_term = cosh_ad + mu_ratio * sinh_ad
    denominator_ = denominator_term**2

    # Calculate the magnitude of the force in the x-direction
    fx = (mu0 / 2.0) * (numerator_ / denominator_)

    # The force has only an x-component
    return (fx, 0, 0), {
        "sin(ay)^2": sin_ay_sq,
        "cosh(ad)": cosh_ad,
        "sinh(ad)": sinh_ad,
        "μ₀/μ": mu_ratio,
        "Numerator term (K₀²sin²(ay))": numerator_,
        "Denominator term [cosh(ad)+(μ₀/μ)sinh(ad)]^2": denominator_,
        "Total pre-factor (μ₀/2)": mu0 / 2.0
    }

def main():
    """
    Main function to define parameters, calculate the force, and print the results.
    """
    # --- Define physical parameters ---
    # Permeability of free space
    mu0 = 4 * np.pi * 1e-7  # H/m

    # Amplitude of the current sheet
    K0 = 100.0  # A/m
    # Spatial variation constant
    a = 50.0   # rad/m
    # Air gap thickness
    d = 0.01   # m (1 cm)
    # Permeability of the magnetic material (e.g., relative permeability of 1000)
    mu_relative = 1000.0
    mu = mu_relative * mu0  # H/m
    # Evaluate at a point of maximum field, where sin(ay) = 1
    y = np.pi / (2 * a)

    print("The derived formula for the force per unit area on the conductor is:")
    print("f/A = (μ₀/2) * [K₀² * sin²(ay)] / [cosh(ad) + (μ₀/μ) * sinh(ad)]² * î_x\n")

    print("--- Using the following physical values ---")
    print(f"K₀ (Current Amplitude) = {K0:.1f} A/m")
    print(f"a (Spatial Constant)    = {a:.1f} rad/m")
    print(f"d (Gap Thickness)       = {d:.3f} m")
    print(f"y (Position)            = {y:.4f} m")
    print(f"μ (Material Perm.)      = {mu:.6f} H/m (relative: {mu_relative})")
    print(f"μ₀ (Free Space Perm.)   = {mu0:.4e} H/m\n")

    # --- Calculation ---
    (force_vector, intermediate_values) = calculate_force_per_area(K0, a, d, y, mu, mu0)

    print("--- Calculating each number in the final equation ---")
    for name, value in intermediate_values.items():
        print(f"Value of {name:<42}: {value:.4e}")

    print("\n--- Final Result ---")
    print(f"The force per unit area vector is: ({force_vector[0]:.4e}, {force_vector[1]}, {force_vector[2]}) N/m^2")
    print(f"So, f/A = {force_vector[0]:.4e} î_x N/m^2")

if __name__ == "__main__":
    main()