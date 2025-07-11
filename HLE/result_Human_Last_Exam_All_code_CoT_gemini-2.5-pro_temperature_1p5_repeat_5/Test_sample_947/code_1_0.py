import numpy as np

def calculate_magnetic_field(Ha, H0, w, D, x, z):
    """
    Calculates the z-component of the magnetic field for a stack of
    superconducting strips in an applied field Ha.

    This formula is valid for Ha > H0 and |x| >> w.

    Args:
        Ha (float): Applied magnetic field.
        H0 (float): Full penetration field (H0 = Jc * d / pi).
        w (float): Half-width of the strips.
        D (float): Stacking interval between strips.
        x (float): x-coordinate for field calculation.
        z (float): z-coordinate for field calculation.

    Returns:
        float: The z-component of the total magnetic field H_z(x, z).
    """
    # The final expression for the total magnetic field H_z is:
    # H_z(x, z) = Ha + H_induced_z(x, z)
    # H_induced_z(x, z) = - (pi^2 * H0 * w^2 / D^2) * (1 - cosh(2*pi*x/D)*cos(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2

    print("--- Calculating components of the equation ---")

    # Prefactor term from the derivation
    prefactor = (np.pi**2 * H0 * w**2) / D**2
    print(f"Prefactor (pi^2 * H0 * w^2 / D^2): {prefactor:.4f}")

    # Arguments for the hyperbolic and trigonometric functions
    arg_x = 2 * np.pi * x / D
    arg_z = 2 * np.pi * z / D
    print(f"Argument for cosh (2*pi*x/D): {arg_x:.4f}")
    print(f"Argument for cos (2*pi*z/D): {arg_z:.4f}")

    # Evaluate the cosh and cos terms
    cosh_term = np.cosh(arg_x)
    cos_term = np.cos(arg_z)
    print(f"cosh(2*pi*x/D): {cosh_term:.4f}")
    print(f"cos(2*pi*z/D): {cos_term:.4f}")
    
    # Numerator of the fraction
    numerator = 1 - cosh_term * cos_term
    print(f"Numerator (1 - cosh*cos): {numerator:.4f}")

    # Denominator of the fraction
    denominator = (cosh_term - cos_term)**2
    if np.isclose(denominator, 0):
        # This occurs at the strip locations (z=nD) in the dipole model
        # and represents a singularity.
        return np.inf if numerator > 0 else -np.inf

    print(f"Denominator ((cosh - cos)^2): {denominator:.4f}")

    # Calculate the field from the induced currents in the stack
    H_induced_z = -prefactor * numerator / denominator
    print(f"Induced Field (H_induced_z): {H_induced_z:.4f}")

    # The total field is the sum of the applied and induced fields
    H_z = Ha + H_induced_z

    return H_z

if __name__ == '__main__':
    # Define example parameters consistent with the problem's assumptions
    # e.g., |x| >> w and Ha > H0
    w_val = 0.1   # Half-width of the strip (e.g., in mm)
    D_val = 1.0   # Stacking distance (e.g., in mm)
    x_val = 2.0   # x-position, satisfying |x| >> w
    z_val = 0.25  # z-position, between strips (z/D = 0.25)
    H0_val = 100.0 # Full penetration field (e.g., in Oe)
    Ha_val = 150.0 # Applied field, satisfying Ha > H0

    print("--- Input Parameters ---")
    print(f"Applied Field (Ha): {Ha_val}")
    print(f"Penetration Field (H0): {H0_val}")
    print(f"Strip Half-width (w): {w_val}")
    print(f"Stacking Distance (D): {D_val}")
    print(f"Position (x, z): ({x_val}, {z_val})\n")

    # Calculate and print the final result
    total_field_z = calculate_magnetic_field(Ha_val, H0_val, w_val, D_val, x_val, z_val)

    print("\n--- Final Result ---")
    # The final result is printed with all its numerical components, as requested.
    # H_z(x,z) = Ha + H_induced
    print(f"The z-component of the magnetic field at ({x_val}, {z_val}) is:")
    print(f"H_z = {Ha_val} + ({total_field_z - Ha_val:.4f})")
    print(f"H_z(x,z) = {total_field_z:.4f}")
    print(f"\nFinal expression: H_z(x, z) = Ha - (pi^2*H0*w^2/D^2) * (1 - cosh(2*pi*x/D)*cos(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2")

<<<H_z(x,z) = Ha - (pi^2*H0*w^2/D^2) * (1 - cosh(2*pi*x/D)*cos(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2>>>