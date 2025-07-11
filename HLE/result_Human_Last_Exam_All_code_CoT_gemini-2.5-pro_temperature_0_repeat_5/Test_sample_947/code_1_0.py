import numpy as np

def calculate_magnetic_field(Jc, d, w, D, Ha, x, z):
    """
    Calculates the magnetic field H_z for a stack of superconducting strips.

    The calculation is based on the critical state model for an infinite stack
    of thin strips, under the approximation that the observation point |x|
    is much larger than the flux penetration front 'a'.

    Args:
        Jc (float): Critical current density.
        d (float): Thickness of the strips.
        w (float): Half-width of the strips.
        D (float): Stacking interval.
        Ha (float): Applied magnetic field.
        x (float): x-coordinate of the observation point.
        z (float): z-coordinate of the observation point.
    """
    # --- Check for full penetration ---
    # The full penetration field Hp is the field at which a = 0.
    Hp = (Jc * d / D) * np.log(np.cosh(np.pi * w / D))
    if Ha >= Hp:
        print(f"Warning: Applied field Ha ({Ha}) is greater than or equal to the full penetration field Hp ({Hp:.4f}).")
        print("The flux front 'a' is 0 (fully penetrated state).")
        sinh2_A = 0.0
    else:
        # Calculate sinh^2(pi*a/D) from Ha
        # cosh(pi*a/D) = cosh(pi*w/D) * exp(-Ha*D/(Jc*d))
        cosh2_A = (np.cosh(np.pi * w / D)**2) * np.exp(-2 * Ha * D / (Jc * d))
        sinh2_A = cosh2_A - 1

    # --- Define complex variables ---
    zeta = x + 1j * z
    Z = np.pi * zeta / D
    W = np.pi * w / D

    # --- Calculate the terms of the equation ---
    # The total field is H_z = H_a + H_z_induced
    # H_z_induced is approximated as Term2 + Term3

    # Term 1: Applied Field
    term1_Ha = Ha

    # Term 2: Field from a fully penetrated strip of width 2w
    # This term is (Jc*d / (2D)) * Re[ln(1 - sinh^2(W)/sinh^2(Z))]
    term2_coeff = Jc * d / (2 * D)
    term2_arg = 1 - (np.sinh(W)**2) / (np.sinh(Z)**2)
    term2_val = term2_coeff * np.real(np.log(term2_arg))

    # Term 3: Correction term due to the shielded region (dipole-like field)
    # This term is (Jc*d / (2D)) * sinh^2(A) * Re[1/sinh^2(Z)]
    term3_coeff = (Jc * d / (2 * D)) * sinh2_A
    term3_geo = np.real(1 / (np.sinh(Z)**2))
    term3_val = term3_coeff * term3_geo

    # Total Field
    Hz_total = term1_Ha + term2_val + term3_val

    # --- Print the equation with numerical components ---
    print("The magnetic field H_z is calculated as: H_z = H_a + H_induced_term2 + H_induced_term3")
    print(f"H_a = {term1_Ha}")
    print(f"H_induced_term2 = {term2_val:.6f}")
    print(f"H_induced_term3 = {term3_val:.6f}")
    print("-" * 30)
    print(f"Total magnetic field H_z at (x={x}, z={z}) is: {Hz_total:.6f}")


if __name__ == '__main__':
    # --- Input Parameters ---
    Jc = 1.0e10  # Critical current density in A/m^2
    d = 1.0e-6   # Strip thickness in m
    w = 1.0e-3   # Strip half-width in m
    D = 2.0e-3   # Stacking interval in m
    Ha = 3000.0  # Applied magnetic field in A/m
    x = 1.0e-2   # x-coordinate in m
    z = 0.0      # z-coordinate in m

    # --- Execute Calculation ---
    calculate_magnetic_field(Jc, d, w, D, Ha, x, z)
