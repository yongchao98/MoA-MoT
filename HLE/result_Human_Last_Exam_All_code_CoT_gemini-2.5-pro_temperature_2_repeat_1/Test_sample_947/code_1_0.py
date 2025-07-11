import numpy as np

def calculate_magnetic_field_z(x, z, Ha, Jc, d, w, D):
    """
    Calculates the z-component of the magnetic field for a stack of
    superconducting strips in the far-field limit |x| >> w.
    The formula is valid for an applied field Ha > H0, where H0 = Jc * d / pi.

    Args:
        x (float): x-coordinate in meters.
        z (float): z-coordinate in meters.
        Ha (float): Applied magnetic field in A/m.
        Jc (float): Critical current density in A/m^2.
        d (float): Strip thickness in meters.
        w (float): Strip half-width in meters.
        D (float): Spacing between strips in meters.

    Returns:
        float: The z-component of the magnetic field, Hz, in A/m.
    """
    # Check the validity of the approximations
    H0 = Jc * d / np.pi
    if Ha <= H0:
        print(f"Warning: The formula is derived for Ha > H0. (Current Ha={Ha:.2f}, H0={H0:.2f})")
    if abs(x) < 5 * w:  # A practical check for the far-field approximation |x| >> w
        print(f"Warning: The formula is an approximation for |x| >> w. (Current |x|={abs(x):.2e}, w={w:.2e})")

    # Term representing the coefficient of the summation part
    # Jc*d*w^2 is the magnitude of the magnetic moment m_z
    prefactor = (Jc * d * w**2 * np.pi) / (2 * D**2)

    # Complex argument for the cosecant function
    zeta = np.pi * (z + 1j * x) / D

    # The summation term S is the real part of csc^2(zeta)
    sin_zeta = np.sin(zeta)
    csc_sq_zeta = 1 / (sin_zeta**2)
    sum_term = np.real(csc_sq_zeta)

    # The total field is H_z = H_a + H_induced
    # H_induced = m_z * pi / (2D^2) * S = - (Jc*d*w^2) * pi / (2D^2) * S
    Hz = Ha - prefactor * sum_term

    # Output the equation with the calculated numbers
    print("The final equation for H_z is:")
    print("H_z = Ha - (Jc * d * w^2 * pi / (2 * D^2)) * Re[csc^2(pi * (z + i*x) / D)]\n")
    print("Substituting the numerical values:")
    print(f"H_z = {Ha} - {prefactor:.6f} * {sum_term:.4f}")
    final_result_str = f"H_z = {Hz:.4f} A/m"
    print(final_result_str)
    
    return Hz

if __name__ == '__main__':
    # --- User-definable parameters in SI units ---
    # Superconductor properties
    Jc = 2e10      # Critical current density in A/m^2
    d = 0.5e-6     # Strip thickness in m
    w = 150e-6     # Strip half-width in m

    # Stack geometry
    D = 5e-3       # Spacing between strips in m

    # External conditions
    Ha = 4000      # Applied field in A/m (must be > H0)

    # Point of interest
    x = 2e-3       # x-coordinate in m
    z = 2.5e-3     # z-coordinate in m (midway between two strips)
    # --- End of parameters ---

    # Perform the calculation and print the results
    final_value = calculate_magnetic_field_z(x, z, Ha, Jc, d, w, D)
    
    # Final answer in the required format
    print(f"\n<<<{final_value:.4f}>>>")
