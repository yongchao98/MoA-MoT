import numpy as np

def calculate_magnetic_field(x, z, Ha, w, D, Jc, d):
    """
    Calculates the magnetic field for a stack of superconducting strips.

    This function calculates the z-component of the magnetic field B_z at a point (x, z)
    for an infinite stack of thin superconducting strips. The calculation is based on the
    Bean critical state model for Ha > H0, where H0 = Jc*d/pi.

    Args:
        x (float): x-coordinate for field calculation (m).
        z (float): z-coordinate for field calculation (m).
        Ha (float): Applied magnetic field (A/m).
        w (float): Half-width of the strips (m).
        D (float): Spacing between the strips (m).
        Jc (float): Critical current density (A/m^2).
        d (float): Thickness of the strips (m).

    Returns:
        float: The calculated magnetic field B_z in Tesla (T), or None if inputs are invalid.
    """
    # Physical constants
    mu_0 = 4 * np.pi * 1e-7  # Permeability of free space (T*m/A)

    # Calculate the characteristic field H0
    H0 = Jc * d / np.pi

    # The model is specified for Ha > H0.
    if Ha <= H0:
        print(f"Error: This model is valid for Ha > H0.")
        print(f"Provided Ha = {Ha:.2f} A/m, but calculated H0 = {H0:.2f} A/m.")
        return None

    # Step 1: Calculate the flux front position 'a' using the isolated strip approximation.
    # a = w / cosh(Ha / H0)
    ha_over_h0 = Ha / H0
    # cosh overflows for large arguments, where a would be effectively zero.
    if ha_over_h0 > 700:
        a = 0.0
    else:
        a = w / np.cosh(ha_over_h0)

    # Step 2: Calculate the components of the final formula.
    # The total field B_z = B_applied + B_induced
    # B_applied = mu_0 * Ha
    # B_induced = (mu_0 * H0 / 4) * log( Numerator / Denominator )
    
    # Check for invalid positions that coincide with current sheet edges
    if np.isclose(abs(x), w) or np.isclose(abs(x), a):
        print(f"Error: The position x={x} is too close to a physical boundary (w={w}, a={a}), where the ideal model has a singularity.")
        return None

    cos_term = np.cos(2 * np.pi * z / D)

    # Numerator part of the log argument
    cosh_num1 = np.cosh(2 * np.pi * (x - a) / D)
    cosh_num2 = np.cosh(2 * np.pi * (x + a) / D)
    numerator = (cosh_num1 - cos_term) * (cosh_num2 - cos_term)

    # Denominator part of the log argument
    cosh_den1 = np.cosh(2 * np.pi * (x - w) / D)
    cosh_den2 = np.cosh(2 * np.pi * (x + w) / D)
    denominator = (cosh_den1 - cos_term) * (cosh_den2 - cos_term)

    if numerator <= 0 or denominator <= 0:
        print(f"Error: Logarithm argument is not positive (N={numerator:.2e}, D={denominator:.2e}). Check input values.")
        return None

    log_term = np.log(numerator / denominator)

    B_applied = mu_0 * Ha
    B_induced = (mu_0 * H0 / 4) * log_term
    B_total = B_applied + B_induced

    # --- Outputting the results and the equation numbers ---
    print("--- Calculation Breakdown ---")
    print("The final equation for the magnetic field is:")
    print("B_z(x,z) = μ₀*Ha + (μ₀*H₀/4) * log(N/D)\n")
    print("where:")
    print("N = [cosh(2π(x-a)/D) - cos(2πz/D)] * [cosh(2π(x+a)/D) - cos(2πz/D)]")
    print("D = [cosh(2π(x-w)/D) - cos(2πz/D)] * [cosh(2π(x+w)/D) - cos(2πz/D)]")
    print("a = w / cosh(Ha / H₀)")
    print("H₀ = Jc * d / π\n")
    
    print("--- Numerical Values for the Equation ---")
    print(f"μ₀ = {mu_0:.3e} T*m/A")
    print(f"Ha = {Ha:.3e} A/m")
    print(f"Jc = {Jc:.3e} A/m^2")
    print(f"d = {d:.3e} m")
    print(f"w = {w:.3e} m")
    print(f"D = {D:.3e} m")
    print(f"x = {x:.3e} m")
    print(f"z = {z:.3e} m\n")

    print("--- Intermediate Calculations ---")
    print(f"H₀ = {Jc:.3e} * {d:.3e} / π = {H0:.3e} A/m")
    print(f"Ha / H₀ = {Ha:.3e} / {H0:.3e} = {ha_over_h0:.3f}")
    print(f"a = {w:.3e} / cosh({ha_over_h0:.3f}) = {a:.3e} m")
    print(f"log(N/D) = log({numerator:.3e} / {denominator:.3e}) = {log_term:.3f}")
    print("\n--- Final Result ---")
    print(f"Applied Field Component (B_applied) = {B_applied:.5e} T")
    print(f"Induced Field Component (B_induced) = {B_induced:.5e} T")
    print(f"Total Magnetic Field (B_z) = {B_total:.5e} T")

    return B_total

if __name__ == '__main__':
    # --- Example Usage ---
    # Define parameters for a typical high-temperature superconductor stack
    Jc_val = 1.0e10   # Critical current density in A/m^2
    d_val = 1.0e-6    # Strip thickness in m
    w_val = 1.0e-3    # Strip half-width in m (2mm wide strip)
    D_val = 5.0e-3    # Spacing between strips in m (D > w)
    Ha_val = 5000.0   # Applied field in A/m
    
    # Point to calculate the field at (outside the strip, between two strips)
    x_val = 2.0 * w_val
    z_val = 0.5 * D_val
    
    print("Calculating magnetic field with example parameters:\n")
    final_B_field = calculate_magnetic_field(x=x_val, z=z_val, Ha=Ha_val, w=w_val, D=D_val, Jc=Jc_val, d=d_val)
    # The function prints all the details.
    
    # <<< B_z >>> This is just a placeholder for a single-value answer format which is not applicable here.
    # The full output is provided by the print statements within the function.
