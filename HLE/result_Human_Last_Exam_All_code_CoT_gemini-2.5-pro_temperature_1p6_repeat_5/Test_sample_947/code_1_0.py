import numpy as np
from scipy.optimize import root_scalar

def solve_for_a(Ha, Jc, d, w, D):
    """
    Solves for the penetration depth 'a' for a given applied field 'Ha'.
    The equation is Ha = (Jc*d/pi) * ln(sinh(pi*w/D) / sinh(pi*a/D)).
    """
    H0 = Jc * d / np.pi
    if Ha <= 0:
        # No applied field, no penetration.
        return w
        
    # The field for full penetration (a=0) is infinite in this model.
    # We can calculate a for any Ha > 0.
    
    # Target value for the ratio of sinh functions
    target_ratio = np.exp(Ha / H0)

    # Function to find the root of. We are solving for 'a'.
    # f(a) = sinh(pi*w/D) / sinh(pi*a/D) - target_ratio = 0
    def equation_for_a(a):
        if a <= 0 or a >= w:
            return np.inf
        # Use log to avoid dealing with very large numbers from sinh
        # ln(sinh(pi*w/D)) - ln(sinh(pi*a/D)) - Ha/H0 = 0
        log_sinh_w = np.log(np.sinh(np.pi * w / D))
        log_sinh_a = np.log(np.sinh(np.pi * a / D))
        return log_sinh_w - log_sinh_a - (Ha / H0)

    # Solve for 'a' in the interval (0, w)
    try:
        sol = root_scalar(equation_for_a, bracket=[1e-9 * w, 0.99999 * w], method='brentq')
        if sol.converged:
            return sol.root
        else:
            raise RuntimeError("Failed to converge for 'a'")
    except ValueError:
        print("Error: The bracket for the root finding is not valid.")
        print(f"f(1e-9*w) = {equation_for_a(1e-9*w)}, f(0.999*w) = {equation_for_a(0.999*w)}")
        return None

def calculate_magnetic_field(x, z, Ha, Jc, d, w, D, a):
    """
    Calculates the magnetic field H_z at position (x, z).
    """
    if a is None:
        return None

    # Helper function for the log argument terms
    def term(val, x_coord, z_coord, D_coord):
        return np.cosh(2 * np.pi * (x_coord + val) / D_coord) - np.cos(2 * np.pi * z_coord / D_coord)

    # Numerator of the log argument
    num1 = term(a, x, z, D)
    num2 = term(-a, x, z, D)
    numerator = num1 * num2

    # Denominator of the log argument
    den1 = term(w, x, z, D)
    den2 = term(-w, x, z, D)
    denominator = den1 * den2

    if numerator <= 0 or denominator <= 0:
        print("Warning: Log argument is non-positive. This may occur at the boundaries.")
        log_term_val = -np.inf
    else:
        log_term_val = np.log(numerator / denominator)

    # Field from the screening currents
    H_induced = (Jc * d / (4 * np.pi)) * log_term_val

    # Total field
    H_total = Ha + H_induced
    
    print("--- Intermediate Calculation Steps ---")
    print(f"Applied Field (Ha): {Ha:.4f} A/m")
    print(f"Penetration depth (a): {a*1e6:.4f} um")
    print(f"Numerator of log argument: {numerator:.4e}")
    print(f"Denominator of log argument: {denominator:.4e}")
    print(f"Value of log term: {log_term_val:.4f}")
    print(f"Induced Field (H_induced): {H_induced:.4f} A/m")
    print("------------------------------------")
    
    return H_total


if __name__ == '__main__':
    # --- Parameters ---
    # Geometric parameters (in meters)
    w = 1.0e-3      # strip half-width
    d = 1.0e-6      # strip thickness
    D = 100.0e-6    # stacking interval

    # Superconducting parameters
    Jc = 3.0e10     # Critical current density (A/m^2)

    # Applied field (A/m)
    Ha = 15000.0

    # Position to calculate field (in meters)
    x = 0.5 * w
    z = 0.25 * D

    # --- Calculations ---
    # Characteristic field H0 = Jc*d/pi
    H0 = Jc * d / np.pi
    print(f"Strip half-width (w): {w*1e3} mm")
    print(f"Stacking interval (D): {D*1e6} um")
    print(f"Characteristic Field (H0 = Jc*d/pi): {H0:.2f} A/m")
    if Ha <= H0:
        print(f"\nWarning: The applied field Ha={Ha:.2f} A/m is not greater than H0={H0:.2f} A/m.")
        print("The formulas are still valid for any Ha > 0, but the problem specifies Ha > H0.")

    print("\n1. Solving for penetration depth 'a'...")
    penetration_depth_a = solve_for_a(Ha, Jc, d, w, D)
    if penetration_depth_a is None:
        print("Could not solve for 'a'. Exiting.")
    else:
        print(f"Solved 'a' = {penetration_depth_a*1e6:.4f} um (for w = {w*1e6} um)")

        # Validate |x| >> a condition
        if abs(x) > 3 * penetration_depth_a:
            print(f"Condition |x| >> a is met for x={x*1e6} um.")
        else:
            print(f"Warning: Condition |x| >> a may not be well-satisfied for x={x*1e6} um.")
        
        print(f"\n2. Calculating magnetic field at (x={x*1e6:.2f} um, z={z*1e6:.2f} um)...")
        Hz = calculate_magnetic_field(x, z, Ha, Jc, d, w, D, penetration_depth_a)
        if Hz is not None:
            print(f"\nFinal Result: The magnetic field H_z is {Hz:.4f} A/m")
