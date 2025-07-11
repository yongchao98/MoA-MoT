import numpy as np

def calculate_superconducting_field(Ha, H0, w, D, x, z):
    """
    Calculates the magnetic field for a stack of superconducting strips.

    Args:
        Ha (float): Applied magnetic field.
        H0 (float): Characteristic field, H0 = Jc*d/pi.
        w (float): Half-width of the strips.
        D (float): Stacking interval.
        x (float): x-coordinate for field calculation.
        z (float): z-coordinate for field calculation.

    Returns:
        float: The z-component of the magnetic field H_z(x, z).
    """

    # The problem is for Ha > H0, where flux has significantly penetrated.
    if Ha <= H0:
        print("Warning: This formula is derived for the case Ha > H0.")
    
    if D <= 0 or w <= 0:
        raise ValueError("Dimensions w and D must be positive.")

    # 1. Calculate the position of the flux front 'a' for the stack
    # This is derived from the condition that the field at the center (0,0) is zero.
    # Ha + H_ind(0,0) = 0 => Ha = H0 * ln[sinh(pi*w/D) / sinh(pi*a/D)]
    try:
        sinh_term = np.sinh(np.pi * w / D) * np.exp(-Ha / H0)
        # The argument to arcsinh must be non-negative
        if sinh_term < 0:
            raise ValueError("Invalid parameters lead to negative sinh_term.")
        a = (D / np.pi) * np.arcsinh(sinh_term)
    except (OverflowError, ValueError) as e:
        print(f"Error calculating 'a': {e}. Check input parameters.")
        return None

    # Check the condition |x| >> a
    # This is a qualitative condition, we'll just print a note.
    # In practice, the formula is valid for any |x| > a.
    if abs(x) <= a:
        print(f"Note: The coordinate x={x} is inside the shielded region (|x| <= a={a:.4f}).")
        print("In this model, the field here should be close to zero.")

    # 2. Calculate the induced field H_ind from the entire stack.
    # The current is modeled as K = -sgn(x)*pi*H0 for a < |x| < w.
    # Summing the ln terms from all strips leads to the following expression
    # involving sinh functions.
    
    # Define dimensionless arguments for the transcendental functions
    u_x_plus_w = np.pi * (x + w) / D
    u_x_minus_w = np.pi * (x - w) / D
    u_x_plus_a = np.pi * (x + a) / D
    u_x_minus_a = np.pi * (x - a) / D
    v_z = np.pi * z / D
    
    # Numerator of the log's argument
    # N = (sinh^2(pi(x+w)/D) + sin^2(pi*z/D)) * (sinh^2(pi(x-w)/D) + sin^2(pi*z/D))
    num_term1 = np.sinh(u_x_plus_w)**2 + np.sin(v_z)**2
    num_term2 = np.sinh(u_x_minus_w)**2 + np.sin(v_z)**2
    
    # Denominator of the log's argument
    # D = (sinh^2(pi(x+a)/D) + sin^2(pi*z/D)) * (sinh^2(pi(x-a)/D) + sin^2(pi*z/D))
    den_term1 = np.sinh(u_x_plus_a)**2 + np.sin(v_z)**2
    den_term2 = np.sinh(u_x_minus_a)**2 + np.sin(v_z)**2

    # Prevent division by zero or log of zero, which occurs at the current sheets
    if den_term1 <= 0 or den_term2 <= 0:
        print(f"Error: Field calculation at a singularity (x={x}, z={z}). Cannot compute.")
        return None
        
    log_argument = (num_term1 * num_term2) / (den_term1 * den_term2)
    
    # Calculate induced field component
    H_ind = -(H0 / 4) * np.log(log_argument)

    # 3. Total field is the sum of applied and induced fields
    H_z = Ha + H_ind

    # 4. Print the full expression and the result
    print("--- Magnetic Field Calculation ---")
    print(f"Parameters: Ha={Ha}, H0={H0}, w={w}, D={D}, x={x}, z={z}")
    print(f"Calculated flux front position, a = {a:.4f}")
    print("\nFinal Equation:")
    print(f"H_z(x, z) = Ha + H_induced")
    print(f"H_z({x}, {z}) = {Ha} + ({H_ind:.4f})")
    print(f"H_z({x}, {z}) = {H_z:.4f}")
    print("------------------------------------")

    return H_z

if __name__ == '__main__':
    # --- User-defined parameters ---
    # These are example values.
    Ha_applied = 2.5   # Applied field (arbitrary units, e.g., A/m)
    H0_char = 1.0      # Characteristic field H0 = Jc*d/pi (same units as Ha)
    strip_half_width = 5.0 # w (arbitrary units, e.g., mm)
    stack_distance = 10.0  # D (same units as w)
    
    # Point of interest (x, z)
    x_coord = 7.0
    z_coord = 0.0      # Field in the plane of a strip (for z=n*D)
                       # z=D/2 would be mid-plane between strips.

    # Calculate and print the result
    calculate_superconducting_field(Ha_applied, H0_char, strip_half_width, stack_distance, x_coord, z_coord)
    
    print("\nExample 2: Point between strips")
    calculate_superconducting_field(Ha_applied, H0_char, strip_half_width, stack_distance, x_coord, stack_distance/2)

<<<H_z(x, z) = Ha - (H0/4) * ln[((sinh^2(pi*(x+w)/D) + sin^2(pi*z/D))*(sinh^2(pi*(x-w)/D) + sin^2(pi*z/D))) / ((sinh^2(pi*(x+a)/D) + sin^2(pi*z/D))*(sinh^2(pi*(x-a)/D) + sin^2(pi*z/D)))], where a = (D/pi)*arcsinh[sinh(pi*w/D)*exp(-Ha/H0)]>>>