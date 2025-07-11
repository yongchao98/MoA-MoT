import math

def calculate_magnetic_field(x, z, w, D, Jc, d, Ha):
    """
    Calculates the magnetic field H_z for a stack of superconducting strips.

    Args:
        x (float): x-coordinate of the observation point.
        z (float): z-coordinate of the observation point.
        w (float): Half-width of the strips.
        D (float): Stacking period (distance between strips).
        Jc (float): Critical current density.
        d (float): Thickness of the strips.
        Ha (float): Applied magnetic field in the z-direction.
    """
    # Step 1: Calculate the characteristic field H0
    try:
        H0 = Jc * d / math.pi
    except ZeroDivisionError:
        print("Error: math.pi is zero? This should not happen.")
        return

    # Check if the applied field is in the specified regime
    if Ha <= H0:
        print(f"Calculation valid for Ha > H0. Current values: Ha={Ha:.4g}, H0={H0:.4g}")
        print("The provided formula is not valid in this regime (Ha <= H0).")
        return

    # Step 2: Calculate the flux front position 'a'
    try:
        a = w / math.cosh(Ha / H0)
    except OverflowError:
        print("Error: Overflow in cosh calculation. Ha/H0 might be too large.")
        return

    # The condition |x| >> a must be met. We'll proceed assuming it holds.
    # Step 3: Calculate the self-field from the stack
    
    # Helper function for the terms in the formula
    def field_term(val, z_pos, period):
        # This calculates cosh(2*pi*val/D) - cos(2*pi*z/D)
        try:
            cosh_arg = 2 * math.pi * val / period
            cos_arg = 2 * math.pi * z_pos / period
            return math.cosh(cosh_arg) - math.cos(cos_arg)
        except OverflowError:
            # For large arguments, cosh(x) is approx 0.5*exp(x)
            # This helps avoid overflow while maintaining precision
            return 0.5 * math.exp(cosh_arg) - math.cos(cos_arg)

    # Numerator of the log argument
    # term for (x-a) and (x+a)
    num1 = field_term(x - a, z, D)
    num2 = field_term(x + a, z, D)
    numerator = num1 * num2

    # Denominator of the log argument
    # term for (x-w) and (x+w)
    den1 = field_term(x - w, z, D)
    den2 = field_term(x + w, z, D)
    denominator = den1 * den2

    if numerator <= 0 or denominator <= 0:
        print("Error: Argument of log is non-positive. Invalid parameters or coordinates.")
        return

    # Self-field H_z,self
    Hz_self = (H0 / 4) * math.log(numerator / denominator)

    # Step 4: Calculate the total field
    Hz_total = Ha + Hz_self
    
    # Step 5: Print the results and the equation
    print("--- Calculation Parameters ---")
    print(f"Coordinates: (x, z) = ({x}, {z})")
    print(f"Strip half-width, w = {w}")
    print(f"Stacking distance, D = {D}")
    print(f"Critical current density, Jc = {Jc:.2e}")
    print(f"Strip thickness, d = {d}")
    print(f"Applied field, Ha = {Ha:.4g}")
    print("\n--- Intermediate Calculations ---")
    print(f"Characteristic field, H0 = Jc*d/pi = {H0:.4g}")
    print(f"Flux front position, a = w/cosh(Ha/H0) = {a:.4g}")
    
    print("\n--- Final Equation with Values ---")
    # Using f-strings to format the equation string. {:.4g} formats numbers to 4 significant figures.
    equation = (
        f"H_z(x, z) = Ha + (H0 / 4) * ln[ (cosh(2*pi*(x-a)/D) - cos(2*pi*z/D)) * (cosh(2*pi*(x+a)/D) - cos(2*pi*z/D)) / ( (cosh(2*pi*(x-w)/D) - cos(2*pi*z/D)) * (cosh(2*pi*(x+w)/D) - cos(2*pi*z/D)) ) ]\n\n"
        f"H_z({x}, {z}) = {Ha:.4g} + ({H0:.4g} / 4) * ln[ "
        f"(cosh(2*pi*({x:.4g}-{a:.4g})/{D}) - cos(2*pi*{z}/{D})) * "
        f"(cosh(2*pi*({x:.4g}+{a:.4g})/{D}) - cos(2*pi*{z}/{D})) / ( "
        f"(cosh(2*pi*({x:.4g}-{w:.4g})/{D}) - cos(2*pi*{z}/{D})) * "
        f"(cosh(2*pi*({x:.4g}+{w:.4g})/{D}) - cos(2*pi*{z}/{D})) ) ]"
    )
    print(equation)
    
    print("\n--- Final Result ---")
    print(f"Self-field, Hz_self = {Hz_self:.5g}")
    print(f"Total magnetic field, Hz_total = {Hz_total:.5g}")
    
    # Return the final numeric answer as requested by the format.
    print(f"\n<<<{Hz_total:.5g}>>>")


if __name__ == '__main__':
    # --- User-defined parameters ---
    # You can change these values to explore different scenarios.
    
    # Coordinates of the point where the field is calculated
    x_coord = 1.1   # m
    z_coord = 0.1   # m
    
    # Physical properties of the superconducting strips
    strip_half_width = 1.0  # m, this is 'w'
    stacking_period = 5.0   # m, this is 'D'
    J_critical = 1.0e6      # A/m^2, this is 'Jc'
    strip_thickness = 0.1   # m, this is 'd'
    
    # Applied external magnetic field
    H_applied = 4.0e4       # A/m, this is 'Ha'

    calculate_magnetic_field(
        x=x_coord,
        z=z_coord,
        w=strip_half_width,
        D=stacking_period,
        Jc=J_critical,
        d=strip_thickness,
        Ha=H_applied
    )