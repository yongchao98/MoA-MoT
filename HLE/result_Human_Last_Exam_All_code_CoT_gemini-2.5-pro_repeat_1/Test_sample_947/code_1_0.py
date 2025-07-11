import math

def calculate_magnetic_field(Jc, d, w, D, Ha, x, z):
    """
    Calculates the magnetic field H=(Hx, Hz) for a stack of superconducting strips.

    The calculation is valid for an applied field Ha > H0 and far from the
    strips' center (|x| >> a).

    Args:
        Jc (float): Critical current density (A/m^2).
        d (float): Thickness of the strips (m).
        w (float): Half-width of the strips (m).
        D (float): Stacking period (m).
        Ha (float): Applied magnetic field in the z-direction (A/m).
        x (float): x-coordinate for field calculation (m).
        z (float): z-coordinate for field calculation (m).
    """
    # Step 1: Define constants and check inputs
    if x == 0:
      # sgn(0) is 0, so Hx will be 0, which is correct.
      # However, to avoid potential division by zero in other models,
      # it's good practice to handle this edge case.
      # The current model is safe at x=0.
      pass

    # Step 2: Calculate characteristic field H0
    H0 = Jc * d / math.pi
    print(f"Step 1: Calculate characteristic field H0 = Jc * d / pi")
    print(f"H0 = {Jc} * {d} / {math.pi:.4f} = {H0:.4f} A/m\n")


    if Ha <= H0:
        print(f"Warning: The applied field Ha={Ha:.4f} is not greater than H0={H0:.4f}.")
        print("The model assumptions might not be met.")
        # Continue calculation for demonstration purposes
        
    # Step 3: Calculate the flux penetration depth 'a'
    # Ha = H0 * ln[sinh(pi*w/2D) / sinh(pi*a/2D)]
    # => sinh(pi*a/2D) = sinh(pi*w/2D) * exp(-Ha/H0)
    try:
        sinh_arg_w = math.pi * w / (2 * D)
        sinh_w_term = math.sinh(sinh_arg_w)
        
        exp_term = math.exp(-Ha / H0)
        
        # Argument for asinh can't be > 1 if a<w, but due to floating point
        # it might be slightly larger, or negative if Ha is very large.
        asinh_arg = sinh_w_term * exp_term
        if asinh_arg < 0:
             print(f"Error: Cannot compute 'a' because argument to asinh is negative ({asinh_arg:.2e}). This can happen for very large Ha/H0.")
             return

        sinh_a_term = math.asinh(asinh_arg)
        a = (2 * D / math.pi) * sinh_a_term
    except (ValueError, OverflowError) as e:
        print(f"Error calculating 'a': {e}")
        return

    print(f"Step 2: Calculate flux front position 'a' from Ha")
    print(f"a = (2*D/pi) * asinh[sinh(pi*w/2D) * exp(-Ha/H0)]")
    print(f"a = (2*{D}/{math.pi:.4f}) * asinh[sinh({sinh_arg_w:.4f}) * exp(-{Ha/H0:.4f})]")
    print(f"a = ({2*D/math.pi:.4f}) * asinh[{sinh_w_term:.4f} * {exp_term:.4f}]")
    print(f"a = ({2*D/math.pi:.4f}) * asinh[{asinh_arg:.4f}] = {a:.4f} m\n")
    
    # Step 4: Calculate the terms for the field expression
    # Amplitude Factor: Amp = H0 * (cosh(2*pi*w/D) - cosh(2*pi*a/D))
    cosh_arg_w = 2 * math.pi * w / D
    cosh_arg_a = 2 * math.pi * a / D
    cosh_w_term = math.cosh(cosh_arg_w)
    cosh_a_term = math.cosh(cosh_arg_a)
    Amplitude = H0 * (cosh_w_term - cosh_a_term)

    # Spatial Decay Factor: Decay = exp(-2*pi*|x|/D)
    decay_arg = -2 * math.pi * abs(x) / D
    Decay = math.exp(decay_arg)

    # Z-dependence Factors
    z_arg = 2 * math.pi * z / D
    CosZ = math.cos(z_arg)
    SinZ = math.sin(z_arg)
    
    # Step 5: Calculate the final field components
    # Note: math.copysign(1.0, x) gives the sign of x (-1 for x<0, 1 for x>=0)
    Hx = -math.copysign(1.0, x) * Amplitude * Decay * SinZ
    Hz_induced = Amplitude * Decay * CosZ
    Hz = Ha + Hz_induced

    # Step 6: Print the results in expression form
    print("Step 3: Calculate field components Hx and Hz using far-field expressions.\n")
    print("--- Magnetic Field Hx ---")
    print("Expression: Hx = -sgn(x) * H0 * (cosh(2*pi*w/D) - cosh(2*pi*a/D)) * exp(-2*pi*|x|/D) * sin(2*pi*z/D)")
    print(f"Hx = -sgn({x}) * {H0:.4f} * (cosh({cosh_arg_w:.4f}) - cosh({cosh_arg_a:.4f})) * exp({decay_arg:.4f}) * sin({z_arg:.4f})")
    print(f"Hx = -({math.copysign(1.0, x):.0f}) * {H0:.4f} * ({cosh_w_term:.4f} - {cosh_a_term:.4f}) * {Decay:.4f} * {SinZ:.4f}")
    print(f"Hx = {Hx:.4f} A/m\n")

    print("--- Magnetic Field Hz ---")
    print("Expression: Hz = Ha + H0 * (cosh(2*pi*w/D) - cosh(2*pi*a/D)) * exp(-2*pi*|x|/D) * cos(2*pi*z/D)")
    print(f"Hz = {Ha} + {H0:.4f} * ({cosh_w_term:.4f} - {cosh_a_term:.4f}) * {Decay:.4f} * {CosZ:.4f}")
    print(f"Hz = {Ha} + ({Hz_induced:.4f})")
    print(f"Hz = {Hz:.4f} A/m")


if __name__ == '__main__':
    # --- Example Parameters ---
    # Material and Geometry
    Jc = 1e10    # Critical current density in A/m^2
    d = 1e-6     # Strip thickness in m
    w = 1e-3     # Strip half-width in m (2w = 2mm)
    D = 5e-3     # Stacking period in m (D = 5mm)

    # Applied Field and Position
    Ha = 4000.0  # Applied field in A/m (H0 for these params is ~3183 A/m)
    x = 1e-2     # x-coordinate (1 cm, |x|>>w)
    z = 1.25e-3  # z-coordinate (z=D/4)
    
    print("--- Input Parameters ---")
    print(f"Jc = {Jc:.2e} A/m^2")
    print(f"d = {d:.2e} m")
    print(f"w = {w:.2e} m")
    print(f"D = {D:.2e} m")
    print(f"Ha = {Ha:.2f} A/m")
    print(f"Position (x, z) = ({x:.2e}, {z:.2e}) m\n")

    calculate_magnetic_field(Jc, d, w, D, Ha, x, z)