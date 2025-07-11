import math

def calculate_superconducting_stack_field(Jc, d, w, D, Ha, x, z):
    """
    Calculates the magnetic field for an infinite stack of superconducting strips.

    The calculation is based on a critical state model for thin strips and uses a
    dipole approximation for the field from each strip, summed over the infinite stack.

    Args:
        Jc (float): Critical current density in A/m^2.
        d (float): Thickness of the strips in m.
        w (float): Half-width of the strips in m.
        D (float): Stacking interval (distance between strips) in m.
        Ha (float): Applied magnetic field in the z-direction in A/m.
        x (float): x-coordinate of the observation point in m.
        z (float): z-coordinate of the observation point in m.
    """

    print("--- Input Parameters ---")
    print(f"Critical current density (Jc): {Jc:.2e} A/m^2")
    print(f"Strip thickness (d): {d:.2e} m")
    print(f"Strip half-width (w): {w:.2e} m")
    print(f"Stacking interval (D): {D:.2e} m")
    print(f"Applied field (Ha): {Ha:.2e} A/m")
    print(f"Observation point (x, z): ({x:.2e} m, {z:.2e} m)")
    print("-" * 26)

    # Step 1: Calculate characteristic field H0
    # H0 = Jc * d / pi
    H0 = Jc * d / math.pi
    print("--- Calculated Properties ---")
    print(f"Characteristic field (H0 = Jc*d/pi): {H0:.4f} A/m")

    if Ha <= H0:
        print("\nWarning: The applied field Ha is not greater than H0.")
        print("The model assumes Ha > H0 for this specific expression.")
        # We can still proceed, but the physical assumptions of the prompt are not met.
        # cos(Ha/H0) will be between cos(0)=1 and cos(1)=0.54, so 'a' is well-defined.
        
    # Step 2: Calculate the position of the flux front 'a'
    # a = w * cos(Ha / H0)
    # This relation is valid for Ha <= H0 * pi/2 (full penetration field)
    # The problem specifies Ha > H0.
    cos_arg = Ha / H0
    if cos_arg > math.pi / 2:
        # Full penetration has occurred (a=0). The model for p_x changes.
        # For Ha > Hp, the strip is fully penetrated, and the dipole moment saturates.
        a = 0.0
        print(f"Flux front position (a): {a:.4f} m (Strip is fully penetrated)")
    else:
        a = w * math.cos(cos_arg)
        print(f"Flux front position (a = w*cos(Ha/H0)): {a:.4f} m")

    # Step 3: Calculate the magnetic dipole moment per unit length 'p_x'
    # p_x = Jc * d * (w^2 - a^2)
    K_c = Jc * d
    p_x = K_c * (w**2 - a**2)
    # An alternative expression using Ha: p_x = (pi * H0) * w^2 * sin^2(Ha/H0)
    # p_x_alt = math.pi * H0 * w**2 * (math.sin(min(cos_arg, math.pi/2)))**2
    print(f"Dipole moment per unit length (p_x = Jc*d*(w^2-a^2)): {p_x:.4e} A*m")

    # Step 4: Calculate the field from the infinite array of dipoles
    # Helper variables for the trigonometric/hyperbolic functions
    u = 2 * math.pi * x / D
    v = 2 * math.pi * z / D

    # Denominator in the field expressions
    cosh_u = math.cosh(u)
    cos_v = math.cos(v)
    denominator = cosh_u - cos_v

    if abs(denominator) < 1e-9:
        print("\nError: Observation point is on a singularity (on top of a strip).")
        print("Hx and Hz are divergent.")
        return

    # Geometric factor for Hx component
    # -(1 - cosh(u)*cos(v)) / (cosh(u) - cos(v))^2
    geo_factor_hx = -(1 - cosh_u * cos_v) / (denominator**2)

    # Geometric factor for Hz component
    # sinh(u)*sin(v) / (cosh(u) - cos(v))^2
    sinh_u = math.sinh(u)
    sin_v = math.sin(v)
    geo_factor_hz = (sinh_u * sin_v) / (denominator**2)
    
    # Induced field components
    # H_ind = -p_x * d(H_current_array)/dx
    Hx_induced = -p_x * (math.pi / D**2) * geo_factor_hx
    Hz_induced = -p_x * (math.pi / D**2) * geo_factor_hz
    # There was a sign error in the derivation thought process, let's trace it back.
    # H_z,dipole = -p_x d/dx(H_z_current) and H_x,dipole = -p_x d/dx(H_x_current)
    # H_z,current ~ sin(v)/(cosh(u)-cos(v)) -> d/dx gives ~ sinh(u)sin(v)/(...)^2
    # H_x,current ~ sinh(u)/(cosh(u)-cos(v)) -> d/dx gives ~ (1-cosh*cos)/(...)^2
    # Let's use the final derived formula which seems correct.
    # H_z,ind = p_x * (pi/D^2) * sinh(u)sin(v)/(cosh(u)-cos(v))^2
    # H_x,ind = -p_x * (pi/D^2) * (1-cosh(u)cos(v))/(cosh(u)-cos(v))^2
    Hx_induced_final = -p_x * (math.pi / D**2) * (1 - cosh_u * cos_v) / denominator**2
    Hz_induced_final = p_x * (math.pi / D**2) * sinh_u * sin_v / denominator**2

    print(f"\n--- Field Components ---")
    print(f"Equation for Hx: Hx_induced")
    print(f"Hx = (-p_x * pi / D^2) * (1 - cosh(2*pi*x/D)*cos(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2")
    print(f"Hx = ({-p_x * math.pi / D**2:.4e}) * ({-(1 - cosh_u * cos_v):.4e}) / ({denominator**2:.4e})")
    
    print(f"\nEquation for Hz: Ha + Hz_induced")
    print(f"Hz = Ha + (p_x * pi / D^2) * (sinh(2*pi*x/D)*sin(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2")
    print(f"Hz = {Ha:.4e} + ({p_x * math.pi / D**2:.4e}) * ({sinh_u * sin_v:.4e}) / ({denominator**2:.4e})")

    # Step 5: Calculate total field
    Hx_total = Hx_induced_final
    Hz_total = Ha + Hz_induced_final

    print("\n--- Final Magnetic Field (A/m) ---")
    print(f"Hx = {Hx_total:.4f}")
    print(f"Hz = {Hz_total:.4f}")
    print("-" * 34)


if __name__ == '__main__':
    # --- User-defined parameters ---
    # Material and geometry properties
    Jc_val = 3e10   # Critical current density [A/m^2] (typical for YBCO)
    d_val = 1e-6     # Strip thickness [m] (1 micron)
    w_val = 1e-3     # Strip half-width [m] (1 mm)
    D_val = 5e-3     # Stacking interval [m] (5 mm)

    # Applied field and observation point
    Ha_val = 12000.0 # Applied field [A/m] (approx 150 Oe), chosen to be > H0
    x_val = 2.0e-3   # x-coordinate [m] (2 mm)
    z_val = 1.25e-3  # z-coordinate [m] (z = D/4)
    
    # --- Calculation ---
    calculate_superconducting_stack_field(Jc_val, d_val, w_val, D_val, Ha_val, x_val, z_val)

    # Let's construct the final answer as an expression string
    # Hx(x,z) = - (pi*Jc*d*(w^2 - (w*cos(Ha/H0))^2)/D^2) * [ 1 - cosh(2*pi*x/D)*cos(2*pi*z/D) ] / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2
    # Hz(x,z) = Ha + (pi*Jc*d*(w^2 - (w*cos(Ha/H0))^2)/D^2) * [ sinh(2*pi*x/D)*sin(2*pi*z/D) ] / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2
    # This is too complex to write out elegantly. The Python code is the clearest representation.
    # The user asks for "the expression for the magnetic field".
    final_expr_Hx = "Hx(x,z) = -(pi*Jc*d*w^2*sin(Ha/(Jc*d/pi))^2/D^2) * (1 - cosh(2*pi*x/D)*cos(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2"
    final_expr_Hz = "Hz(x,z) = Ha + (pi*Jc*d*w^2*sin(Ha/(Jc*d/pi))^2/D^2) * (sinh(2*pi*x/D)*sin(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2"
    
    # <<< Hx(x,z) = -(pi*Jc*d*w^2*sin(pi*Ha/(Jc*d))^2/D^2) * (1 - cosh(2*pi*x/D)*cos(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2 ; Hz(x,z) = Ha + (pi*Jc*d*w^2*sin(pi*Ha/(Jc*d))^2/D^2) * (sinh(2*pi*x/D)*sin(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2 >>>
    # This is also messy. Let's provide a slightly more compact version.
    H0_str = "Jc*d/pi"
    p_x_str = f"pi*H0*w^2*sin(Ha/H0)^2" # simplified p_x = Jc*d*w^2*sin(Ha/H0)^2
    final_answer = f"<<< Hx(x,z) = -({p_x_str}) * (pi/D^2) * (1 - cosh(2*pi*x/D)*cos(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2; Hz(x,z) = Ha + ({p_x_str}) * (pi/D^2) * (sinh(2*pi*x/D)*sin(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2 >>>"
