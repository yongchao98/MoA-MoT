import math

def calculate_magnetic_field():
    """
    Calculates the magnetic field for a stack of superconducting strips.

    The function uses placeholder values for physical quantities and calculates
    the magnetic field based on the critical state model for a fully penetrated
    stack of thin strips. It then prints the full equation with intermediate
    values and the final result.
    """
    # --- Input Parameters (using SI units) ---
    # Jc: Critical current density in A/m^2
    Jc = 1.0e10
    # d: Thickness of the strip in m
    d = 1.0e-6
    # w: Half-width of the strip in m
    w = 1.0e-3
    # D: Stacking interval (distance between strips) in m
    D = 0.2e-3
    # Ha: Applied magnetic field in A/m
    Ha = 5000.0
    # x, z: Coordinates of the point to calculate the field at, in m
    x = 0.5e-3  # |x| >> w is approximated as x being a few times w
    z = 0.0     # Field on the mid-plane between two strips

    # --- Step 1: Calculate the full penetration field H0 ---
    H0 = (Jc * d) / math.pi

    # --- Validate Conditions ---
    if Ha <= H0:
        print(f"Calculation not valid: Applied field Ha ({Ha:.2f} A/m) must be greater than H0 ({H0:.2f} A/m).")
        return
    if abs(x) < w:
        # The condition is |x| >> a, where a is the flux front position.
        # For Ha > H0, a=w. So we check |x| > w.
        print(f"Warning: The condition |x| >> w might not be well satisfied (x={x}, w={w}). Result accuracy may be affected.")


    # --- Step 2: Calculate the magnetic moment per unit length m_z ---
    try:
        atanh_val = math.atanh(H0 / Ha)
    except ValueError:
        print("Error: Invalid argument for atanh. H0/Ha must be < 1.")
        return
        
    m_z = - (math.pi * w**2 / 2) * (Ha - H0 * atanh_val)

    # --- Step 3: Calculate the field from the superconductors, H_sc_z ---
    arg_x = (2 * math.pi * x) / D
    arg_z = (2 * math.pi * z) / D

    cosh_val = math.cosh(arg_x)
    cos_val = math.cos(arg_z)

    # Numerator and Denominator of the main expression for H_sc_z
    numerator = 1 - cosh_val * cos_val
    denominator = (cosh_val - cos_val)**2
    
    # Pre-factor for the H_sc_z expression
    prefactor = (2 * m_z * math.pi) / (D**2)
    
    if abs(denominator) < 1e-12:
        print("Error: Denominator is zero, cannot calculate field.")
        return

    H_sc_z = prefactor * numerator / denominator

    # --- Step 4: Calculate the total magnetic field H_z ---
    H_z = Ha + H_sc_z

    # --- Print the results ---
    print("Calculation of the magnetic field H_z(x, z) for a stack of superconducting strips.")
    print("\n--- Parameters ---")
    print(f"Jc = {Jc:.1e} A/m^2, d = {d:.1e} m, w = {w:.1e} m, D = {D:.1e} m")
    print(f"Ha = {Ha:.1f} A/m")
    print(f"Point (x, z) = ({x:.1e} m, {z:.1e} m)")
    
    print("\n--- Equation ---")
    print("H_z(x,z) = Ha + H_sc_z")
    print("H_sc_z = (2 * m_z * pi / D^2) * (1 - cosh(2*pi*x/D) * cos(2*pi*z/D)) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2")
    
    print("\n--- Step-by-Step Calculation ---")
    print(f"1. H0 = Jc * d / pi = ({Jc:.1e} * {d:.1e}) / {math.pi:.4f} = {H0:.4f} A/m")
    print(f"2. m_z = -(pi * w^2 / 2) * (Ha - H0 * atanh(H0/Ha))")
    print(f"   m_z = -({math.pi:.4f} * ({w:.1e})^2 / 2) * ({Ha:.1f} - {H0:.4f} * atanh({H0/Ha:.4f})) = {m_z:.4e} A*m")
    print(f"3. H_sc_z = (2 * {m_z:.4e} * {math.pi:.4f} / ({D:.1e})^2) * ({numerator:.4e}) / ({denominator:.4e})")
    print(f"   H_sc_z = ({prefactor:.4e}) * ({numerator / denominator:.4e}) = {H_sc_z:.4f} A/m")
    print(f"4. H_z = Ha + H_sc_z = {Ha:.4f} + {H_sc_z:.4f}")

    print("\n--- Final Result ---")
    print(f"H_z = {H_z:.4f} A/m")
    
    # Return final answer in specified format
    # The user did not ask for a single value, but the calculation result.
    # The final numerical value will be returned.
    print(f"\n<<<{H_z:.4f}>>>")

if __name__ == '__main__':
    calculate_magnetic_field()