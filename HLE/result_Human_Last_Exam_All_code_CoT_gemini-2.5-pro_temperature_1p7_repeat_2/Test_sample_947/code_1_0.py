import math

def calculate_superconducting_field():
    """
    Calculates the magnetic field for a stack of superconducting strips.

    This function implements the formula for the magnetic field H(x) on the
    surface of a thin superconducting strip, which is part of an infinite
    stack. The calculation is performed under the high-field approximation
    (Ha > H0) where the flux penetration is nearly complete (a -> 0).

    The formula is:
    H(x) = Ha - (Jc * d / pi) * arccosh(sinh(pi * w / D) / sinh(pi * |x| / D))
    """

    # --- 1. Define physical parameters with example values (in SI units) ---
    # Jc: Critical current density in Amperes per square meter (A/m^2)
    Jc = 1.0e10
    # d: Thickness of the strip in meters (m)
    d = 1.0e-6
    # w: Half-width of the strip in meters (m)
    w = 1.0e-3
    # D: Stacking interval (center-to-center distance) in meters (m)
    D = 5.0e-3
    # Ha: Applied magnetic field in Amperes per meter (A/m)
    Ha = 20000.0
    # x: Position along the strip width to calculate the field (m)
    # This must satisfy 0 < |x| < w
    x = 0.5e-3
    
    # --- 2. Print the formula and parameters ---
    print("This script calculates the magnetic field H(x) for a stack of superconducting strips.")
    print("The governing equation for the high-field limit is:")
    print("H(x) = Ha - (Jc * d / pi) * arccosh(sinh(pi * w / D) / sinh(pi * |x| / D))\n")
    
    print("--- Given Parameters ---")
    print(f"Applied Field (Ha) = {Ha} A/m")
    print(f"Critical Current Density (Jc) = {Jc:.1e} A/m^2")
    print(f"Strip Thickness (d) = {d:.1e} m")
    print(f"Strip Half-Width (w) = {w:.1e} m")
    print(f"Stacking Interval (D) = {D:.1e} m")
    print(f"Position of Interest (x) = {x:.1e} m")
    print(f"Value of pi = {math.pi}\n")
    
    # --- 3. Step-by-step calculation ---
    print("--- Calculation Breakdown ---")

    # Check validity of x
    if not (0 < abs(x) < w):
        print(f"Error: The position x={x} must be within the strip boundaries (0 < |x| < {w}).")
        return

    # Calculate the pre-factor H0 = Jc * d / pi
    H0 = (Jc * d) / math.pi
    print(f"1. Calculate the characteristic field H0 = (Jc * d / pi):")
    print(f"   H0 = ({Jc:.1e} * {d:.1e} / {math.pi:.5f}) = {H0:.2f} A/m")

    # Calculate the argument of the top sinh
    arg_sinh_top = (math.pi * w) / D
    val_sinh_top = math.sinh(arg_sinh_top)
    print(f"\n2. Calculate the numerator inside arccosh: sinh(pi * w / D)")
    print(f"   sinh({math.pi:.5f} * {w:.1e} / {D:.1e}) = sinh({arg_sinh_top:.5f}) = {val_sinh_top:.5f}")

    # Calculate the argument of the bottom sinh
    arg_sinh_bot = (math.pi * abs(x)) / D
    val_sinh_bot = math.sinh(arg_sinh_bot)
    print(f"\n3. Calculate the denominator inside arccosh: sinh(pi * |x| / D)")
    print(f"   sinh({math.pi:.5f} * {abs(x):.1e} / {D:.1e}) = sinh({arg_sinh_bot:.5f}) = {val_sinh_bot:.5f}")

    # Calculate the argument of arccosh
    if val_sinh_top <= val_sinh_bot:
        print(f"\nError: The argument for arccosh must be > 1. This requires |x| < w.")
        return
    arg_acosh = val_sinh_top / val_sinh_bot
    print(f"\n4. Calculate the argument of arccosh:")
    print(f"   arg = {val_sinh_top:.5f} / {val_sinh_bot:.5f} = {arg_acosh:.5f}")
    
    # Calculate the value of arccosh
    val_acosh = math.acosh(arg_acosh)
    print(f"\n5. Calculate the arccosh value:")
    print(f"   arccosh({arg_acosh:.5f}) = {val_acosh:.5f}")

    # Calculate the total shielding field
    H_shielding = H0 * val_acosh
    print(f"\n6. Calculate the total shielding field term H_shielding = H0 * arccosh(...)")
    print(f"   H_shielding = {H0:.2f} * {val_acosh:.5f} = {H_shielding:.2f} A/m")

    # --- 4. Final Result ---
    H_final = Ha - H_shielding
    print("\n--- Final Result ---")
    print("The final equation with values is:")
    print(f"H({x}) = {Ha} - {H_shielding:.2f}")
    print(f"H({x}) = {H_final:.2f} A/m")
    
    # The final numerical answer in the required format
    print(f"\n<<<H({x}) = {H_final:.2f}>>>")

# Execute the function
calculate_superconducting_field()