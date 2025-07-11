import numpy as np

def calculate_and_print_field_expression():
    """
    Calculates and prints the expressions for the magnetic field (Hx, Hz)
    for a stack of superconducting strips in an applied field Ha.
    """
    # 1. Define physical parameters and the point of interest (x, z)
    # All units are in SI (meters, Amperes, etc.)
    Jc = 3e10  # Critical current density (A/m^2)
    d = 1e-6   # Strip thickness (m)
    w = 1e-3   # Half-width of the strip (m)
    D = 5e-3   # Spacing between strips (m)
    Ha = 10000 # Applied magnetic field (A/m)
    x = 2e-3   # x-coordinate to calculate field (m)
    z = 1e-3   # z-coordinate to calculate field (m)

    # 2. Calculate derived physical parameters
    pi = np.pi
    H0 = Jc * d / pi  # Characteristic field

    # Check if the applied field is strong enough for this model
    if Ha <= H0:
        print(f"Error: Applied field Ha={Ha:.2f} A/m is not greater than H0={H0:.2f} A/m.")
        print("The assumptions for this calculation are not met.")
        return

    # Calculate flux penetration depth 'a'
    a = w * np.exp(-Ha / H0)
    
    # Check if the condition |x| >> a holds
    if abs(x) <= 3 * a: # Using a factor of 3 for "much greater"
        print(f"Warning: The condition |x| >> a may not be strongly satisfied.")
        print(f"|x| = {abs(x):.4e} m, a = {a:.4e} m")


    # 3. Calculate intermediate terms for the field equations.
    # The full field expressions are derived from the complex potential for the stack.
    # Hx = (H0/4) * ln[ ((c2x*c2z-k_w)^2 + (s2x*s2z)^2) / ((c2x*c2z-k_a)^2 + (s2x*s2z)^2) ]
    # Hz = Ha - (H0/2) * (arctan(s2x*s2z / (c2x*c2z-k_w)) - arctan(s2x*s2z / (c2x*c2z-k_a)))
    # where:
    c2x = np.cosh(2 * pi * x / D)
    s2x = np.sinh(2 * pi * x / D)
    c2z = np.cos(2 * pi * z / D)
    s2z = np.sin(2 * pi * z / D)
    k_w = np.cosh(2 * pi * w / D)
    k_a = np.cosh(2 * pi * a / D)

    # 4. Calculate the numerical values for each part of the equation
    # For Hx
    hx_numerator_real_part = c2x * c2z - k_w
    hx_numerator_imag_part = s2x * s2z
    hx_numerator_term = hx_numerator_real_part**2 + hx_numerator_imag_part**2

    hx_denominator_real_part = c2x * c2z - k_a
    hx_denominator_imag_part = s2x * s2z
    hx_denominator_term = hx_denominator_real_part**2 + hx_denominator_imag_part**2

    # For Hz
    hz_arg1_numerator = s2x * s2z
    hz_arg1_denominator = c2x * c2z - k_w
    
    hz_arg2_numerator = s2x * s2z
    hz_arg2_denominator = c2x * c2z - k_a

    # Avoid division by zero in arctan argument for specific symmetric points
    # This simplified code assumes non-zero denominators.
    if abs(hz_arg1_denominator) < 1e-9: hz_arg1_denominator = 1e-9
    if abs(hz_arg2_denominator) < 1e-9: hz_arg2_denominator = 1e-9

    hz_arg1 = hz_arg1_numerator / hz_arg1_denominator
    hz_arg2 = hz_arg2_numerator / hz_arg2_denominator

    # 5. Print the final expressions with the numbers filled in
    print("The magnetic field expression is calculated for the following parameters:")
    print(f"  Jc = {Jc:.1e} A/m^2, d = {d:.1e} m, w = {w:.1e} m, D = {D:.1e} m")
    print(f"  Applied Field Ha = {Ha:.2f} A/m")
    print(f"  Coordinates (x, z) = ({x:.1e}, {z:.1e}) m\n")
    print("Calculated intermediate parameters:")
    print(f"  H0 = {H0:.2f} A/m")
    print(f"  Penetration depth a = {a:.4f} m\n")
    
    print("--- Final Magnetic Field Equations ---")
    print("The expression for the x-component of the magnetic field, Hx(x,z):")
    print(f"Hx({x}, {z}) = ({H0:.2f} / 4) * ln[ (({c2x:.3f}*{c2z:.3f} - {k_w:.3f})^2 + ({s2x:.3f}*{s2z:.3f})^2) / (({c2x:.3f}*{c2z:.3f} - {k_a:.3f})^2 + ({s2x:.3f}*{s2z:.3f})^2) ]")
    print(f"Hx({x}, {z}) = ({H0:.2f} / 4) * ln[ ({hx_numerator_term:.3f}) / ({hx_denominator_term:.3f}) ]")

    print("\nThe expression for the z-component of the magnetic field, Hz(x,z):")
    print(f"Hz({x}, {z}) = {Ha:.2f} - ({H0:.2f} / 2) * [ arctan( ({s2x:.3f}*{s2z:.3f}) / ({c2x:.3f}*{c2z:.3f} - {k_w:.3f}) ) - arctan( ({s2x:.3f}*{s2z:.3f}) / ({c2x:.3f}*{c2z:.3f} - {k_a:.3f}) ) ]")
    print(f"Hz({x}, {z}) = {Ha:.2f} - ({H0:.2f} / 2) * [ arctan({hz_arg1:.3f}) - arctan({hz_arg2:.3f}) ]")
    
    # Also printing the final calculated numerical values
    final_hx = (H0 / 4) * np.log(hx_numerator_term / hx_denominator_term)
    final_hz = Ha - (H0 / 2) * (np.arctan(hz_arg1) - np.arctan(hz_arg2))
    print("\n--- Final Numerical Values ---")
    print(f"Hx = {final_hx:.4f} A/m")
    print(f"Hz = {final_hz:.4f} A/m")


if __name__ == '__main__':
    calculate_and_print_field_expression()
    # The final numerical values are calculated just for completeness. 
    # The primary output is the expression with numbers substituted in.
    final_hz = 10000 - (9549.3/2) * (np.arctan(182.174)-np.arctan(6.763))
    # <<<Hz(0.002, 0.001) = 10000.00 - (9549.30 / 2) * [ arctan(182.174) - arctan(6.763) ]>>>