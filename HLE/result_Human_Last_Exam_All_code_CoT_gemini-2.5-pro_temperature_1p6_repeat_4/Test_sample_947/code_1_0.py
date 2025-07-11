def print_magnetic_field_expression():
    """
    This function prints the derived symbolic expression for the magnetic field
    H_z(x, z) for a stack of superconducting strips in an applied field Ha,
    valid in the far-field approximation |x| >> w.
    """

    print("Based on the standard model for a stack of superconducting strips, the z-component of the magnetic field H_z(x, z) is the sum of the applied field (Ha) and the field induced by the screening currents.\n")
    print("The final expression, valid for the far-field region (|x| >> w), is presented below.")
    print("----------------------------------------------------------------------------------\n")

    # The final equation is printed step-by-step to clearly show all parts,
    # as requested by the prompt to "output each number in the final equation".
    
    print("The final equation for the magnetic field H_z(x, z) is:")
    
    # We print the formula string by building it from its components.
    
    # Term 1: Applied Field
    part1 = "Ha"
    
    # Term 2: Amplitude of the induced field
    # The amplitude consists of several symbolic parts.
    part2_H0 = "H0"
    part2_sinh1 = "sinh((2 * pi * w) / D)"
    # The second sinh term contains the penetration depth 'a'
    part2_sinh2 = "sinh((2 * pi / D) * a)"
    part2_amplitude = f"({part2_H0} * ({part2_sinh1} - {part2_sinh2}))"

    # Term 3: Exponential decay with x
    part3_exp = "exp((-2 * pi * |x|) / D)"

    # Term 4: Periodic variation with z
    part4_cos = "cos((2 * pi * z) / D)"

    # Full expression for Hz
    print(f"H_z(x, z) = {part1} - {part2_amplitude} * {part3_exp} * {part4_cos}\n")

    print("where 'a' is the penetration half-width, which depends on the applied field Ha:\n")

    # Expression for 'a'
    part_a_cosh1 = "cosh((pi * w) / D)"
    part_a_cosh2 = "cosh(Ha / H0)"
    print(f"a = (D / pi) * acosh({part_a_cosh1} / {part_a_cosh2})\n")
    
    print("And the parameters are:")
    print("  Ha: Applied magnetic field")
    print("  x, z: Spatial coordinates")
    print("  w: Half-width of the strips")
    print("  D: Stacking interval between strips")
    print("  d: Thickness of each strip")
    print("  Jc: Critical current density")
    print("  H0: Characteristic field, defined as H0 = (Jc * d) / pi")
    print("  pi: Mathematical constant pi (approx. 3.14159)")
    print("  exp, sinh, cosh, acosh, cos: Standard mathematical functions")


if __name__ == '__main__':
    print_magnetic_field_expression()
    # The derived mathematical expression for the magnetic field is:
    final_expression = "H_z(x,z) = Ha - ( (Jc*d/pi) * (sinh(2*pi*w/D) - sinh(2*acosh(cosh(pi*w/D)/cosh(pi*Ha/(Jc*d))))) ) * exp(-2*pi*abs(x)/D) * cos(2*pi*z/D)"
    # The final result is the formula printed above. For a single "answer", we provide the symbolic formula string.
    # The format below requires a single value. Since this is symbolic, the full string is the answer.
    # To conform to the format, I'll provide the formula for the induced field amplitude C(Ha), which is the core of the derivation.
    final_answer_string = "H0 * (sinh(2*pi*w/D) - sinh(2*acosh(cosh(pi*w/D) / cosh(Ha/H0))))"
    # Given the complexity, returning the final complex formula string.
    print(f"\n<<<{final_expression}>>>")
