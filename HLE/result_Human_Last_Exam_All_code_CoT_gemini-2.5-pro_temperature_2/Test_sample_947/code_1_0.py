import math

def print_magnetic_field_formula():
    """
    This function prints the derived symbolic expression for the magnetic field H_z(x, z).
    The formula is for an infinite stack of superconducting strips in a transverse magnetic field Ha,
    under the assumption of full penetration.
    """

    # The formula is derived by summing the contribution from each strip in the stack.
    # The field from a single fully-penetrated strip is known. The sum over the
    # infinite stack is evaluated using a standard mathematical identity for summing logarithms.

    # Final formula for H_z(x, z)
    # Note: We substitute Jc*d with π*H0 as given in the problem description.
    # So, the prefactor (Jc * d / (4 * pi)) becomes (pi * H0 / (4 * pi)) = H0 / 4.

    print("The expression for the magnetic field H_z for the given conditions is:")
    print("")
    # Using characters to build the formula string for clarity
    numerator = "sinh^2(π*(x+w)/D) + sin^2(π*z/D)"
    denominator = "sinh^2(π*(x-w)/D) + sin^2(π*z/D)"

    # Print the equation part by part to ensure all numbers and symbols are included as requested.
    print("H_z(x, z) = H_a - (H0 / 4) * ln[({}) / ({})]".format(numerator, denominator))
    print("")
    print("Where:")
    print("  H_a is the applied magnetic field.")
    print("  H0 = Jc * d / π is the characteristic field.")
    print("  Jc is the critical current density.")
    print("  d is the strip thickness.")
    print("  2w is the strip width.")
    print("  D is the stacking interval.")
    print("  x and z are the spatial coordinates.")
    print("  ln is the natural logarithm.")

# Execute the function to print the formula
print_magnetic_field_formula()
