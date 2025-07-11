import math

def print_magnetic_field_expression():
    """
    Calculates and prints the mathematical expressions for the magnetic field
    components H_x and H_z for a stack of superconducting strips.
    
    The expressions are valid under the conditions:
    - Applied field Ha > H0 = Jc*d/pi (full penetration)
    - Far-field approximation |x| >> w, D
    """

    # The expressions are derived symbolically. This code prints the final formulas.
    # The common factor in the induced field expressions
    # Note: We substitute Jc*d = pi*H0
    common_factor_str = "(2 * pi**2 * H0 * w**2 / D**2)"
    
    # Exponential decay term, common to both components
    exp_decay_str = "exp(-2 * pi * |x| / D)"

    # H_x component expression
    # This component arises solely from the induced currents.
    hx_expression = f"Hx(x, z) = - {common_factor_str} * sign(x) * sin(2*pi*z/D) * {exp_decay_str}"

    # H_z component expression
    # This is the sum of the applied field Ha and the induced z-component.
    hz_expression = f"Hz(x, z) = Ha + {common_factor_str} * {exp_decay_str}"
    
    print("The derived expressions for the magnetic field components Hx and Hz are:")
    print("-" * 70)
    print("In these expressions:")
    print("  Ha: Applied magnetic field in the z-direction.")
    print("  H0: Characteristic field, H0 = Jc * d / pi.")
    print("  w: Half-width of the strips.")
    print("  D: Spacing between the strips.")
    print("  x, z: Coordinates where the field is measured.")
    print("  |x|: Absolute value of x.")
    print("  sign(x): Sign function of x (+1 for x>0, -1 for x<0).")
    print("  pi: The mathematical constant " + str(math.pi) + "...")
    print("-" * 70)
    
    print("\nThe x-component of the magnetic field is:")
    print(hx_expression)
    
    print("\nThe z-component of the magnetic field is:")
    print(hz_expression)

# Execute the function to print the expressions
print_magnetic_field_expression()

# Final answer representation as requested.
# The answer is the set of expressions derived. Since a single output is expected,
# we format them into one string.
final_answer_str = "Hx(x, z) = - (2 * pi**2 * H0 * w**2 / D**2) * sign(x) * sin(2*pi*z/D) * exp(-2*pi*|x|/D); Hz(x, z) = Ha + (2 * pi**2 * H0 * w**2 / D**2) * exp(-2*pi*|x|/D)"
# The prompt is a bit ambiguous about how to represent the "answer" so the code prints a readable version,
# and this provides a single string format of the result. For a single-letter or single-number response this is not suitable.
# Thus, I will return the expressions.
# <<<'Hx(x, z) = - (2 * pi**2 * H0 * w**2 / D**2) * sign(x) * sin(2*pi*z/D) * exp(-2*pi*|x|/D)', 'Hz(x, z) = Ha + (2 * pi**2 * H0 * w**2 / D**2) * exp(-2*pi*|x|/D)'>>>