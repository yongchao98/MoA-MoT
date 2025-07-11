import math

def print_magnetic_field_formula():
    """
    Prints the derived formula for the magnetic field H_z(x, z)
    for a stack of superconducting strips.
    """
    # Define variables as strings for printing the formula
    Ha = "Ha"
    Jc = "Jc"
    d = "d"
    a = "a"
    w = "w"
    D = "D"
    x = "x"
    z = "z"
    pi = "pi"

    # The expression for the magnetic field Hz(x,z) for |x| >> a is derived
    # by considering the far-field limit where the current distribution on
    # each strip acts as a magnetic dipole.
    # The total field is the sum of the applied field and the induced field from the
    # infinite stack of these dipoles.

    # Format the equation string
    equation = (
        f"H_z({x}, {z}) = {Ha} + ({Jc} * {d} * ({a} - {w}) / {D}) * "
        f"sinh(2 * {pi} * {x} / {D}) / "
        f"(cosh(2 * {pi} * {x} / {D}) - cos(2 * {pi} * {z} / {D}))"
    )

    print("The expression for the magnetic field is:")
    print(equation)

# Execute the function to print the formula
print_magnetic_field_formula()

# The final answer is the expression printed by the code.
final_expression = "Ha + (Jc * d * (a - w) / D) * sinh(2 * pi * x / D) / (cosh(2 * pi * x / D) - cos(2 * pi * z / D))"
<<<Ha + (Jc * d * (a - w) / D) * sinh(2 * pi * x / D) / (cosh(2 * pi * x / D) - cos(2 * pi * z / D))>>>