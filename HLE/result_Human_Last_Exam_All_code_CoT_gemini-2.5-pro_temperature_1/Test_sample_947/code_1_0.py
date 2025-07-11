import math

def display_magnetic_field_expression():
    """
    This function prints the analytical expression for the magnetic field H_z(x, z)
    for a stack of superconducting strips under the specified conditions.
    """

    # Define the variables for clarity in the equation string.
    # These are symbolic representations.
    Ha = "Ha"   # Applied magnetic field
    Jc = "Jc"   # Critical current density
    d = "d"     # Strip thickness
    w = "w"     # Strip half-width
    D = "D"     # Separation between strips
    x = "x"     # x-coordinate
    z = "z"     # z-coordinate
    pi = "π"    # The mathematical constant pi

    # Construct the numerator and denominator of the field contribution from the strips
    numerator_strips = f"{pi}*{Jc}*{d}*{w}² * sinh(2*{pi}*|{x}|/{D}) * sin(2*{pi}*{z}/{D})"
    denominator_strips = f"2*{D}² * (sinh²({pi}*|{x}|/{D}) + sin²({pi}*|{z}|/{D}))²"

    # The final expression for the total magnetic field H_z
    # The negative sign comes from the direction of the induced magnetic moment.
    equation = f"H_z({x}, {z}) = {Ha} - ({numerator_strips}) / ({denominator_strips})"

    print("The expression for the z-component of the magnetic field H_z(x, z) is:")
    print(equation)

if __name__ == '__main__':
    display_magnetic_field_expression()