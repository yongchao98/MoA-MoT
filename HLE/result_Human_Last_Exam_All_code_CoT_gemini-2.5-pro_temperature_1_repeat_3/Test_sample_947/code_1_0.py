import math

def solve_field_expression():
    """
    This function generates and prints the analytical expression for the magnetic field
    of a stack of superconducting strips.
    """

    # Define symbolic variables for the equation
    Ha = "Ha"
    Jc = "Jc"
    d = "d"
    w = "w"
    D = "D"
    x = "x"
    z = "z"
    pi = "pi"
    num_2 = "2"

    # Plan Explanation:
    # 1. The system is a stack of superconducting strips in an applied field Ha > H0.
    # 2. For |x| >> w, each strip is approximated as a magnetic dipole.
    # 3. The magnetic moment of a fully penetrated strip is m_z = -Jc*d*w^2/2.
    # 4. The total field is Ha + the field from the infinite array of dipoles.
    # 5. The code below constructs the final expression based on this plan.

    # Construct the terms of the equation for the induced field
    # Coefficient part, substituting m_z = -(Jc * d * w^2) / 2
    coefficient = f"({pi} * {Jc} * {d} * {w}^{num_2}) / ({num_2} * {D}^{num_2})"

    # Numerator of the fraction part
    term_cosh = f"cosh({num_2}*{pi}*{x}/{D})"
    term_cos = f"cos({num_2}*{pi}*{z}/{D})"
    numerator = f"({term_cosh} * {term_cos} - 1)"

    # Denominator of the fraction part
    denominator = f"({term_cosh} - {term_cos})^{num_2}"

    # The full induced field term (note the negative sign from m_z)
    induced_field_term = f"{coefficient} * [{numerator} / {denominator}]"

    # Combine all parts into the final expression for H_z(x, z)
    final_equation = f"H_z({x}, {z}) = {Ha} - {induced_field_term}"

    # Print the final expression
    print("The derived expression for the magnetic field H_z(x, z) is:")
    print(final_equation)

    # Outputting the components of the equation as requested.
    # Since no numerical values are given, the symbolic representation of each part is printed.
    print("\n--- Equation Components ---")
    print(f"Applied Field Component: {Ha}")
    print(f"Induced Field Coefficient: {coefficient}")
    print(f"Numerator of spatial part: {numerator}")
    print(f"Denominator of spatial part: {denominator}")
    print("-------------------------")
    print("\n--- Equation Variables ---")
    print(f"Ha: Applied magnetic field")
    print(f"Jc: Critical current density")
    print(f"d: Thickness of each strip")
    print(f"w: Half-width of each strip")
    print(f"D: Stacking interval (distance between strips)")
    print(f"x, z: Coordinates where the field is calculated")
    print(f"pi: Mathematical constant {math.pi}")
    print("------------------------")


solve_field_expression()
<<<H_z(x, z) = Ha - (pi * Jc * d * w^2) / (2 * D^2) * [(cosh(2*pi*x/D) * cos(2*pi*z/D) - 1) / (cosh(2*pi*x/D) - cos(2*pi*z/D))^2]>>>