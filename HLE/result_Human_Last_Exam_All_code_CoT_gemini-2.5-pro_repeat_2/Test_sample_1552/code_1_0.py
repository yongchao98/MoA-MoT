from fractions import Fraction

def find_second_heat_kernel_coefficient():
    """
    Calculates the coefficients of the b_2 term in the heat kernel expansion
    for a massless gauged Dirac spinor field in 4 dimensions.

    The local coefficient b_2(x) has the general form:
    b_2(x) = N * (c1 * R^2 + c2 * R_uv*R^uv + c3 * R_uvrs*R^uvrs)
             + c4 * tr(F_uv*F^uv)
    where N is the dimension of the gauge group representation.
    """

    # The calculation is based on the standard formula for the a_2 coefficient:
    # a_2(x, P) = (1/360) * tr[ (5*R^2 - 2*R_ab^2 + 2*R_abcd^2)*Id
    #                          + 60*R*E + 180*E^2 + 30*Omega_ab*Omega^ab ]
    # for an operator P = -Delta - E.
    # We adapt this for P = D^2.

    # Numerators of the coefficients for each term, derived from the trace calculations.
    # The common denominator for all terms is 360.
    denominator = 360

    # Coefficient for N * R^2
    c_R_sq_num = 5 * 4 - 60 + 180 / 4
    # Coefficient for N * R_uv*R^uv
    c_Ric_sq_num = -2 * 4
    # Coefficient for N * R_uvrs*R^uvrs
    c_Riem_sq_num = 2 * 4 - 30 / 2
    # Coefficient for tr(F_uv*F^uv)
    c_F_sq_num = 180 * 4 + 30 * 4

    # Create rational numbers for the final coefficients
    c1 = Fraction(int(c_R_sq_num), denominator)
    c2 = Fraction(int(c_Ric_sq_num), denominator)
    c3 = Fraction(int(c_Riem_sq_num), denominator)
    c4 = Fraction(int(c_F_sq_num), denominator)

    print("The second coefficient in the heat kernel expansion, b_2(x), is given by the equation:")
    print("-" * 80)
    # To avoid potential formatting issues with superscripts, we use descriptive text.
    equation = (
        f"b_2(x) = N * ( {c1} * R^2 + ({c2}) * R_uv*R^uv + ({c3}) * R_uvrs*R^uvrs )"
        f" + ({c4}) * tr(F_uv*F^uv)"
    )
    print(equation)
    print("-" * 80)
    print("Where:")
    print("  N: Dimension of the gauge group representation.")
    print("  R: Ricci scalar curvature.")
    print("  R_uv: Ricci curvature tensor.")
    print("  R_uvrs: Riemann curvature tensor.")
    print("  F_uv: Gauge field strength tensor.")
    print("  tr: Trace over the gauge group generators.")
    
    print("\nThe individual numerical coefficients are:")
    print(f"  Coefficient of N * R^2: {c1}")
    print(f"  Coefficient of N * R_uv*R^uv: {c2}")
    print(f"  Coefficient of N * R_uvrs*R^uvrs: {c3}")
    print(f"  Coefficient of tr(F_uv*F^uv): {c4}")

if __name__ == '__main__':
    find_second_heat_kernel_coefficient()