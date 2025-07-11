import sympy

def find_heat_kernel_coefficient():
    """
    This function calculates the local density of the second Seeley-DeWitt
    coefficient, a_2(x), for a massless gauged Dirac spinor field in 4D.

    The calculation is based on standard formulas from the heat kernel expansion
    literature. The Dirac operator squared (D^2) is a Laplace-type operator,
    and its a_2 coefficient can be expressed as a sum of terms involving
    the curvature of the spacetime manifold and the gauge field.

    The formula used is (from I.G. Avramidi's work):
    a_2(x) = (1/360) * Tr_V[-10*C^2 - 13*E] - (1/12) * Tr_V(Omega^2)
    where:
    - Tr_V is the trace over the vector bundle (spinor and gauge spaces).
    - C^2 is the square of the Weyl tensor.
    - E is the integrand of the Gauss-Bonnet term.
    - Omega is the total curvature of the connection on the bundle.

    We will express this in the basis of squared Riemann tensor (R_abcd^2),
    squared Ricci tensor (R_ab^2), squared Ricci scalar (R^2), and the
    squared gauge field strength (F^2).
    """

    # Define symbolic variables using sympy
    R_abcd2 = sympy.Symbol('R_μνρσ^2')  # Squared Riemann tensor
    R_ab2 = sympy.Symbol('R_μν^2')      # Squared Ricci tensor
    R2 = sympy.Symbol('R^2')            # Squared Ricci scalar
    F2 = sympy.Symbol('tr(F_μν*F^μν)') # Squared gauge field strength term
    N_F = sympy.Symbol('N_F')           # Dimension of the gauge representation

    # Part 1: Contribution from the generic formula for a minimal Laplacian
    # a_2_generic = (1/360) * Tr_V[-10*C^2 - 13*E]
    # We convert C^2 and E to the R_abcd2, R_ab2, R2 basis.
    # C^2 = R_abcd2 - 2*R_ab2 + (1/3)*R2
    # E = R_abcd2 - 4*R_ab2 + R2
    # Tr_V = 4 * N_F for a Dirac spinor in 4D.

    coeff_factor_1 = (4 * N_F) / 360

    # Coefficient of C^2 is -10
    term_C2 = -10 * (R_abcd2 - 2*R_ab2 + sympy.Rational(1, 3)*R2)

    # Coefficient of E is -13
    term_E = -13 * (R_abcd2 - 4*R_ab2 + R2)

    part1 = coeff_factor_1 * (term_C2 + term_E)

    # Part 2: Contribution from the curvature of the connection
    # a_2_omega = -(1/12) * Tr_V(Omega^2)
    # For a gauged Dirac operator, Tr_V(Omega^2) = 2*N_F*R_abcd2 + 4*F2

    part2_R = -sympy.Rational(1, 12) * (2 * N_F * R_abcd2)
    part2_F = -sympy.Rational(1, 12) * (4 * F2)

    # Combine all parts
    a2_density = sympy.simplify(part1 + part2_R + part2_F)

    # Expand the expression to show all the terms and coefficients clearly
    a2_expanded = sympy.expand(a2_density)

    # Extract coefficients for printing
    coeff_R_abcd2 = a2_expanded.coeff(R_abcd2)
    coeff_R_ab2 = a2_expanded.coeff(R_ab2)
    coeff_R2 = a2_expanded.coeff(R2)
    coeff_F2 = a2_expanded.coeff(F2)

    print("The local density of the second heat kernel coefficient, a_2(x), is given by the formula:")
    print(f"a_2(x) = ({coeff_R_abcd2}) * R_μνρσ^2 + ({coeff_R_ab2}) * R_μν^2 + ({coeff_R2}) * R^2 + ({coeff_F2}) * tr(F_μν*F^μν)")

    print("\nWhere:")
    print(" - N_F is the dimension of the gauge representation of the spinor field.")
    print(" - R_μνρσ^2, R_μν^2, R^2 are the squares of the Riemann tensor, Ricci tensor, and Ricci scalar, respectively.")
    print(" - tr(F_μν*F^μν) is the trace of the squared gauge field strength, which for a semi-simple gauge group is Tr(T^a T^b)F_μν^a F_b^μν.")


if __name__ == '__main__':
    find_heat_kernel_coefficient()