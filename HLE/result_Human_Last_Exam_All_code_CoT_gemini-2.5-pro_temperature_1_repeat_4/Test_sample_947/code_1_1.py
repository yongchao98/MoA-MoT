import sympy

def calculate_magnetic_field_expression():
    """
    This script symbolically derives the expression for the magnetic field H(x)
    for an infinite stack of superconducting strips.
    """
    # Step 1: Define all the necessary symbols for the equation.
    # H(x): The total magnetic field as a function of x
    # Ha: The applied magnetic field
    # H0: The characteristic field for full penetration, defined as Jc*d/pi
    # M0: The magnetic dipole moment per unit length of a single strip
    # w: The half-width of the strip
    # d: The thickness of the strip
    # D: The stacking interval between strips
    # x: The lateral position from the center of the stack
    # Jc: The critical current density
    # H_sc: The magnetic field from the screening currents
    H = sympy.Function('H')(x)
    Ha, H0, w, D, x, Jc, d, M0 = sympy.symbols('H_a H_0 w D x J_c d M_0', real=True, positive=True)
    H_sc = sympy.Symbol('H_sc')

    # Step 2: Define the expression for the magnetic moment per unit length (M0).
    # This is based on the critical state model for a thin strip for Ha > H0.
    M0_expr = sympy.Rational(1, 2) * Ha * w**2 * sympy.tanh(Ha / H0)**2

    # Step 3: Define the expression for the magnetic field from the screening currents (H_sc).
    # This is derived by summing the fields of an infinite lattice of 2D dipoles.
    H_sc_expr = -(M0 * sympy.pi / (2 * D**2)) * (sympy.csch(sympy.pi * x / D))**2

    # Step 4: The total magnetic field H(x) is the sum of the applied field and the screening field.
    H_total_expr = Ha + H_sc_expr

    # Step 5: Substitute the expression for the magnetic moment M0 into the total field equation.
    H_final_expr = H_total_expr.subs(M0, M0_expr)

    # Print the final, complete expression for the total magnetic field H(x).
    # sympy.pprint is used for a more readable, formatted output.
    print("The final expression for the total magnetic field H(x) is:")
    sympy.pprint(sympy.Eq(H, H_final_expr), use_unicode=True)


if __name__ == '__main__':
    calculate_magnetic_field_expression()