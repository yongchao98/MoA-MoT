import sympy

def solve_mass_problem():
    """
    Calculates the squared mass of the sixth degree of freedom in a modified
    theory of linearized gravity.
    """
    m_sq = sympy.Symbol('m^2')

    print("Step 1: Define coefficients for the traceless tensor part (phi).")
    # Lagrangian for phi: L_phi = A1 * (d_phi)^2 - B1 * phi^2
    # From derivation, L_phi = (1/2)*(d_phi)^2 - (m^2/2)*phi^2
    A1 = sympy.Rational(1, 2)
    B1 = m_sq / 2
    print(f"Kinetic coefficient A1 = {A1}")
    print(f"Mass term coefficient B1 = {B1}")

    M5_sq = B1 / A1
    print("\nStep 2: Calculate the squared mass for the 5 d.o.f.")
    print(f"M_5^2 = B1 / A1 = ({B1}) / ({A1}) = {M5_sq}")
    print("This matches the given information in the problem.")

    print("\nStep 3: Define coefficients for the scalar part (h_bar).")
    # Lagrangian for h_bar: L_h = A2 * (d_h)^2 - B2 * h^2
    # From derivation, L_h = (-1/8)*(d_h)^2 + (3*m^2/8)*h^2
    # So B2 = -3*m^2/8
    A2 = sympy.Rational(-1, 8)
    B2 = -(3 * m_sq / 8)
    print(f"Kinetic coefficient A2 = {A2}")
    print(f"Mass term coefficient B2 = {B2}")

    M6_sq = B2 / A2
    print("\nStep 4: Calculate the squared mass for the 6th d.o.f.")
    # The final equation for the mass M6_sq:
    print(f"M_6^2 = B2 / A2 = ({B2}) / ({A2})")
    print(f"Final Result: The squared mass of the sixth degree of freedom is {M6_sq}.")

solve_mass_problem()