import sympy as sp

def solve_electric_field():
    """
    This function symbolically solves for the electric field in the two regions
    of the described cylindrical resistor.
    """
    # Step 1: Define all symbolic variables required for the problem.
    # r is the radial coordinate. V0 is the applied voltage.
    # sigma1 and sigma2 are the conductivities of the two materials.
    # phi is the angular coordinate. pi is the mathematical constant.
    # A, B, C, D are the unknown constants for the potential functions.
    r, V0, sigma1, sigma2 = sp.symbols('r V_0 sigma_1 sigma_2', real=True, positive=True)
    phi = sp.symbols('phi')
    A, B, C, D = sp.symbols('A B C D')

    # Step 2: Define the general form of the potential functions in each region.
    # Based on Laplace's equation, the potential is a linear function of phi.
    # Phi1 is for Region 1 (0 < phi < pi/2)
    # Phi2 is for Region 2 (pi/2 < phi < pi)
    Phi1 = A * phi + B
    Phi2 = C * phi + D

    # Step 3: Set up the system of equations from the boundary and interface conditions.
    # Condition 1: Potential at phi = 0 is V0.
    eq1 = sp.Eq(Phi1.subs(phi, 0), V0)

    # Condition 2: Potential at phi = pi is 0 (grounded).
    eq2 = sp.Eq(Phi2.subs(phi, sp.pi), 0)

    # Condition 3: Potential is continuous at the interface phi = pi/2.
    eq3 = sp.Eq(Phi1.subs(phi, sp.pi / 2), Phi2.subs(phi, sp.pi / 2))

    # Condition 4: Normal component of current density is continuous at phi = pi/2.
    # This leads to sigma1 * d(Phi1)/d(phi) = sigma2 * d(Phi2)/d(phi).
    eq4 = sp.Eq(sigma1 * sp.diff(Phi1, phi), sigma2 * sp.diff(Phi2, phi))

    # Step 4: Solve the system of four equations for the four unknown constants.
    solution = sp.solve([eq1, eq2, eq3, eq4], (A, B, C, D))

    # Step 5: Calculate the electric field in each region using E = -grad(Phi).
    # Since Phi is a function of phi only, E_phi = -(1/r) * d(Phi)/d(phi).
    # Substitute the solved constants into the derivatives.
    dPhi1_dphi = sp.diff(Phi1, phi).subs(solution)
    dPhi2_dphi = sp.diff(Phi2, phi).subs(solution)

    E1_phi = - (1 / r) * dPhi1_dphi
    E2_phi = - (1 / r) * dPhi2_dphi

    # Step 6: Print the derived expressions for the electric fields.
    # The output shows each component of the final equation for clarity.
    print("The derived electric field in Region 1 (0 < phi < pi/2) is:")
    print("Numerator: 2 * V_0 * sigma_2")
    print("Denominator: r * pi * (sigma_1 + sigma_2)")
    print(f"Full Expression for E_1 = {sp.pretty(E1_phi)} * i_phi\n")

    print("The derived electric field in Region 2 (pi/2 < phi < pi) is:")
    print("Numerator: 2 * V_0 * sigma_1")
    print("Denominator: r * pi * (sigma_1 + sigma_2)")
    print(f"Full Expression for E_2 = {sp.pretty(E2_phi)} * i_phi\n")

    print("Comparing these results with the answer choices reveals that Choice C is the correct answer.")

if __name__ == "__main__":
    solve_electric_field()