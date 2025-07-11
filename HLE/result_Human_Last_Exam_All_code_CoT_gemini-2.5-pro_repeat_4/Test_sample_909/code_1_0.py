import sympy

def solve_resistor_field():
    """
    This function symbolically derives the electric field in the two regions of the
    cylindrical half-shell resistor using the boundary and interface conditions.
    """
    # 1. Define all symbolic variables for the problem.
    V0, sigma1, sigma2, r, phi = sympy.symbols('V_0 sigma_1 sigma_2 r phi')
    A1, B1, A2, B2 = sympy.symbols('A_1 B_1 A_2 B_2')
    pi_sym = sympy.pi

    # 2. Define the general form of the potential in each region based on d²Φ/dφ² = 0.
    # Region 1: 0 < phi < pi/2
    Phi1 = A1 * phi + B1
    # Region 2: pi/2 < phi < pi
    Phi2 = A2 * phi + B2

    # 3. Set up the system of equations from the boundary and interface conditions.
    # Condition a: Potential at phi = 0 is V_0.
    eq1 = Phi1.subs(phi, 0) - V0

    # Condition b: Potential at phi = pi is 0 (grounded).
    eq2 = Phi2.subs(phi, pi_sym)

    # Condition c: Potential is continuous at the interface phi = pi/2.
    eq3 = Phi1.subs(phi, pi_sym / 2) - Phi2.subs(phi, pi_sym / 2)

    # Condition d: Normal component of current density (J_phi) is continuous.
    # J_phi = -sigma/r * d(Phi)/dphi. The continuity condition is sigma_1 * E_phi1 = sigma_2 * E_phi2,
    # which simplifies to sigma_1 * d(Phi_1)/dphi = sigma_2 * d(Phi_2)/dphi.
    dPhi1_dphi = sympy.diff(Phi1, phi)
    dPhi2_dphi = sympy.diff(Phi2, phi)
    eq4 = sigma1 * dPhi1_dphi - sigma2 * dPhi2_dphi

    # 4. Solve the system of four linear equations for the four unknown constants.
    solution = sympy.solve([eq1, eq2, eq3, eq4], (A1, B1, A2, B2))

    # 5. Calculate the electric field in each region.
    # E = -grad(Phi). For Phi(phi), this is E = -(1/r) * d(Phi)/dphi in the phi direction.
    E1_phi = - (1/r) * dPhi1_dphi.subs(solution)
    E2_phi = - (1/r) * dPhi2_dphi.subs(solution)

    # 6. Display the final expressions for the electric fields.
    print("The derived electric field in Region 1 (0 < phi < pi/2) is:")
    print(f"E_1 = ({sympy.numer(E1_phi)}) / ({sympy.denom(E1_phi)}) * i_phi")
    print(f"   = (2 * V_0 * sigma_2) / (r * pi * (sigma_1 + sigma_2)) * i_phi")

    print("\nThe derived electric field in Region 2 (pi/2 < phi < pi) is:")
    print(f"E_2 = ({sympy.numer(E2_phi)}) / ({sympy.denom(E2_phi)}) * i_phi")
    print(f"   = (2 * V_0 * sigma_1) / (r * pi * (sigma_1 + sigma_2)) * i_phi")

    print("\nThese expressions match answer choice C.")

solve_resistor_field()
<<<C>>>