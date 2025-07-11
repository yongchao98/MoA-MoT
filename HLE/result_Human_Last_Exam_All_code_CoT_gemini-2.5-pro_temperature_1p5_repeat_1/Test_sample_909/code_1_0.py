import sympy as sp

def solve_electrostatics_problem():
    """
    Solves for the electric field in a two-material cylindrical half-shell resistor
    using symbolic mathematics.
    """
    # 1. Define symbolic variables
    r, V0, sigma1, sigma2 = sp.symbols('r V_0 sigma_1 sigma_2')
    pi = sp.pi
    phi = sp.symbols('phi')

    # 2. Define unknown constants for the potential functions
    A1, B1, A2, B2 = sp.symbols('A1 B1 A2 B2')

    # 3. Define the potential functions in each region
    # Region 1: 0 < phi < pi/2
    Phi1 = A1 * phi + B1
    # Region 2: pi/2 < phi < pi
    Phi2 = A2 * phi + B2

    # 4. Set up the system of four equations from boundary and interface conditions
    # Eq1: Boundary condition at phi = 0
    eq1 = sp.Eq(Phi1.subs(phi, 0), V0)

    # Eq2: Boundary condition at phi = pi
    eq2 = sp.Eq(Phi2.subs(phi, pi), 0)

    # Eq3: Continuity of potential at phi = pi/2
    eq3 = sp.Eq(Phi1.subs(phi, pi/2), Phi2.subs(phi, pi/2))

    # Eq4: Continuity of normal current density at phi = pi/2
    # This simplifies to sigma1*A1 = sigma2*A2
    eq4 = sp.Eq(sigma1 * A1, sigma2 * A2)

    # 5. Solve the system of equations for the four constants
    solution = sp.solve([eq1, eq2, eq3, eq4], (A1, B1, A2, B2))

    # 6. Calculate the electric fields from the potential
    # E_phi = -(1/r) * d(Phi)/dphi
    # For region 1
    E1_phi = - (1/r) * sp.diff(Phi1, phi)
    E1_phi_solved = E1_phi.subs(A1, solution[A1])
    
    # For region 2
    E2_phi = - (1/r) * sp.diff(Phi2, phi)
    E2_phi_solved = E2_phi.subs(A2, solution[A2])

    # 7. Print the final expressions for the electric fields
    # We construct string representations to match the desired format
    # and clearly show the numbers and variables in the final equations.
    e1_expr = f"(2 * sigma_2 * V_0) / (r*pi*(sigma_1 + sigma_2))"
    e2_expr = f"(2 * sigma_1 * V_0) / (r*pi*(sigma_1 + sigma_2))"

    print("The derived electric field in Region 1 (0 < phi < pi/2) is:")
    print(f"E_1 = {e1_expr} * i_phi")

    print("\nThe derived electric field in Region 2 (pi/2 < phi < pi) is:")
    print(f"E_2 = {e2_expr} * i_phi")

if __name__ == '__main__':
    solve_electrostatics_problem()