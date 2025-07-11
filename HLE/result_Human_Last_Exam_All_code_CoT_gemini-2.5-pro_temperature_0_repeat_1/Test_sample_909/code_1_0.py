import sympy

def solve_resistor_field():
    """
    This function symbolically derives the electric field in the two regions of the
    cylindrical half-shell resistor.
    """
    # Define the symbolic variables for the problem
    phi, r = sympy.symbols('phi r', real=True, positive=True)
    C1, C2, D1, D2 = sympy.symbols('C1 C2 D1 D2')
    V0, s1, s2 = sympy.symbols('V0 sigma_1 sigma_2', real=True, positive=True)
    pi = sympy.pi

    # Define the general solutions for the potential in each region
    # Phi(phi) = C*phi + D
    Phi1 = C1 * phi + D1
    Phi2 = C2 * phi + D2

    # Set up the system of equations based on boundary and interface conditions

    # 1. Boundary condition at phi = 0: Potential is V0
    eq1 = sympy.Eq(Phi1.subs(phi, 0), V0)

    # 2. Boundary condition at phi = pi: Potential is 0
    eq2 = sympy.Eq(Phi2.subs(phi, pi), 0)

    # 3. Interface condition at phi = pi/2: Potential is continuous
    eq3 = sympy.Eq(Phi1.subs(phi, pi/2), Phi2.subs(phi, pi/2))

    # 4. Interface condition at phi = pi/2: Normal component of current density is continuous
    # J_phi = -sigma * (1/r) * d(Phi)/d(phi). Continuity of J_phi implies continuity of sigma*d(Phi)/d(phi).
    # d(Phi1)/d(phi) is C1, and d(Phi2)/d(phi) is C2.
    eq4 = sympy.Eq(s1 * sympy.diff(Phi1, phi), s2 * sympy.diff(Phi2, phi))

    # Solve the system of four linear equations for the four constants
    solution = sympy.solve([eq1, eq2, eq3, eq4], (C1, C2, D1, D2))

    # The electric field is E = -grad(Phi). For Phi(phi), E = -(1/r) * d(Phi)/d(phi) * phi_hat
    # The magnitude in the phi direction is E_phi = -C/r
    E1_phi = -solution[C1] / r
    E2_phi = -solution[C2] / r

    # --- Output the results ---
    print("The derived electric field in Region 1 (0 < phi < pi/2) is:")
    print("E1 = (Numerator / Denominator) * r_hat")
    
    # Break down the expression for E1
    num1 = 2 * s2 * V0
    den1 = r * pi * (s1 + s2)
    print(f"Numerator: {num1}")
    print(f"Denominator: {den1}")
    print(f"Final Expression for E1: ({E1_phi}) * phi_hat\n")

    print("The derived electric field in Region 2 (pi/2 < phi < pi) is:")
    print("E2 = (Numerator / Denominator) * r_hat")

    # Break down the expression for E2
    num2 = 2 * s1 * V0
    den2 = r * pi * (s1 + s2)
    print(f"Numerator: {num2}")
    print(f"Denominator: {den2}")
    print(f"Final Expression for E2: ({E2_phi}) * phi_hat\n")
    
    print("Comparing these results with the given options, the correct choice is C.")

if __name__ == '__main__':
    solve_resistor_field()