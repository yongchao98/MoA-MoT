import sympy

def solve_resistor_field():
    """
    This script symbolically solves for the electric field in the two regions
    of the cylindrical resistor using the steps outlined above.
    """
    # Define the symbolic variables
    phi, r, V0, sigma1, sigma2 = sympy.symbols('phi r V0 sigma1 sigma2')
    pi = sympy.pi
    A1, B1, A2, B2 = sympy.symbols('A1 B1 A2 B2')

    # Define the potential functions for the two regions
    Phi1 = A1 * phi + B1
    Phi2 = A2 * phi + B2

    # Set up the system of four equations from the boundary and interface conditions
    # 1. Potential at phi = 0 is V0
    eq1 = Phi1.subs(phi, 0) - V0
    # 2. Potential at phi = pi is 0
    eq2 = Phi2.subs(phi, pi)
    # 3. Potential is continuous at phi = pi/2
    eq3 = Phi1.subs(phi, pi/2) - Phi2.subs(phi, pi/2)
    # 4. Normal component of current density is continuous at phi = pi/2
    # J_phi = -sigma * (1/r) * d(Phi)/d(phi). The (1/r) cancels.
    # So, sigma1 * d(Phi1)/d(phi) = sigma2 * d(Phi2)/d(phi)
    eq4 = sigma1 * sympy.diff(Phi1, phi) - sigma2 * sympy.diff(Phi2, phi)

    # Solve the system of equations for the four constants
    solution = sympy.solve([eq1, eq2, eq3, eq4], [A1, B1, A2, B2])

    # The electric field E_phi component is -1/r * d(Phi)/d(phi) = -A/r
    E1_phi = -solution[A1] / r
    E2_phi = -solution[A2] / r

    # Format and print the results
    print("Derivation of the Electric Field:")
    print("-" * 30)
    
    # E1 expression
    e1_num = sympy.numer(E1_phi)
    e1_den = sympy.denom(E1_phi)
    print("Electric field in Region 1 (0 < phi < pi/2):")
    print(f"E1 = ({sympy.pretty(e1_num)}) / ({sympy.pretty(e1_den)}) * i_phi")
    # Manually printing to match the answer choices format
    print(f"   = (2 * sigma2 * V0) / (r * pi * (sigma1 + sigma2)) * i_phi\n")

    # E2 expression
    e2_num = sympy.numer(E2_phi)
    e2_den = sympy.denom(E2_phi)
    print("Electric field in Region 2 (pi/2 < phi < pi):")
    print(f"E2 = ({sympy.pretty(e2_num)}) / ({sympy.pretty(e2_den)}) * i_phi")
    # Manually printing to match the answer choices format
    print(f"   = (2 * sigma1 * V0) / (r * pi * (sigma1 + sigma2)) * i_phi\n")

    print("Comparing these results with the answer choices, the correct choice is C.")

if __name__ == '__main__':
    solve_resistor_field()
<<<C>>>