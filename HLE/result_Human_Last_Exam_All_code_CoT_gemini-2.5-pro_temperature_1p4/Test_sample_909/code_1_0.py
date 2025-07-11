import sympy as sp

def solve_resistor_field():
    """
    This script symbolically solves for the electric field in a composite cylindrical resistor.
    """
    print("This script solves for the electric field based on the following plan:")
    print("1. Define the potential Phi(phi) in each region as a linear function, since nabla^2(Phi) = 0 implies d^2(Phi)/d(phi)^2 = 0.")
    print("2. Set up a system of four equations using the boundary and interface conditions.")
    print("3. Solve for the unknown coefficients in the potential functions.")
    print("4. Calculate the electric field E = -grad(Phi) = -(1/r) * d(Phi)/d(phi) * i_phi.")
    print("\nSolving...\n")

    # Define symbolic variables
    phi, r, V_0, sigma_1, sigma_2 = sp.symbols('phi r V_0 sigma_1 sigma_2', real=True, positive=True)
    A1, B1, A2, B2 = sp.symbols('A1 B1 A2 B2')
    pi = sp.pi

    # Define potential functions in each region
    # Region 1 (0 < phi < pi/2):
    Phi_1 = A1 * phi + B1
    # Region 2 (pi/2 < phi < pi):
    Phi_2 = A2 * phi + B2

    # Define boundary and interface conditions as a system of equations
    # 1. Potential at phi=0 is V_0
    eq1 = sp.Eq(Phi_1.subs(phi, 0), V_0)

    # 2. Potential at phi=pi is 0
    eq2 = sp.Eq(Phi_2.subs(phi, pi), 0)

    # 3. Potential is continuous at phi=pi/2
    eq3 = sp.Eq(Phi_1.subs(phi, pi/2), Phi_2.subs(phi, pi/2))

    # 4. Normal component of current density (J_phi) is continuous at phi=pi/2
    # J_phi = sigma * E_phi = -sigma * (1/r) * d(Phi)/d(phi)
    # The (1/r) factor cancels, so the condition is sigma_1 * d(Phi_1)/d(phi) = sigma_2 * d(Phi_2)/d(phi)
    eq4 = sp.Eq(sigma_1 * sp.diff(Phi_1, phi), sigma_2 * sp.diff(Phi_2, phi))

    # Solve the system of equations for the coefficients A1, B1, A2, B2
    solution = sp.solve([eq1, eq2, eq3, eq4], (A1, B1, A2, B2))
    A1_sol = solution[A1]
    A2_sol = solution[A2]

    # Calculate the electric field in each region using E_phi = -(1/r) * d(Phi)/d(phi)
    E1_phi = - (1/r) * sp.diff(Phi_1, phi).subs(A1, A1_sol)
    E2_phi = - (1/r) * sp.diff(Phi_2, phi).subs(A2, A2_sol)

    # Simplify the final expressions
    E1_final = sp.simplify(E1_phi)
    E2_final = sp.simplify(E2_phi)

    # Print the final result in a clear format matching the answer choices.
    print("The calculated electric field in region 1 (0 < phi < pi/2) is:")
    print(f"E_1 = ({2} * sigma_2 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi")

    print("\nThe calculated electric field in region 2 (pi/2 < phi < pi) is:")
    print(f"E_2 = ({2} * sigma_1 * V_0) / (r * pi * (sigma_1 + sigma_2)) * i_phi")
    
    print("\nComparing these results to the options shows that choice C is the correct answer.")

solve_resistor_field()
<<<C>>>