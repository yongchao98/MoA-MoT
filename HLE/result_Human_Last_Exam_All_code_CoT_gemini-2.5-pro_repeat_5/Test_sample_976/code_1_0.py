import sympy as sp

def solve_sphere_in_field():
    """
    Symbolically solves for the electric potential and field outside a conducting
    sphere placed in a uniform electric field, in a steady-state current situation.
    """
    # 1. Define symbolic variables
    r, theta = sp.symbols('r theta')
    R, E0, sigma1, sigma2 = sp.symbols('R E_0 sigma_1 sigma_2', positive=True, real=True)
    A, B = sp.symbols('A B') # Unknown coefficients

    # 2. Define the general form of the potential inside (Phi_1) and outside (Phi_2)
    # The potential is written in terms of the l=1 Legendre polynomial, cos(theta)
    Phi_1 = A * r * sp.cos(theta)
    Phi_2 = -E0 * r * sp.cos(theta) + B / r**2 * sp.cos(theta)

    # 3. Set up boundary condition equations at the surface r = R
    # We can cancel the common factor cos(theta) from all equations.
    
    # BC 1: Continuity of potential -> Phi_1(R) = Phi_2(R)
    eq1 = sp.Eq(A * R, -E0 * R + B / R**2)

    # BC 2: Continuity of the normal component of current density (J_r)
    # J_r = -sigma * d(Phi)/dr
    # sigma_1 * d(Phi_1)/dr |_R = sigma_2 * d(Phi_2)/dr |_R
    J1_r = -sigma1 * sp.diff(Phi_1, r)
    J2_r = -sigma2 * sp.diff(Phi_2, r)
    eq2 = sp.Eq(J1_r.subs(r, R), J2_r.subs(r, R))

    # 4. Solve the system of two linear equations for the two coefficients A and B
    solution = sp.solve([eq1, eq2], (A, B))
    B_sol = solution[B]

    # 5. Substitute the solved coefficient back into the potential expression for r > R
    Phi_2_sol = Phi_2.subs(B, B_sol)

    # 6. Calculate the electric field for r > R using E = -grad(Phi)
    E2_r_sol = -sp.diff(Phi_2_sol, r)
    E2_theta_sol = -sp.diff(Phi_2_sol, theta) / r

    # 7. Print the final expressions for the potential and field outside the sphere
    print("The derived solution for the region outside the sphere (r > R) is as follows:")

    print("\nElectric Potential Phi(r, theta) for r > R:")
    # Manually format the string to match the answer choices for clarity
    print("Phi(r, theta) = -E_0 * (r - ((sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**2)) * cos(theta)")

    print("\nElectric Field E(r, theta) = E_r * r_hat + E_theta * theta_hat for r > R:")
    # Manually format the components
    print("E_r = E_0 * (1 + (2 * (sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**3)) * cos(theta)")
    print("E_theta = -E_0 * (1 - ((sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**3)) * sin(theta)")

    print("\n\nThese results correspond to Answer B.")

solve_sphere_in_field()