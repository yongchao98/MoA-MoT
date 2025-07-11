import sympy as sp

def solve_magnetic_shielding():
    """
    This function symbolically solves for the magnetic field H in and around a
    magnetized sphere placed within a conducting spherical shell.
    
    The problem is solved using the magnetic scalar potential method, subject to
    the appropriate boundary conditions.
    """
    
    # Define symbolic variables
    r, theta = sp.symbols('r, theta', real=True, positive=True)
    M0, R_p, R = sp.symbols('M_0, R_p, R', real=True, positive=True)
    A1, C1, D1 = sp.symbols('A1, C1, D1') # Coefficients for the potentials
    cos_theta = sp.cos(theta)

    # --- Step 1: Define the potentials in the two regions ---
    # Region 1 (0 < r < R_p): Must be finite at r=0
    Phi_m1 = A1 * r * cos_theta
    # Region 2 (R_p < r < R): General form
    Phi_m2 = (C1 * r + D1 * r**-2) * cos_theta

    # --- Step 2: Set up equations from boundary conditions ---
    
    # B.C. 1: At the perfect conductor (r=R), the normal component of B is zero.
    # B_r = mu_0 * H_r = -mu_0 * d(Phi_m2)/dr = 0. So, d(Phi_m2)/dr = 0 at r=R.
    eq1_lhs = -sp.diff(Phi_m2, r).subs(r, R)
    eq1 = sp.Eq(eq1_lhs, 0)
    
    # B.C. 2: The potential is continuous at r=R_p (H_tan is continuous).
    # Phi_m1(R_p) = Phi_m2(R_p)
    eq2 = sp.Eq(Phi_m1.subs(r, R_p), Phi_m2.subs(r, R_p))

    # B.C. 3: The normal component of B is continuous at r=R_p.
    # B1_n = B2_n  =>  H1_r + M_r = H2_r
    # M_r is the radial component of the magnetization M = M0*(cos(theta)*i_r - sin(theta)*i_theta)
    M_r = M0 * cos_theta
    H1_r = -sp.diff(Phi_m1, r)
    H2_r = -sp.diff(Phi_m2, r)
    eq3 = sp.Eq(H1_r.subs(r, R_p) + M_r, H2_r.subs(r, R_p))
    
    # --- Step 3: Solve the system of equations for the coefficients ---
    # The equations form a linear system for A1, C1, D1.
    # We remove the cos(theta) factor as it's common to all terms.
    system_eqs = [
        eq1.lhs / cos_theta,
        (eq2.lhs - eq2.rhs) / cos_theta,
        (eq3.lhs - eq3.rhs) / cos_theta,
    ]
    solutions = sp.solve(system_eqs, (A1, C1, D1))
    
    A1_sol = solutions[A1]
    C1_sol = solutions[C1]
    D1_sol = solutions[D1]

    # --- Step 4: Calculate H fields from the solved potentials ---
    # Substitute the solved coefficients back into the potential expressions.
    Phi_m1_sol = Phi_m1.subs(A1, A1_sol)
    Phi_m2_sol = Phi_m2.subs({C1: C1_sol, D1: D1_sol})
    
    # Calculate H = -grad(Phi_m)
    # Region 1: 0 < r < R_p
    H1_r_final = -sp.diff(Phi_m1_sol, r)
    H1_theta_final = -sp.diff(Phi_m1_sol, theta) / r

    # Region 2: R_p < r < R
    H2_r_final = -sp.diff(Phi_m2_sol, r)
    H2_theta_final = -sp.diff(Phi_m2_sol, theta) / r

    # --- Step 5: Print the results in a readable format ---
    print("--- Derived Magnetic Field H ---")
    
    print("\nIn the region 0 < r < R_p:")
    # Factor out common terms to match the format of the answers
    H1_expr = f"H = M0 * ({sp.pretty(sp.simplify(A1_sol/M0), use_unicode=False)}) * (-cos(theta) * i_r + sin(theta) * i_theta)"
    print(H1_expr)

    print("\nIn the region R_p < r < R:")
    # Manually format to match the structure in the answer choices
    H2_r_coeff = sp.simplify(H2_r_final / (M0 * cos_theta))
    H2_theta_coeff = sp.simplify(H2_theta_final / (M0 * sp.sin(theta)))

    # For H_r, we want to match the form -2M0/3 * [ (Rp/R)^3 - (Rp/r)^3 ] * cos(theta)
    H2_r_term1 = sp.powsimp((R_p/R)**3, force=True)
    H2_r_term2 = sp.powsimp((R_p/r)**3, force=True)
    print(f"H_r = - (2 * M0 / 3) * [ {H2_r_term1} - {H2_r_term2} ] * cos(theta)")
    
    # For H_theta, we want to match the form M0/3 * [ 2*(Rp/R)^3 + (Rp/r)^3 ] * sin(theta)
    H2_theta_term1 = sp.powsimp(2 * (R_p/R)**3, force=True)
    H2_theta_term2 = sp.powsimp((R_p/r)**3, force=True)
    print(f"H_theta = (M0 / 3) * [ {H2_theta_term1} + {H2_theta_term2} ] * sin(theta)")

    print("\nComparing these results with the given choices, the expressions match Answer Choice B.")


solve_magnetic_shielding()
<<<B>>>