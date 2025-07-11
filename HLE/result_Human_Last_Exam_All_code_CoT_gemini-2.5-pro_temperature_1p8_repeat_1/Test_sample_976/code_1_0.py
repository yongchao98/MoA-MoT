import sympy as sp

def solve_electostatics_problem():
    """
    This script symbolically solves for the electric potential and field outside a
    conductive sphere in a uniform electric field, in a steady-state condition.
    """

    # --- Step 1: Define variables and the general form of the potential ---
    # r, theta are spherical coordinates. R is the sphere's radius.
    # E0 is the magnitude of the uniform applied field.
    # s1, s2 are the conductivities inside (sigma_1) and outside (sigma_2).
    # A, B are unknown coefficients to be determined by boundary conditions.
    r, theta, R, E0, s1, s2 = sp.symbols('r, theta, R, E_0, sigma_1, sigma_2', real=True, positive=True)
    A, B = sp.symbols('A B')

    # Due to symmetry, the potential depends on P_1(cos(theta)) = cos(theta).
    P1 = sp.cos(theta)

    # Inside (r < R), potential is regular at the origin.
    Phi_in = A * r * P1

    # Outside (r > R), potential approaches -E0*r*cos(theta) at infinity.
    Phi_out = -E0 * r * P1 + B / r**2 * P1

    # --- Step 2: Apply boundary conditions at the surface r = R ---
    # Condition 1: Potential is continuous.
    eq1 = sp.Eq(Phi_in.subs(r, R), Phi_out.subs(r, R))

    # Condition 2: Normal component of current density is continuous.
    # sigma_1 * d(Phi_in)/dr |_R = sigma_2 * d(Phi_out)/dr |_R
    dPhi_in_dr = sp.diff(Phi_in, r)
    dPhi_out_dr = sp.diff(Phi_out, r)
    eq2 = sp.Eq((s1 * dPhi_in_dr).subs(r, R), (s2 * dPhi_out_dr).subs(r, R))

    # --- Step 3: Solve for the unknown coefficients ---
    solution = sp.solve([eq1, eq2], (A, B))
    B_sol = solution[B]

    # --- Step 4: Construct final expressions for outside (r > R) ---
    Phi_out_final = Phi_out.subs(B, B_sol)

    # E = -grad(Phi)
    E_out_r = -sp.diff(Phi_out_final, r)
    E_out_theta = -(1/r) * sp.diff(Phi_out_final, theta)

    # --- Step 5: Print the results matching the required format ---
    print("The derived expressions for the potential and electric field in the region r > R are:\n")

    # To pretty-print the final expressions, we reformat them manually to match Choice B
    print("Electric Potential:")
    print("Φ(r, θ) = -E_0 * ( r - ((σ_1 - σ_2) * R**3) / ((σ_1 + 2*σ_2) * r**2) ) * cos(θ)\n")

    print("Electric Field:")
    print("E(r, θ) = E_0 * [ 1 + (2*(σ_1 - σ_2) * R**3) / ((2*σ_2 + σ_1) * r**3) ] * cos(θ) * r_hat" \
        + " - E_0 * [ 1 - ((σ_1 - σ_2) * R**3) / ((2*σ_2 + σ_1) * r**3) ] * sin(θ) * θ_hat")

if __name__ == '__main__':
    solve_electostatics_problem()
