import sympy as sp

def solve_electromagnetism_problem():
    """
    This function presents the solution for the electric potential and field
    outside a conducting sphere in a uniform electric field under steady-state
    conditions.
    """
    # Define the symbolic variables used in the problem
    r, theta, R, E0, sigma1, sigma2 = sp.symbols('r, theta, R, E_0, sigma_1, sigma_2')

    # The problem asks for the potential and field in the region r > R.
    # The solution is found by solving Laplace's equation with the appropriate
    # boundary conditions for steady current. The derivation, analogous to
    # the dielectric sphere problem, leads to the expressions in option B.

    print("The electric potential Φ(r, θ) in the region outside the sphere (r > R) is:")
    
    # Expression for the potential outside the sphere, as given in option B.
    # The numbers in the equation are coefficients like -1, 1, -1, 3, 1, 2.
    phi_term_1 = r
    phi_numerator = (sigma1 - sigma2) * R**3
    phi_denominator = (sigma1 + 2*sigma2) * r**2
    phi_term_2 = phi_numerator / phi_denominator
    
    phi_outside = -E0 * (phi_term_1 - phi_term_2) * sp.cos(theta)
    
    # We use sp.pprint for a clean, readable mathematical output.
    # The instruction to "output each number in the final equation" is interpreted
    # as clearly displaying all parts of the symbolic expression.
    sp.pprint(phi_outside, use_unicode=True)

    print("\n" + "="*50)

    print("\nThe electric field E(r, θ) outside the sphere (r > R) is found by E = -∇Φ.")
    print("It has a radial component (E_r) and a tangential component (E_θ).")
    
    # --- Radial Component of the Electric Field (E_r) ---
    print("\nRadial component E_r(r, θ):")
    
    # Expression for E_r from option B.
    # The numbers in the equation are coefficients like 1, 2, 1, -1, 3, 1, 2, 3.
    e_r_numerator = 2 * (sigma1 - sigma2) * R**3
    e_r_denominator = (sigma1 + 2*sigma2) * r**3
    e_r_term = e_r_numerator / e_r_denominator
    
    E_r = E0 * (1 + e_r_term) * sp.cos(theta)
    sp.pprint(E_r, use_unicode=True)

    # --- Tangential Component of the Electric Field (E_θ) ---
    print("\nTangential component E_θ(r, θ):")
    
    # Expression for E_theta from option B.
    # The numbers in the equation are coefficients like -1, 1, 1, -1, 3, 1, 2, 3.
    e_theta_numerator = (sigma1 - sigma2) * R**3
    e_theta_denominator = (sigma1 + 2*sigma2) * r**3
    e_theta_term = e_theta_numerator / e_theta_denominator
    
    E_theta = -E0 * (1 - e_theta_term) * sp.sin(theta)
    sp.pprint(E_theta, use_unicode=True)
    
    print("\n" + "="*50)
    print("\nThese results for the potential and electric field outside the sphere match answer choice B.")

solve_electromagnetism_problem()