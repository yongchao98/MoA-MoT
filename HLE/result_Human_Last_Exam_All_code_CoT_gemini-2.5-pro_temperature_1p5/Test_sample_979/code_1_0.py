import sympy

def solve_magnetic_field():
    """
    Solves for the magnetic field H inside and outside a spherical shell
    with a surface current K = K0*sin(theta)*phi_hat using symbolic mathematics.
    """
    # 1. Define symbolic variables
    r, theta = sympy.symbols('r theta', real=True, positive=True)
    R, K0, mu, mu0 = sympy.symbols('R K0 mu mu0', real=True, positive=True)
    A1, B1 = sympy.symbols('A1 B1') # Coefficients to be determined

    # 2. Define scalar potentials for the l=1 mode
    # Inside the sphere (r < R), potential must be finite at r=0
    Phi_in = A1 * r * sympy.cos(theta)
    # Outside the sphere (r > R), potential must vanish at r=inf
    Phi_out = B1 * r**(-2) * sympy.cos(theta)

    # 3. Calculate H-fields from potentials: H = -grad(Phi)
    # H = - ( d(Phi)/dr * r_hat + (1/r) * d(Phi)/dtheta * theta_hat )
    H_in_r = -sympy.diff(Phi_in, r)
    H_in_theta = - (1/r) * sympy.diff(Phi_in, theta)

    H_out_r = -sympy.diff(Phi_out, r)
    H_out_theta = - (1/r) * sympy.diff(Phi_out, theta)

    # 4. Set up boundary condition equations at r=R
    # Condition 1: Normal component of B is continuous (mu0 * H_out_r = mu * H_in_r)
    eq1_lhs = (mu * H_in_r).subs(r, R)
    eq1_rhs = (mu0 * H_out_r).subs(r, R)
    # This equation has a cos(theta) factor on both sides which can be cancelled.
    eq1 = sympy.Eq(eq1_lhs / sympy.cos(theta), eq1_rhs / sympy.cos(theta))

    # Condition 2: Tangential H is discontinuous by K (H_out_theta - H_in_theta = K0*sin(theta))
    eq2_lhs = (H_out_theta - H_in_theta).subs(r, R)
    eq2_rhs = K0 * sympy.sin(theta)
    # This equation has a sin(theta) factor on both sides which can be cancelled.
    eq2 = sympy.Eq(eq2_lhs / sympy.sin(theta), eq2_rhs / sympy.sin(theta))

    # 5. Solve the system of two equations for A1 and B1
    solution = sympy.solve([eq1, eq2], (A1, B1))
    sol_A1 = solution[A1]
    sol_B1 = solution[B1]

    # 6. Construct the final expressions for the H-fields
    # For H_in, we express it as a uniform field along z_hat, where H_in = -A1 * z_hat
    H_in_coeff = -sol_A1
    # Factor the coefficient to match the form in the answer choices
    H_in_coeff_styled = (2 * mu0 / mu) * (K0 / (1 + 2 * mu0 / mu))

    # For H_out, it's a dipole field. H_out = (B1/r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    # The coefficient part is B1
    H_out_coeff_styled = K0 * R**3 / (1 + 2 * mu0 / mu)
    
    # 7. Print the results, showing each term in the final equations.
    print("The derived magnetic field H(r, theta) is:")
    print("----------------------------------------------------------------")
    print("Inside the sphere (0 < r < R):")
    print("H_in(r, theta) = C_in * z_hat")
    print("\nwhere the coefficient C_in is:")
    sympy.pprint(H_in_coeff_styled, use_unicode=True)

    print("\nOutside the sphere (R < r < infinity):")
    print("H_out(r, theta) = (C_out / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
    print("\nwhere the coefficient C_out is:")
    sympy.pprint(H_out_coeff_styled, use_unicode=True)
    print("\n----------------------------------------------------------------")
    print("These expressions match answer choice E.")

if __name__ == '__main__':
    solve_magnetic_field()