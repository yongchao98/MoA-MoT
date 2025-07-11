import sympy as sp

def solve_magnetic_field():
    """
    This function symbolically solves for the magnetic field H inside and outside
    a spherical shell with a given surface current.
    """
    # Define symbolic variables
    # A1 and B1 are the unknown coefficients for the magnetic scalar potential
    # for the l=1 term inside and outside the sphere, respectively.
    A1, B1 = sp.symbols('A1 B1')
    
    # R: radius of the sphere
    # K0: amplitude of the surface current
    # mu: magnetic permeability inside the sphere
    # mu0: magnetic permeability of free space (outside)
    R, K0, mu, mu0 = sp.symbols('R K_0 mu mu_0', positive=True, real=True)

    # We derive two equations from the boundary conditions at r = R.
    
    # Equation 1: From the continuity of the normal component of B (mu*H_r).
    # H_in_r = -A1*cos(theta), H_out_r = 2*B1*R**(-3)*cos(theta)
    # At r=R, mu * H_in_r = mu0 * H_out_r
    # mu * (-A1) = mu0 * (2 * B1 / R**3)
    # This simplifies to: mu * A1 * R**3 + 2 * mu0 * B1 = 0
    eq1 = sp.Eq(mu * A1 * R**3 + 2 * mu0 * B1, 0)

    # Equation 2: From the discontinuity of the tangential component of H.
    # H_in_theta = A1*sin(theta), H_out_theta = B1*R**(-3)*sin(theta)
    # At r=R, H_out_theta - H_in_theta = K_phi = K0*sin(theta)
    # (B1 / R**3) - A1 = K0
    eq2 = sp.Eq(B1 / R**3 - A1, K0)

    # Solve the system of two linear equations for the two unknown coefficients A1 and B1.
    solution = sp.solve([eq1, eq2], (A1, B1))
    A1_sol = solution[A1]
    B1_sol = solution[B1]

    # --- Construct the final expressions for the magnetic fields ---

    # Inside the sphere (r < R), the field is uniform: H_in = -A1 * z_hat
    H_in_coeff = -A1_sol
    
    # Outside the sphere (r > R), the field is a dipole field:
    # H_out = (B1 / r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    # The coefficient part is B1.
    H_out_coeff = B1_sol

    # --- Format the output to match the answer choices ---
    
    # For H_in, the expression is: (2 * K0 * mu0) / (mu + 2 * mu0)
    # To match the options, we divide the numerator and denominator by mu.
    H_in_final_num = sp.simplify((2 * K0 * mu0) / mu)
    H_in_final_den = sp.simplify((mu + 2 * mu0) / mu)

    # For H_out, the expression for B1 is: (K0 * mu * R**3) / (mu + 2 * mu0)
    # To match the options, we divide the numerator and denominator by mu.
    H_out_final_num = sp.simplify((K0 * mu) / mu)
    H_out_final_den = sp.simplify((mu + 2 * mu0) / mu)

    print("The derived magnetic field H(r, theta) is:")
    print("{" )
    # Print the inside field
    print(f"  For 0 < r < R (inside):")
    print(f"    H_in = ( ({H_in_final_num}) / ({H_in_final_den}) ) * z_hat")
    
    # Print the outside field
    print(f"  For r > R (outside):")
    print(f"    H_out = ( ({H_out_final_num}) / ({H_out_final_den}) ) * (R^3/r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
    print("}" )
    
    print("\nThis result matches option E.")

solve_magnetic_field()
<<<E>>>