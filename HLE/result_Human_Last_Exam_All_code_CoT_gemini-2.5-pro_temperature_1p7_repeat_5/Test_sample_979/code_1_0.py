import sys

def solve_electromagnetism_problem():
    """
    This function outlines the step-by-step derivation of the magnetic field
    and prints the final correct answer.
    """
    
    # Using a placeholder for mu_0 to avoid Unicode issues in some terminals
    mu0_char = "u0"
    mu_char = "u"
    
    # Step 1 & 2: Define potentials based on symmetry
    # Inside (r < R), potential must be regular at r=0. Due to sin(theta) symmetry, only l=1 term is needed.
    # Phi_in = A * r * cos(theta)
    # Outside (r > R), potential must vanish at infinity.
    # Phi_out = B * (1/r^2) * cos(theta)

    # Step 3: Apply boundary conditions to get equations for coefficients A and B
    # Let H_in = -grad(Phi_in) and H_out = -grad(Phi_out)
    
    # BC1: mu * H_r_in = mu0 * H_r_out at r=R
    # mu * (-A * cos(theta)) = mu0 * (2 * B * (1/R^3) * cos(theta))
    # -> Eq 1: mu * A = -2 * mu0 * B / R^3
    
    # BC2: H_theta_out - H_theta_in = K0 * sin(theta) at r=R
    # (B * (1/R^3) * sin(theta)) - (A * sin(theta)) = K0 * sin(theta)
    # -> Eq 2: B / R^3 - A = K0
    
    # Step 4: Solve the system of equations
    # From Eq 2: B / R^3 = A + K0
    # Substitute into Eq 1:
    # mu * A = -2 * mu0 * (A + K0)
    # mu * A = -2 * mu0 * A - 2 * mu0 * K0
    # A * (mu + 2*mu0) = -2 * mu0 * K0
    # A = - (2 * mu0 * K0) / (mu + 2*mu0)
    
    # Solve for B:
    # B = R^3 * (A + K0) = R^3 * (K0 - (2 * mu0 * K0) / (mu + 2*mu0))
    # B = K0 * R^3 * (1 - 2*mu0 / (mu + 2*mu0))
    # B = K0 * R^3 * ((mu + 2*mu0 - 2*mu0) / (mu + 2*mu0))
    # B = (mu * K0 * R^3) / (mu + 2*mu0)

    # Step 5: Determine the fields
    
    # Inside field (r < R):
    # H_in is uniform. Phi_in = A * z, so H_in = -grad(A*z) = -A * z_hat.
    # H_in_z = -A = (2 * mu0 * K0) / (mu + 2*mu0)
    # We can rewrite the coefficient by dividing numerator and denominator by mu:
    # H_in_z = (2 * K0 * (mu0/mu)) / (1 + 2*(mu0/mu))
    
    h_in_str = f"H_in = ( (2 * {mu0_char}) / {mu_char} ) / ( 1 + (2 * {mu0_char})/{mu_char} ) * K0 * z_hat"

    # Outside field (r > R):
    # H_out = -grad(B * cos(theta) / r^2) which gives a dipole field.
    # H_out = (B / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    # Coefficient: B = (mu * K0 * R^3) / (mu + 2*mu0)
    # Divide numerator and denominator by mu:
    # B = (K0 * R^3) / (1 + 2*mu0/mu)
    
    h_out_str = f"H_out = ( K0 / (1 + (2*{mu0_char})/{mu_char}) ) * (R^3/r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)"

    # Step 6: Final Answer Formulation
    print("Based on the derivation, the magnetic fields are:")
    print("-" * 50)
    
    print("For the region inside the sphere (0 < r < R):")
    print(f"    H(r, theta) = ( (2 * \u03BC\u2080 / \u03BC) / (1 + 2 * \u03BC\u2080 / \u03BC) ) * K\u2080 * z_hat")
    print("\nWhich is a uniform field in the z-direction.")

    print("\nFor the region outside the sphere (R < r < \u221E):")
    print(f"    H(r, theta) = ( K\u2080 / (1 + 2 * \u03BC\u2080 / \u03BC) ) * (R\u00B3/r\u00B3) * (2*cos(\u03B8)*r_hat + sin(\u03B8)*\u03B8_hat)")
    print("\nThis corresponds to a magnetic dipole field.")
    print("-" * 50)
    print("This result matches option E.")


solve_electromagnetism_problem()

# The final result in the format of the options:
print("\nFinal Answer:")
print(
"""
H(r, theta) =
    /  2 * u0/u                 
    |  ---------- * K0 * z_hat                                               ,  0 < r < R
    <  1 + 2*u0/u
    |
    |      K0                 R^3
    \  ---------- * (2*cos(theta)*r_hat + sin(theta)*theta_hat) * ---       ,  R < r < infinity
       1 + 2*u0/u              r^3
"""
)
<<<E>>>