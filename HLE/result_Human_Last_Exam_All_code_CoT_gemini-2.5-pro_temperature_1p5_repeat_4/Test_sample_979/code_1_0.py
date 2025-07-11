import sympy as sp

def solve_magnetic_field_problem():
    """
    This script outlines the step-by-step derivation to find the magnetic field H
    for a spherical shell with a surface current.
    """
    
    print("Step 1: Define potentials from Laplace's Equation.")
    print("Inside (r < R), potential must be finite at r=0: Phi_in = A1 * r * cos(theta)")
    print("Outside (r > R), potential must be zero at infinity: Phi_out = B1 * r**(-2) * cos(theta)")
    print("-" * 20)

    print("Step 2: Define magnetic field H = -grad(Phi).")
    print("H_in = -A1 * z_hat (a uniform field)")
    print("H_out = (B1 / r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat) (a dipole field)")
    print("-" * 20)

    print("Step 3: Apply Boundary Conditions at r = R.")
    print("BC 1 (normal B): mu * H_in_r = mu0 * H_out_r")
    print("   => mu * (-A1*cos(theta)) = mu0 * (2*B1/R**3 * cos(theta))")
    print("   => Equation (I): A1 = -2 * mu0 * B1 / (mu * R**3)")
    print("BC 2 (tangential H): H_out_theta - H_in_theta = K0 * sin(theta)")
    print("   => (B1/R**3 * sin(theta)) - (A1*sin(theta)) = K0*sin(theta)")
    print("   => Equation (II): B1/R**3 - A1 = K0")
    print("-" * 20)

    print("Step 4: Solve for coefficients A1 and B1.")
    print("Substituting (I) into (II):")
    print("B1/R**3 - (-2*mu0*B1/(mu*R**3)) = K0")
    print("B1/R**3 * (1 + 2*mu0/mu) = K0  =>  B1 = (K0 * R**3) / (1 + 2*mu0/mu)")
    print("A1 = - (2*mu0 / (mu*R**3)) * B1 = -(2*mu0*K0) / (mu*(1 + 2*mu0/mu))")
    print("-" * 20)

    print("Step 5: State the final fields and compare with options.")
    
    print("Final Expression for H_in(r, theta) (r < R):")
    # Using H_in = -A1 * z_hat
    print("H_in = - (-(2*mu0*K0) / (mu*(1 + 2*mu0/mu))) * z_hat")
    print("H_in = (2 * mu0 / mu) * (K0 / (1 + 2*mu0/mu)) * z_hat")
    print("")

    print("Final Expression for H_out(r, theta) (r > R):")
    # Using H_out expression with B1
    print("H_out = ((K0 * R**3) / (1 + 2*mu0/mu)) / r**3 * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
    print("H_out = (K0 / (1 + 2*mu0/mu)) * (R**3 / r**3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
    print("-" * 20)
    
    print("Comparing these derived results with the provided choices, the expressions match option E.")

solve_magnetic_field_problem()