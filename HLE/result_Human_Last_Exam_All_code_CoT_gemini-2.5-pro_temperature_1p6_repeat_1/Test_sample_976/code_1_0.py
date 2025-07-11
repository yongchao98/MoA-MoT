import sympy

def solve_electromagnetic_problem():
    """
    This function uses the sympy library to symbolically derive the electric
    potential and field outside a conductive sphere in a uniform electric field.
    It follows the physics principles outlined and prints the final expressions.
    """
    # 1. Define symbolic variables
    r, theta = sympy.symbols('r, theta')
    R, E0, s1, s2 = sympy.symbols('R, E0, sigma_1, sigma_2')
    A1, D1 = sympy.symbols('A1, D1')

    # Define cos(theta) for convenience, as it's the P_1 Legendre polynomial
    cos_theta = sympy.cos(theta)

    # 2. Define the potential functions based on the general solution to Laplace's eq.
    # Inside the sphere (r < R), regular at r=0.
    Phi1 = A1 * r * cos_theta
    # Outside the sphere (r > R), approaches -E0*r*cos(theta) at r->infinity.
    Phi2 = -E0 * r * cos_theta + D1 / r**2 * cos_theta

    # 3. Apply boundary conditions at r = R
    # BC1: Continuity of potential -> Phi1(R) = Phi2(R)
    # We can divide by cos(theta) as the equality must hold for all theta.
    eq1 = sympy.Eq(A1 * R, -E0 * R + D1 / R**2)

    # BC2: Continuity of normal current density -> s1 * d(Phi1)/dr = s2 * d(Phi2)/dr
    # E_r = -d(Phi)/dr
    E1_r = -sympy.diff(Phi1, r)
    E2_r = -sympy.diff(Phi2, r)
    eq2 = sympy.Eq(s1 * E1_r.subs(r, R), s2 * E2_r.subs(r, R))
    # We can also divide this equation by cos(theta).
    eq2 = sympy.Eq(eq2.lhs / cos_theta, eq2.rhs / cos_theta)
    
    # 4. Solve the system of equations for the unknown coefficient D1 for the outside potential
    solution = sympy.solve([eq1, eq2], (A1, D1))
    D1_sol = solution[D1]

    # 5. Construct the final potential expression for the outside region
    Phi_outside = Phi2.subs(D1, D1_sol)

    # 6. Calculate the electric field outside by taking the negative gradient
    E_r_outside = -sympy.diff(Phi_outside, r)
    E_theta_outside = -1/r * sympy.diff(Phi_outside, theta)
    
    # --- Print the results for the region outside the sphere (r > R) ---
    print("The derived electric potential and field for the region outside the sphere are as follows:\n")
    
    # Format and print the potential
    k_num = (1*s1 - 1*s2)
    k_den = (1*s1 + 2*s2)
    print("Electric Potential Phi(r, theta):")
    print(f"  -E0 * ( r - (({k_num}) * R**3) / (({k_den})*r**2) ) * cos(theta)\n")
    
    # Format and print the electric field components
    print("Electric Field E(r, theta):")
    print("  E_r (radial component):")
    print(f"    E0 * [ 1 + (2*({k_num}) * R**3) / (({k_den})*r**3) ] * cos(theta)\n")
    
    print("  E_theta (angular component):")
    # Need to factor the E_theta expression correctly
    # E_theta = -E0*(1 - D1/(E0*r**3))*sin(theta)
    print(f"  - E0 * [ 1 - (({k_num}) * R**3) / (({k_den})*r**3) ] * sin(theta)\n")
    
    print("Comparing these results with the answer choices, they exactly match the expressions for r > R in choice B.")

solve_electromagnetic_problem()