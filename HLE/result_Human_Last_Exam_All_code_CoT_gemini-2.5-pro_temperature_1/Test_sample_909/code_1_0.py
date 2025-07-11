import sympy

def solve_resistor_field():
    """
    This function symbolically solves for the electric field in the two regions
    of the described cylindrical half-shell resistor.
    """
    # Define the symbolic variables
    V0, sigma1, sigma2, r = sympy.symbols('V_0 sigma_1 sigma_2 r')
    A1, B1, A2, B2 = sympy.symbols('A_1 B_1 A_2 B_2')
    phi = sympy.Symbol('phi')
    pi = sympy.pi

    # Step 1: Define the potential functions for the two regions
    # The general solution to Laplace's equation d²Φ/dφ² = 0 is Φ(φ) = Aφ + B
    Phi_1 = A1 * phi + B1
    Phi_2 = A2 * phi + B2

    # Step 2: Set up the equations based on the boundary conditions
    # BC 1: Potential at phi=0 is V0
    eq1 = sympy.Eq(Phi_1.subs(phi, 0), V0)

    # BC 2: Potential at phi=pi is 0 (grounded)
    eq2 = sympy.Eq(Phi_2.subs(phi, pi), 0)

    # BC 3: Potential is continuous at the interface phi=pi/2
    eq3 = sympy.Eq(Phi_1.subs(phi, pi/2), Phi_2.subs(phi, pi/2))

    # BC 4: Normal component of current density is continuous at phi=pi/2
    # J = σE, E_φ = -(1/r) * dΦ/dφ. dΦ₁/dφ = A₁, dΦ₂/dφ = A₂
    # σ₁ * (-A₁/r) = σ₂ * (-A₂/r)  =>  σ₁A₁ = σ₂A₂
    eq4 = sympy.Eq(sigma1 * A1, sigma2 * A2)

    # Step 3: Solve the system of four linear equations for the four constants
    solution = sympy.solve([eq1, eq2, eq3, eq4], (A1, B1, A2, B2))

    # Extract the solutions for A1 and A2, which are needed for the E-field
    A1_sol = solution[A1]
    A2_sol = solution[A2]

    # Step 4: Calculate the electric field in each region
    # E_φ = - (1/r) * dΦ/dφ = -A/r
    E1_field = -A1_sol / r
    E2_field = -A2_sol / r

    # Step 5: Print the final expressions for the electric fields
    print("The electric field in Region 1 (0 < phi < pi/2) is:")
    sympy.pprint(sympy.Eq(sympy.Symbol('E_1'), E1_field * sympy.Symbol('i_phi')), use_unicode=False)
    
    print("\nThe electric field in Region 2 (pi/2 < phi < pi) is:")
    sympy.pprint(sympy.Eq(sympy.Symbol('E_2'), E2_field * sympy.Symbol('i_phi')), use_unicode=False)
    
solve_resistor_field()
<<<C>>>