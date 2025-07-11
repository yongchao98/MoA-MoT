import sympy as sp

def solve_polymer_force():
    """
    This function derives the force law for a thermally isolated,
    freely jointed polymer chain using symbolic mathematics.
    """
    # Step 1: Define the symbolic variables for the problem.
    # n: number of segments
    # l: length of each segment
    # x: end-to-end separation of the polymer
    # E_0: kinetic energy at zero extension (x=0)
    # k_B: Boltzmann's constant
    # E_x: kinetic energy at extension x
    n, l, x, E_0, k_B = sp.symbols('n l x E_0 k_B', positive=True, real=True)
    E_x = sp.Symbol('E_x', positive=True, real=True)

    # Define the configurational and kinetic entropy components (ignoring additive constants).
    # S_config is based on the Gaussian distribution for the end-to-end distance.
    S_config = - (3 * k_B * x**2) / (2 * n * l**2)
    # S_kinetic is based on the phase space volume for 3n kinetic degrees of freedom.
    S_kinetic = (3 * n * k_B / 2) * sp.log(E_x)
    
    # The total entropy S(E, x) is the sum of the two.
    S_total = S_kinetic + S_config

    # Step 2: Apply the adiabatic condition (constant entropy).
    # The entropy S(E(x), x) is equal to the initial entropy S(E(0), 0).
    # S_initial is the total entropy when x=0 and the energy is E_0.
    S_initial = S_total.subs({E_x: E_0, x: 0})
    
    # We form an equation based on the conservation of entropy.
    entropy_conservation_eq = sp.Eq(S_total, S_initial)
    
    # Solve this equation for E_x to find how energy changes with extension.
    E_x_solution = sp.solve(entropy_conservation_eq, E_x)[0]

    # Step 3: Calculate the force f.
    # The force is given by f = -(∂S/∂x)_E / (∂S/∂E)_x.
    # First, calculate the partial derivatives of the total entropy.
    dS_dx = sp.diff(S_total, x)
    dS_dE = sp.diff(S_total, E_x)

    # Now, calculate the force in terms of the energy E_x at extension x.
    force_intermediate = -dS_dx / dS_dE

    # Step 4: Express the force in terms of the initial energy E(0).
    # Substitute the expression for E_x (from E_x_solution) into the force equation.
    final_force_law = force_intermediate.subs(E_x, E_x_solution)
    
    # Simplify the final expression.
    final_force_law = sp.simplify(final_force_law)

    # Print the final result in a readable format.
    print("The force law f(x) for the thermally isolated polymer is:")
    # The final equation shows all variables and numeric constants.
    f_x = sp.Symbol('f(x)')
    display_equation = sp.Eq(f_x, final_force_law)
    sp.pprint(display_equation, use_unicode=True)

# Execute the function to find and print the solution.
solve_polymer_force()