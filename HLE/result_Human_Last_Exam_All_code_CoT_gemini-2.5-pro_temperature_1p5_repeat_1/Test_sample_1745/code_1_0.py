import sympy as sp

def solve_edl_potential():
    """
    This function derives the expression for the Electrical Double-Layer (EDL)
    potential distribution psi(y) in a parallel-plate microchannel using sympy.

    The governing equation is the linearized Poisson-Boltzmann equation:
    d^2(psi)/dy^2 - k^2 * psi = 0

    With the following boundary conditions:
    1. At y=0 (bottom plate): psi(0) = zeta_a1 = zeta_1 * (1 + beta * k)
    2. At y=H (top plate): psi(H) = 0 (since zeta_2 = 0)
    """

    # 1. Define the symbolic variables
    y, H, k, beta, zeta_1 = sp.symbols('y H k beta zeta_1', real=True, positive=True)
    psi = sp.Function('psi')

    # 2. Define the slip-dependent zeta potential at the bottom wall
    zeta_a1 = zeta_1 * (1 + beta * k)

    # 3. Define the Ordinary Differential Equation (ODE)
    ode = sp.Eq(psi(y).diff(y, 2) - k**2 * psi(y), 0)

    # 4. Define the boundary conditions
    # In sympy dsolve, these are provided as a dictionary of initial conditions (ics)
    ics = {psi(0): zeta_a1, psi(H): 0}

    # 5. Solve the ODE with the given boundary conditions
    # The general solution is C1*exp(k*y) + C2*exp(-k*y)
    # dsolve will find C1 and C2 based on the ics
    solution = sp.dsolve(ode, psi(y), ics=ics)

    # The solution is returned as an Equality, e.g., Eq(psi(y), expression)
    # We extract the expression for psi(y) (the right-hand side)
    psi_expression = solution.rhs

    # 6. Simplify the expression to the standard hyperbolic form
    # sympy.simplify() is powerful and can convert the exponential form to hyperbolic
    simplified_psi = sp.simplify(psi_expression)

    # 7. Print the final result in a readable format
    print("Derivation of the EDL potential distribution psi(y):")
    print("-" * 50)
    print("Governing Equation: d^2(psi)/dy^2 - k^2*psi = 0")
    print(f"Boundary Condition 1 (y=0): psi(0) = {zeta_a1}")
    print("Boundary Condition 2 (y=H): psi(H) = 0")
    print("-" * 50)
    print("The final expression for the potential distribution is:")
    
    # Pretty print the final equation
    final_equation = sp.Eq(psi(y), simplified_psi)
    sp.pprint(final_equation, use_unicode=True)
    
    # Print the equation in a single line format as well for clarity
    print("\nIn single-line format:")
    print(f"psi(y) = {simplified_psi}")


# Execute the function to find and print the expression
solve_edl_potential()

# For the final answer format, we extract the mathematical formula.
# We'll represent the simplified expression here as requested.
# The expression is: (zeta_1*(1 + beta*k) * sinh(k*(H - y))) / sinh(k*H)