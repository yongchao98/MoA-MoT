import sympy as sp

def solve_edl_potential():
    """
    This script solves the linearized Poisson-Boltzmann equation for the
    Electric Double-Layer (EDL) potential in a parallel-plate microchannel.

    The problem considers a case with asymmetric zeta potentials, where
    the potential at the top wall is zero.
    """

    # 1. Define the symbolic variables
    # y: spatial coordinate across the channel height
    # H: height of the microchannel
    # k: Debye-Hückel parameter
    # z_1: Base zeta potential at the bottom wall (y=0)
    # beta: Slip coefficient
    y, H, k, z_1, beta = sp.symbols('y H k z_1 beta')

    # psi is the EDL potential, a function of y
    psi = sp.Function('psi')

    # 2. Formulate the governing differential equation
    # d^2(psi(y))/dy^2 - k^2 * psi(y) = 0
    ode = sp.Eq(psi(y).diff(y, y) - k**2 * psi(y), 0)

    # 3. Define the boundary conditions
    # At the bottom wall (y=0), the zeta potential is z_{a1} = z_1(1 + beta*k)
    # At the top wall (y=H), the zeta potential z_2 is 0.
    zeta_1_slip = z_1 * (1 + beta * k)
    
    # Dictionary for the initial/boundary conditions for dsolve
    ics = {
        psi(0): zeta_1_slip,
        psi(H): 0
    }

    # 4. Solve the ODE with the given boundary conditions
    # sympy's dsolve handles the boundary value problem automatically.
    solution_eq = sp.dsolve(ode, psi(y), ics=ics)

    # The result from dsolve is an equation, psi(y) = expression.
    # We extract the expression part (the right-hand side).
    potential_expression = solution_eq.rhs

    # 5. Print the final expression for the potential distribution.
    # The output is formatted for clarity.
    # The final print statement shows all the variables in the equation as requested.
    print("The governing ODE is: d^2(psi)/dy^2 - k^2*psi = 0")
    print(f"Boundary condition at y=0: psi(0) = {zeta_1_slip}")
    print("Boundary condition at y=H: psi(H) = 0")
    print("\n" + "="*50)
    print("The final expression for the Electrical Double-Layer potential distribution is:")
    # We create a final equation psi(y) = ... to print it nicely
    final_equation = sp.Eq(psi(y), potential_expression)
    print(final_equation)
    
    # Printing the LaTeX form for the final answer block.
    # print(sp.latex(potential_expression))
    
if __name__ == '__main__':
    solve_edl_potential()
