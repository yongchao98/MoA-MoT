import sympy
from sympy import symbols, Function, Eq, cosh, sinh, solve, simplify

def solve_edl_potential():
    """
    This function derives the expression for the Electrical Double-Layer (EDL)
    potential distribution psi(y) in a parallel-plate microchannel
    using the specified boundary conditions.
    """
    # 1. Define the symbolic variables
    # y: spatial coordinate across the channel height
    # k: Debye-Hückel parameter
    # H: height of the microchannel
    # z1: zeta potential parameter for the bottom surface
    # beta: slip length
    # C1, C2: constants of integration
    y, H, k, z1, beta = symbols('y H k z_1 beta')
    C1, C2 = symbols('C1 C2')
    
    # Let's assume the physical parameters are positive real numbers
    H, k, z1, beta = symbols('H k z_1 beta', real=True, positive=True)

    # 2. Define the general solution for psi(y)
    # The general solution to d^2(psi)/dy^2 = k^2*psi is psi(y) = C1*cosh(k*y) + C2*sinh(k*y)
    psi_y = C1 * cosh(k * y) + C2 * sinh(k * y)
    
    # 3. Define the boundary conditions
    # BC1: at y = H/2, psi = 0
    # BC2: at y = -H/2, psi = z1 * (1 + beta * k)
    
    # Substitute y in the general solution to create equations from BCs
    eq1 = Eq(psi_y.subs(y, H / 2), 0)
    eq2 = Eq(psi_y.subs(y, -H / 2), z1 * (1 + beta * k))

    # 4. Solve the system of linear equations for the constants C1 and C2
    constants = solve([eq1, eq2], (C1, C2))
    
    # The 'solve' function returns a dictionary of the solutions
    # e.g., {C1: expression_for_C1, C2: expression_for_C2}

    # 5. Substitute the constants back into the general solution
    final_psi = psi_y.subs(constants)

    # 6. Simplify the final expression
    # sympy can simplify the expression using hyperbolic identities.
    simplified_psi = simplify(final_psi)
    
    # 7. Print the final result
    print("The derived expression for the Electrical double-layer potential distribution is:")
    
    # To make the output clear, we create and print a Sympy equation psi(y) = result
    psi_func = Function('psi')
    final_equation = Eq(psi_func(y), simplified_psi)
    print(final_equation)

# Execute the function to find and print the expression
solve_edl_potential()

<<<z_1*(1 + beta*k)*sinh(k*(H/2 - y))/sinh(k*H)>>>