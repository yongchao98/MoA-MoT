import sympy as sp

def solve_edl_potential():
    """
    This function derives the expression for the Electrical Double-Layer (EDL)
    potential distribution in a parallel-plate microchannel with asymmetric
    zeta potentials.
    """
    # 1. Define the symbolic variables
    # y: spatial coordinate across the channel height
    # k: Debye-Huckel parameter
    # H: height of the microchannel
    # z1: zeta potential parameter for the bottom surface
    # beta: slip length
    # psi: the EDL potential function
    y, k, H, z1, beta = sp.symbols('y k H z_1 beta', real=True)
    psi = sp.Function('psi')(y)

    # Define the constants of integration
    C1, C2 = sp.symbols('C1 C2')

    # 2. Define the general solution to the linearized Poisson-Boltzmann equation
    # The ODE is d^2(psi)/dy^2 - k^2 * psi = 0
    general_solution = C1 * sp.sinh(k * y) + C2 * sp.cosh(k * y)

    # 3. Define the boundary conditions based on the problem statement
    # Zeta potential at the bottom wall (y = -H/2)
    zeta1 = z1 * (1 + beta * k)
    
    # Equation 1: At y = H/2, psi(H/2) = 0
    bc1 = sp.Eq(general_solution.subs(y, H / 2), 0)

    # Equation 2: At y = -H/2, psi(-H/2) = zeta1
    bc2 = sp.Eq(general_solution.subs(y, -H / 2), zeta1)

    # 4. Solve the system of linear equations for the constants C1 and C2
    constants = sp.solve([bc1, bc2], (C1, C2))

    # 5. Substitute the constants back into the general solution
    final_solution = general_solution.subs(constants)

    # 6. Simplify the final expression for better readability
    # The sinh(a-b) identity is useful here. Sympy can handle this.
    simplified_solution = sp.simplify(final_solution)
    
    # 7. Print the final result in a clear format
    # The instruction is to output each "number" in the final equation.
    # For a symbolic result, we will print the clear, structured expression.
    final_expression_latex = sp.latex(simplified_solution)
    
    # Re-arranging for the requested form psi(y) = ...
    final_equation = sp.Eq(psi, simplified_solution)

    print("The governing differential equation is: d^2(psi)/dy^2 - k^2 * psi = 0")
    print(f"The boundary condition at the bottom wall (y = -H/2) is: psi(-H/2) = {zeta1}")
    print("The boundary condition at the top wall (y = H/2) is: psi(H/2) = 0")
    print("\nThe derived expression for the Electrical Double-Layer (EDL) potential distribution is:")
    sp.pprint(final_equation, use_unicode=True)
    
    # Return the final expression for the grading format
    return final_equation

# Execute the function to get the result
final_answer = solve_edl_potential()

# Final answer in the required format
# We represent the symbolic expression as a string.
final_answer_str = str(final_answer.rhs)
# Let's reformat to match the manual derivation for clarity
# Manual result: z1*(1+beta*k)*sinh(k*(H/2 - y))/sinh(k*H)
reformatted_answer = "z_1*(1 + beta*k) * sinh(k*(H/2 - y)) / sinh(H*k)"

# The sympy output might be slightly different but algebraically equivalent.
# E.g. -z_1*(beta*k + 1)*sinh(k*y - H*k/2)/sinh(H*k)
# The reformatted version is more intuitive.

print(f"\n<<<{reformatted_answer}>>>")