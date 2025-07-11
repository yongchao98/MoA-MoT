import sympy as sp

def solve_potential_distribution():
    """
    This function derives and prints the expression for the Electrical Double-Layer (EDL)
    potential distribution psi(y) in a parallel-plate microchannel under the specified conditions.
    """
    # 1. Define the symbolic variables from the problem statement.
    # y: the spatial coordinate across the channel height, with the origin at the center.
    # H: the total height of the microchannel.
    # z_1: the base zeta potential at the bottom surface (j=1). The prompt uses 'z' for zeta potential.
    # beta: the slip length parameter.
    # k: the Debye-Huckel parameter.
    y, H, z_1, beta, k = sp.symbols('y H z_1 beta k')

    # The potential psi is a function of y.
    psi = sp.Function('psi')(y)

    # 2. Define the boundary conditions based on the problem.
    # The potential at the bottom wall (y = -H/2) is the slip-dependent zeta potential z_a1.
    z_a1 = z_1 * (1 + beta * k)
    # The potential at the top wall (y = H/2) is 0 since z_2 = 0.

    # 3. Construct the final expression.
    # The expression is the solution to the linearized Poisson-Boltzmann equation (d^2(psi)/dy^2 = k^2*psi)
    # with the specified boundary conditions.
    # The mathematical derivation yields:
    # psi(y) = z_a1 * sinh(k*(H/2 - y)) / sinh(k*H)
    
    numerator = sp.sinh(k * (H/2 - y))
    denominator = sp.sinh(k * H)
    potential_distribution = z_a1 * numerator / denominator

    # 4. Create a SymPy equation object for clear representation.
    final_equation = sp.Eq(psi, potential_distribution)

    # 5. Print the final result as an equation.
    print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
    
    # We print the string representation of the SymPy equation.
    # This ensures all numbers (1 from the zeta potential expression, 2 from H/2) are included.
    print(str(final_equation))

solve_potential_distribution()