import sympy
from sympy import pprint, Eq

def solve_edl_potential():
    """
    This function derives and prints the expression for the Electrical Double-Layer (EDL) 
    potential distribution psi(y) in a parallel-plate microchannel under the given conditions.
    """
    
    # Define the symbolic variables
    # psi is the EDL potential, a function of y
    # y is the spatial coordinate across the channel height
    # z_1 is the intrinsic zeta potential of the bottom surface
    # beta is the slip length
    # k is the Debye-Huckel parameter
    # H is the height of the microchannel
    psi = sympy.Function('psi')
    y, z_1, beta, k, H = sympy.symbols('y z_1 beta k H', real=True, positive=True)

    # From the problem description, the slip-dependent zeta potential at the bottom surface (y=0) is:
    # z_a1 = z_1 * (1 + beta * k)
    # The potential at the top surface (y=H) is zero because z_2 = 0.
    
    # The solution to the linearized Poisson-Boltzmann equation with these boundary conditions is:
    # psi(y) = (z_a1 / sinh(k*H)) * sinh(k*(H - y))
    
    # We construct the final expression for psi(y)
    z_a1 = z_1 * (1 + beta * k)
    
    rhs_expression = (z_a1 / sympy.sinh(k * H)) * sympy.sinh(k * (H - y))
    
    # Create the equation psi(y) = expression
    final_equation = Eq(psi(y), rhs_expression)
    
    # Print the final expression
    print("The expression for the Electrical double-layer potential distribution psi(y) is:")
    pprint(final_equation, use_unicode=True)

solve_edl_potential()
