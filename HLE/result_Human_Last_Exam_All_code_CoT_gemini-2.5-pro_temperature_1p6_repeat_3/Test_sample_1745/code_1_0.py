import math

def print_potential_distribution_expression():
    """
    This function prints the derived expression for the Electrical double-layer
    potential distribution psi(y).

    The expression is derived from the linearized Poisson-Boltzmann equation
    d^2(psi)/dy^2 = k^2 * psi, with the boundary conditions:
    1. At the bottom plate (y=0): psi(0) = z1*(1 + beta*k)
    2. At the top plate (y=H): psi(H) = 0

    The symbols represent:
    psi(y): Electrical potential at position y
    z1: Zeta potential at the bottom surface (y=0) without slip
    beta: Slip length
    k: Debye-Huckel parameter
    H: Height of the microchannel
    y: Vertical position in the channel (from 0 to H)
    sinh: Hyperbolic sine function
    """

    # The derived mathematical expression is:
    # psi(y) = z1 * (1 + beta * k) * sinh(k * (H - y)) / sinh(k * H)
    
    # We will print this expression as a string, clearly showing all terms
    # as requested, including the number 1.
    
    expression_string = "psi(y) = z1*(1 + beta*k)*sinh(k*(H - y))/sinh(k*H)"
    
    print("The final expression for the Electrical double-layer potential distribution is:")
    print(expression_string)

# Execute the function to print the result
print_potential_distribution_expression()