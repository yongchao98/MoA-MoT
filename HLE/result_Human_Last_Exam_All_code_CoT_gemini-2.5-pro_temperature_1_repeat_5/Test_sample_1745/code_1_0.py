import sympy

def electrical_potential_distribution():
    """
    This function symbolically derives and prints the expression for the 
    Electrical Double-Layer (EDL) potential distribution, psi(y),
    in an asymmetric microchannel.
    """
    
    # 1. Define the symbolic variables required for the equation.
    # y: the spatial coordinate across the channel height
    # z_1: the base zeta potential value at the bottom surface
    # beta: the slip length parameter
    # k: the Debye-Huckel parameter (inverse of EDL thickness)
    # H: the height of the microchannel
    y, z_1, beta, k, H = sympy.symbols('y z_1 beta k H', real=True)
    
    # Define psi as a symbolic function of y
    psi = sympy.Function('psi')

    # 2. Define the expression for the zeta potential at the bottom surface (y=-H/2).
    # The problem gives the slip-dependant zeta potential as z_j(1 + beta*k).
    # For the bottom surface (j=1), this is zeta_1.
    zeta_1 = z_1 * (1 + beta * k)
    
    # 3. Formulate the final expression for the potential distribution psi(y).
    # This expression is the solution to the linearized Poisson-Boltzmann equation
    # d^2(psi)/dy^2 = k^2*psi with boundary conditions:
    # psi(y=H/2) = 0 and psi(y=-H/2) = zeta_1.
    # The solution is: psi(y) = zeta_1 * sinh(k*(H/2 - y)) / sinh(k*H)
    
    # Construct the numerator and denominator of the solution
    numerator = sympy.sinh(k * (H / 2 - y))
    denominator = sympy.sinh(k * H)
    
    # Combine the parts to form the right-hand side of the equation
    rhs = zeta_1 * numerator / denominator
    
    # 4. Create and print the final equation psi(y) = ...
    # The sympy.Eq function creates a symbolic equality.
    final_equation = sympy.Eq(psi(y), rhs)
    
    print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
    
    # The print function will output the formatted equation. 
    # This fulfills the requirement to "output each number in the final equation" 
    # by showing all its symbolic components.
    print(final_equation)

# Execute the function to display the result
electrical_potential_distribution()