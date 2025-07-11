import sympy

def derive_force_law():
    """
    This function uses the sympy library to derive and display the force law
    for a thermally isolated, freely jointed polymer chain.

    The derivation follows these steps:
    1. Define the symbolic variables: x, l, n, E_0.
    2. Define the energy E(x) as a function of extension x, based on thermodynamic principles.
       E(x) = E_0 * exp((3 * x**2) / (2 * n**2 * l**2))
    3. Calculate the force F by differentiating the energy with respect to x.
       F = -dE/dx
    4. Print the resulting symbolic expression for the force.
    """
    
    # Define the symbols for the variables in the problem
    x, l, n, E_0 = sympy.symbols('x l n E_0', real=True, positive=True)

    # Define the expression for the energy E(x) as a function of the end-to-end distance x
    # E_0 is the kinetic energy at x=0
    energy_expression = E_0 * sympy.exp((3 * x**2) / (2 * n**2 * l**2))

    # Calculate the force F by taking the negative derivative of the energy with respect to x
    force_expression = -sympy.diff(energy_expression, x)

    # Print the final force law
    # The output shows the attractive force exerted by the polymer ends on each other.
    print("The derived force law F(x) is:")
    # Using str() provides a clean, single-line output
    print(str(force_expression))

# Execute the function to see the result
derive_force_law()