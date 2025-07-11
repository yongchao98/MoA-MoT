import sympy

def calculate_amplitude_at_L():
    """
    Calculates and displays the symbolic formula for the electric field amplitude
    at the end of the slab.
    """
    # Define the symbols for the variables in the problem
    A = sympy.Symbol('A')      # Initial amplitude of the electric wave
    alpha = sympy.Symbol('alpha') # Rate of change of the slab's properties
    L = sympy.Symbol('L')      # Length of the slab
    c = sympy.Symbol('c')      # Speed of light in vacuum
    A_L = sympy.Symbol('A(L)') # Amplitude at the end of the slab, as a function of L

    # The derived formula for the amplitude at x=L
    final_amplitude_expr = A * sympy.exp(-alpha * L / c)

    # Create and print the equation
    # The final equation shows how the amplitude at L, A(L), depends on the initial
    # amplitude A, the slab length L, the material property alpha, and the speed of light c.
    equation = sympy.Eq(A_L, final_amplitude_expr)

    print("The final equation for the amplitude of the electric field at x=L is:")
    
    # Use sympy.pprint for a more readable mathematical output
    sympy.pprint(equation, use_unicode=True)

if __name__ == '__main__':
    calculate_amplitude_at_L()