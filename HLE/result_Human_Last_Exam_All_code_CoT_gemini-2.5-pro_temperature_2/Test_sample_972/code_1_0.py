import sympy

def calculate_output_amplitude():
    """
    Calculates and prints the symbolic expression for the electric field amplitude
    at the rightmost boundary of the slab.

    The derivation shows that the output amplitude A_out is related to the
    input amplitude A_in by the formula:
    A_out = A_in * exp(-alpha * L / c)
    where:
    A_in is the initial amplitude of the wave.
    alpha is the rate of change of the slab's properties.
    L is the length of the slab.
    c is the speed of light in vacuum.
    """

    # Define symbolic variables
    A_in = sympy.Symbol('A')
    alpha = sympy.Symbol('alpha')
    L = sympy.Symbol('L')
    c = sympy.Symbol('c')

    # Define the expression for the output amplitude
    A_out = A_in * sympy.exp(-alpha * L / c)

    # Print the equation
    # We create the symbol A_out for prettier printing
    A_out_symbol = sympy.Symbol("A_out(L)")
    final_equation = sympy.Eq(A_out_symbol, A_out)
    
    print("The amplitude of the electric field at the rightmost boundary (x=L) is:")
    sympy.pprint(final_equation, use_unicode=True)

    # Print the equation in a linear format for clarity
    print("\nIn linear format:")
    print(f"{A_out_symbol} = {A_in} * exp(-({alpha}*{L})/{c})")
    
if __name__ == '__main__':
    calculate_output_amplitude()