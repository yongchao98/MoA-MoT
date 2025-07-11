import sympy

def calculate_final_amplitude():
    """
    This function calculates and prints the symbolic amplitude of an electromagnetic wave
    after passing through a slab with time-varying properties.
    The derivation leads to the final formula which is constructed and printed here.
    """

    # Define the symbolic variables used in the problem
    A = sympy.Symbol('A')          # Incident amplitude
    alpha = sympy.Symbol('alpha')  # Time-variation coefficient of the material
    L = sympy.Symbol('L')          # Length of the slab
    c = sympy.Symbol('c')          # Speed of light in vacuum

    # From the physics derivation, the initial amplitude A is amplified (or attenuated)
    # by a factor of exp(alpha * L / c) as it propagates through the slab.
    final_amplitude = A * sympy.exp(alpha * L / c)

    # Print the final result in a clear, human-readable format,
    # showing each component of the final equation as requested.
    print("The amplitude of the electric field at the rightmost boundary (x=L) is given by the formula:")
    print(f"Final Amplitude = {A} * exp(({alpha} * {L}) / {c})")

# Execute the function to display the result
calculate_final_amplitude()