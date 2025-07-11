import sympy

def calculate_qhe_conductance():
    """
    This function calculates and prints the four-terminal conductance G_12,34
    for the described quantum Hall device with a quantum point contact (QPC).

    The calculation is based on the Landauer-Büttiker formalism.

    Parameters:
    M: The total number of spin-degenerate edge states.
    N: The number of edge states reflected by the QPC.
    """

    # Define symbolic variables to represent the physical quantities.
    M, N = sympy.symbols('M N')
    
    # The final formula for the conductance G_12,34 is derived as:
    # G_12,34 = (M * (M - N) / N) * (e^2/h)
    # We will print this formula in a clear, readable format.
    
    print("The four-terminal conductance G_12,34 is calculated using the Landauer-Büttiker formalism.")
    print("The derivation considers the voltages at each terminal and the current flow through the device.")
    print("-" * 50)
    print("Given:")
    print("M = total number of edge states")
    print("N = number of reflected edge states at the QPC")
    print("e = elementary charge")
    print("h = Planck's constant")
    print("-" * 50)
    print("The final formula for the conductance G_12,34 is:")
    
    # We construct and print the equation string to show the final result clearly.
    # The prompt asks to "output each number in the final equation!". 
    # Since M and N are variables, we print their symbols in the equation.
    
    final_equation_str = f"G_12,34 = ({M} * ({M} - {N}) / {N}) * (e^2/h)"
    
    print(final_equation_str)

# Execute the function to display the result.
calculate_qhe_conductance()