def calculate_qpc_conductance():
    """
    This function prints the derived formula for the four-terminal conductance
    G_12,34 of the specified quantum Hall device.
    """
    
    # M and N are symbolic variables representing the number of edge states.
    # M: total number of spin-degenerate edge states.
    # N: number of reflected edge states by the QPC.
    
    # The derived formula is G = (e^2/h) * [M * (M - N)] / N.
    # To explicitly show the numbers in the equation, we can expand the numerator:
    # Numerator = M^2 - M*N
    # Denominator = N
    # The equation can be written as G = (e^2/h) * (c1*M^2 + c2*M*N) / (d1*N)
    
    # We define the numerical coefficients.
    c1 = 1
    c2 = -1
    d1 = 1
    
    # The symbol G0 is used for the conductance quantum, e^2/h.
    G0_symbol = "e^2/h"
    
    print("The final equation for the four-terminal conductance G_12,34 is derived as:")
    
    # The f-string substitutes the variables with their numerical values.
    # This format explicitly outputs each number in the final equation.
    print(f"G_12,34 = ({G0_symbol}) * (({c1})*M^2 + ({c2})*M*N) / (({d1})*N)")
    
    print("\nThis expression can be simplified to its more common form:")
    print(f"G_12,34 = ({G0_symbol}) * M*(M - N)/N")

# Execute the function to print the result.
calculate_qpc_conductance()