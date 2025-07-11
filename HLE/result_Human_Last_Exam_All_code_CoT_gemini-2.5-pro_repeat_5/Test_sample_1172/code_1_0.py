import math

def print_final_expression():
    """
    This function prints the derived expression for the change in mutual inductance.
    """
    
    # Define the components of the formula as strings
    numerator = "mu_0 * h**2"
    denominator = "2 * pi * R1**2"
    
    # The final expression
    final_expression = f"Delta_M = ({numerator}) / ({denominator})"
    
    print("The final expression for the change in mutual inductance per unit length (Delta_M) is:")
    print(final_expression)
    print("\nWhere:")
    print("  - Delta_M: The change in mutual inductance (M2 - M1)")
    print("  - mu_0: The permeability of free space (a physical constant)")
    print("  - h: The separation distance between the wires in each circuit")
    print("  - pi: The mathematical constant pi")
    print("  - R1: The inner radius of the concentrator shell")
    
    # As requested, explicitly output the number(s) in the final equation
    print("\nThe final equation contains the following number:")
    print("  - The number 2 in the denominator of the expression.")

# Execute the function to display the result
print_final_expression()