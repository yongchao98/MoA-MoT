import sympy

def display_partition_function():
    """
    This function displays the derived partition function Z and its components.
    """

    # Define the symbols used in the equation
    Z, beta, mu = sympy.symbols('Z β μ')
    
    # Define the components of the final equation
    # The equation is Z = 1 / (1 - e^(βμ))
    numerator = 1
    denominator_term_1 = 1
    
    # Explain the equation
    print("The derived partition function Z for a single bosonic state is given by the formula:")
    print("Z = numerator / (term1 - term2)\n")
    
    print("The components of the equation are:")
    # Print each number/component of the final equation
    print(f"The numerator is: {numerator}")
    print(f"The first term in the denominator is: {denominator_term_1}")
    print("The second term in the denominator is e^(βμ)\n")
    
    # Construct the final equation using the components
    final_equation_str = f"{Z} = {numerator} / ({denominator_term_1} - exp({beta}{mu}))"
    
    # Print the final, formatted equation
    print("The final equation is:")
    print(final_equation_str)
    
    print("\nWhere:")
    print("Z: The Grand Canonical Partition Function")
    print("β: Inverse temperature (1 / (k_B * T))")
    print("μ: Chemical Potential (must be negative for this result)")
    print("exp: The exponential function")

if __name__ == "__main__":
    display_partition_function()
