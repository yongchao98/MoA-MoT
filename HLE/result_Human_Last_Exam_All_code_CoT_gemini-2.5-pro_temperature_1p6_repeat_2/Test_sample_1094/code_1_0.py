import sympy

def calculate_normalized_ac_loss():
    """
    Calculates and prints the symbolic expression for the normalized AC loss
    in a superconducting elliptic bar for a transport current i = Im/Ic < 1.
    """

    # Define the normalized current 'i' as a symbolic variable
    i = sympy.Symbol('i')

    # The formula for normalized loss is L = 2 * [(1-i)ln(1-i) + (1+i)ln(1+i) - i^2]
    # We use sympy.log for the natural logarithm (ln)
    normalized_loss = 2 * ((1 - i) * sympy.log(1 - i) + (1 + i) * sympy.log(1 + i) - i**2)

    # Print the final result
    print("The normalized AC loss per cycle per unit length, 2*pi*Q/(mu_0*Ic^2), is given by the following function of i:")
    
    # Use sympy.pretty_print for a clear mathematical representation
    sympy.pretty_print(normalized_loss)
    
    # The instruction asks to output each number in the final equation.
    # We can represent the equation as a string to show all coefficients and terms clearly.
    # The numbers in the final equation are 2, 1, -1, 1, 1, -1.
    print("\nThe equation in plain text is:")
    print(f"2*((1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i**2)")


if __name__ == '__main__':
    calculate_normalized_ac_loss()