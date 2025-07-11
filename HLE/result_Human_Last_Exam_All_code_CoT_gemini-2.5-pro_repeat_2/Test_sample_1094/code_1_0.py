import sympy as sp

def calculate_normalized_loss_expression():
    """
    Calculates and prints the symbolic expression for the normalized AC loss
    in a superconductor with i = Im/Ic < 1.
    """
    # Define the symbol for the normalized current i = Im/Ic
    i = sp.Symbol('i')

    # The normalized loss per cycle is given by the expression:
    # 2 * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2]
    # Note: sp.log is the natural logarithm (ln)
    # The expression contains the numbers 2, 1, and -1 (from -i**2) as requested.
    loss_expression = 2 * ((1 - i) * sp.log(1 - i) + (1 + i) * sp.log(1 + i) - i**2)

    # Print the final result in a clear format
    # The left-hand side represents the quantity we want to find.
    # The right-hand side is the derived function of i.
    print("The normalized AC loss per cycle, 2*pi*Q/(mu_0*Ic^2), as a function of i is:")
    
    # Using sympy's pretty print for a more readable mathematical format
    # sp.pprint() prints directly to the console.
    # We can also use a string representation for more control.
    lhs = sp.sympify("2*pi*Q/(mu_0*Ic**2)")
    eq = sp.Eq(lhs, loss_expression)
    
    # Print the equation using standard print on the string representation
    print(str(eq))

if __name__ == '__main__':
    calculate_normalized_loss_expression()