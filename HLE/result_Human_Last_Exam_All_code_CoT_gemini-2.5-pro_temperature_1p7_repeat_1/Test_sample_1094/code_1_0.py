import sympy

def display_normalized_loss_formula():
    """
    This function uses the sympy library to symbolically represent and print
    the normalized AC loss for a superconductor in the low-current limit (i < 1).

    The general formula for AC loss per cycle per unit length (Q) is:
    Q = (mu_0 * Ic**2 / pi) * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i**2]
    where i = Im/Ic.

    This script derives and prints the normalized form: 2*pi*Q / (mu_0 * Ic**2).
    """
    # Define 'i' as a symbolic variable, representing Im/Ic
    # We add constraints to reflect the physical reality of the problem
    i = sympy.Symbol('i', real=True, positive=True)

    # The core part of the loss formula, let's call it f(i)
    # f(i) = (1-i)*ln(1-i) + (1+i)*ln(1+i) - i**2
    # Using sympy.log for the natural logarithm (ln)
    f_i = (1 - i) * sympy.log(1 - i) + (1 + i) * sympy.log(1 + i) - i**2

    # The full expression for normalized loss is 2 * f(i)
    normalized_loss = 2 * f_i

    # Print the resulting equation for the user
    print("The normalized AC loss, 2*pi*Q/(mu_0*Ic^2), as a function of i=Im/Ic is:")
    # Using sympy.pretty_print for a more readable output if available,
    # otherwise default to standard print.
    try:
        sympy.pretty_print(normalized_loss)
    except ImportError:
        print(normalized_loss)

    # Let's also print the equation with all the numbers in a single line as requested
    print("\nFinal equation in a single line:")
    # The expansion makes sure all terms and coefficients are clearly visible.
    final_equation = sympy.expand(normalized_loss, force=True)
    # The numbers in the final equation are: 2, -2, 1, -1
    print(final_equation)


if __name__ == '__main__':
    display_normalized_loss_formula()
