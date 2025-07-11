def print_closed_form():
    """
    Prints the coefficients of the closed-form expression for the sum.
    The closed form is of the structure: 4^n * P(n) / D
    where P(n) is a polynomial in n and D is an integer denominator.
    """

    # The polynomial is P(n) = 63n^5 + 245n^4 + 355n^3 + 235n^2 + 70n + 8
    coeffs = [63, 245, 355, 235, 70, 8]
    denominator = 8

    # Format the equation string
    poly_str = " + ".join([f"{coeffs[i]}*n^{5-i}" for i in range(len(coeffs)-1)])
    poly_str += f" + {coeffs[-1]}"

    equation = f"S_n = 4^n * ({poly_str}) / {denominator}"

    print("The closed form for the sum is:")
    print(equation)
    print("\nHere are the numbers in the final equation:")
    print("Coefficients of the numerator polynomial P(n) from n^5 to n^0:")
    print(coeffs)
    print("Denominator D:")
    print(denominator)

print_closed_form()