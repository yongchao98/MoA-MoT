def display_exact_formula():
    """
    This function prints the derived exact symbolic formula for l_k(n)
    as requested by the user. The derivation steps are outlined above.
    The function prints out the final equation, showing each numerical constant.
    """

    # These are the numerical coefficients and constants in the derived formula:
    # l_k(n) = (c1/c2)*ln(n + c3) - k^2 * (c4 - c5/n)
    c1 = 1
    c2 = 2
    c3 = 1
    c4 = 2
    c5 = 1

    # We construct and print the formula string.
    formula_string = f"l_k(n) = ({c1}/{c2})*ln(n + {c3}) - k^2 * ({c4} - {c5}/n)"

    print("The exact value of the function l_k(n) is given by the formula:")
    print(formula_string)


if __name__ == "__main__":
    display_exact_formula()