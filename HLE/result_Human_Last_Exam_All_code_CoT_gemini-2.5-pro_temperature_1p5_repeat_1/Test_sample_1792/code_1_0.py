def solve_ordinal_expression():
    """
    This function prints the simplified ordinal expression and its coefficients.
    The derivation is based on the rules of ordinal arithmetic and the Continuum Hypothesis.
    """

    # Representation of the ordinals as strings for printing
    w_2 = "w_2"
    w_1 = "w_1"
    w = "w"

    # The derived coefficients (alpha values)
    alpha_1 = "w_1"
    alpha_2 = "1"
    alpha_3 = "0"
    alpha_4 = "0"

    # Construct the final expression in Cantor Normal Form
    # using f-strings for clear formatting.
    # Parentheses are used around coefficients for readability.
    expression = (f"{w_2} * ({alpha_1}) + "
                  f"{w_1} * ({alpha_2}) + "
                  f"{w} * ({alpha_3}) + "
                  f"{alpha_4}")

    print("The original expression simplifies to the following form:")
    print(expression)
    
    print("\nThe coefficients are:")
    print(f"alpha_1 = {alpha_1}")
    print(f"alpha_2 = {alpha_2}")
    print(f"alpha_3 = {alpha_3}")
    print(f"alpha_4 = {alpha_4}")

solve_ordinal_expression()