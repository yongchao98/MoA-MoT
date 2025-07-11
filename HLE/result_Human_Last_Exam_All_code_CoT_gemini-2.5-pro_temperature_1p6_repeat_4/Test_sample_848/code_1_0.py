import math

def solve_limit():
    """
    This function calculates the integer part of the specified limit.
    The problem reduces to finding solutions to a Diophantine equation,
    which fall into two families. The number of solutions up to N
    can be estimated using the growth rates of these solution sequences.
    """

    # The problem asks for the integer part of 10^4 * L, where L is the limit.
    # The limit L is given by the expression: 2/ln(alpha_13) + 2/ln(alpha_5)

    # First root is from the characteristic equation r^2 - 13r + 1 = 0, for the m=13 family.
    alpha_13 = (13 + math.sqrt(165)) / 2

    # Second root is from r^2 - 5r + 1 = 0, for the m=5 family.
    alpha_5 = (5 + math.sqrt(21)) / 2

    # Calculate the natural logarithms of these roots.
    ln_alpha_13 = math.log(alpha_13)
    ln_alpha_5 = math.log(alpha_5)

    # As requested, we output the numbers used in the final calculation.
    # The equation is: 10000 * (2/ln_alpha_13 + 2/ln_alpha_5)
    print("The final calculation is based on the following values:")
    print(f"  Factor: 10000")
    print(f"  Constant for the first family: 2")
    print(f"  Logarithm of the first root (ln(alpha_13)): {ln_alpha_13}")
    print(f"  Constant for the second family: 2")
    print(f"  Logarithm of the second root (ln(alpha_5)): {ln_alpha_5}")

    # Calculate the limit value
    limit_value = (2 / ln_alpha_13) + (2 / ln_alpha_5)

    # Multiply by 10^4
    final_result = 10000 * limit_value

    # The final answer is the integer part of this result.
    print("\nThe integer part of the final result is:")
    print(int(final_result))

solve_limit()
