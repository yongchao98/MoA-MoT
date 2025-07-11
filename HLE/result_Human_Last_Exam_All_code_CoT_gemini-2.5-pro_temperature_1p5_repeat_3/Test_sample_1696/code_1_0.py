def solve():
    """
    This function prints the closed-form formula for P(n).
    """
    n_str = "n"
    num_part_1 = "(2 * pi)"
    num_part_2 = f"({n_str} * ({n_str} - 1))"
    num_part_3 = "4"
    den_part = f"sqrt({n_str}!)"

    print("The closed-form formula for P(n) is:")
    
    final_equation = f"P(n) = ({num_part_1} ^ ({num_part_2} / {num_part_3})) / {den_part}"
    
    print(final_equation)

    # To fulfill the requirement to output each number in the equation.
    # Note that the formula involves the constant pi and the variable n.
    # The numbers in the formula are 2, 1, 4.
    print("\nBreaking down the formula:")
    print("P(n) = ( (2 * pi) ^ ( (n * (n - 1)) / 4 ) ) / sqrt(n!)")
    print("The numbers in the equation are:")
    print("In the base of the exponent: 2")
    print("In the exponent's numerator: 1")
    print("In the exponent's denominator: 4")

solve()