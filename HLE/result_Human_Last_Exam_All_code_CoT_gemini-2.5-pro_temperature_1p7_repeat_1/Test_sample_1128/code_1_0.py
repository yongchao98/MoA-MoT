import math
import sys

def count_cool_strings():
    """
    Calculates the number of "cool strings" of maximal length for n symbols.
    """
    try:
        # We read n from standard input.
        n_input = input("Enter the number of symbols (n): ")
        n = int(n_input)
        if n < 1:
            print("Error: n must be a positive integer.", file=sys.stderr)
            return
    except (ValueError, EOFError):
        # If input is not provided or invalid, we use a default value for demonstration.
        n = 4
        print(f"Invalid input. Using default value n = {n}")
        
    # The number of cool strings of maximal length 3n is given by the formula (2n)! / (n+1)!.
    # This formula simplifies to the product of integers from n + 2 to 2n.

    # For n=1, the product range is empty, so the result is 1.
    if n == 1:
        result = 1
        equation_str = "1"
    else:
        # Calculate the result by multiplying the numbers in the specified range.
        result = 1
        # The range is [n+2, 2n]. The upper bound for Python's range() is exclusive.
        for i in range(n + 2, 2 * n + 1):
            result *= i
        
        # Build the equation string for printing.
        product_terms = [str(i) for i in range(n + 2, 2 * n + 1)]
        equation_str = " * ".join(product_terms)

    print(f"\nFor n = {n}, the number of cool strings is calculated as:")
    # As requested, printing each number in the final equation.
    print(f"{equation_str} = {result}")

if __name__ == '__main__':
    count_cool_strings()