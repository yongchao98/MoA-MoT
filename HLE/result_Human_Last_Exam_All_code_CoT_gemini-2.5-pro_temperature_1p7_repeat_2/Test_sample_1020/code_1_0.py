import math
import sys

def calculate_complexity(n, d):
    """
    Calculates the complexity for computing the polynomial f_n,d based on a
    standard construction for non-commutative polynomials.

    Args:
      n: The number of available indices.
      d: The degree of the polynomial.

    Returns:
      The complexity value.
    """
    sum_of_widths = 0
    for i in range(1, d):
        try:
            term = math.comb(n, i)
        except ValueError:
            # This happens if n < i, which shouldn't occur given the problem constraints
            term = 0
        sum_of_widths += term

    complexity = 2 + sum_of_widths
    return complexity

def main():
    """
    Main function to get user input and print the complexity calculation.
    """
    # You can change these default values for n and d
    n = 10
    d = 5

    # Check for correct problem constraints.
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
        
    if not (2 <= d <= n):
        print(f"Error: Input values n={n}, d={d} do not satisfy 2 <= d <= n.")
        return

    if d % 2 != 1:
        print(f"Error: The value of d ({d}) must be odd.")
        return

    print(f"Calculating the complexity for n = {n} and d = {d}:")
    print("The complexity is given by the formula: 2 + C(n, 1) + C(n, 2) + ... + C(n, d-1)")

    # Calculate and print each term of the sum
    sum_of_widths = 0
    equation_parts = ["2"]
    
    for i in range(1, d):
        width = math.comb(n, i)
        print(f"Term for i={i}: C({n}, {i}) = {width}")
        sum_of_widths += width
        equation_parts.append(str(width))

    complexity = 2 + sum_of_widths
    
    equation_str = " + ".join(equation_parts)

    print("\nFinal Complexity Calculation:")
    print(f"Complexity = {equation_str}")
    print(f"Total Complexity = {complexity}")

if __name__ == "__main__":
    main()