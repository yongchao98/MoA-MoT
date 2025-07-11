import math
import sys

def calculate_P(n):
    """
    Calculates the value of P(n) for a given positive integer n
    using the derived closed-form formula.

    The formula is: P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Input must be a non-negative integer.", file=sys.stderr)
        return

    # The formula is derived mathematically.
    # The derived closed-form formula for P(n) is:
    # P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)
    
    print(f"The closed-form formula for P(n) is: P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)")
    print(f"\nCalculating for n = {n}:")

    # Handle the base case n=0, where 0! = 1 and the exponent is 0.
    if n == 0:
        result = 1.0
        val_n_factorial = 1
        val_n_minus_1 = -1 # for display
    else:
        # The exponent term
        exponent = n * (n - 1) / 4
        
        # The numerator
        numerator = (2 * math.pi) ** exponent
        
        # The denominator, using math.factorial for n!
        val_n_factorial = math.factorial(n)
        denominator = math.sqrt(val_n_factorial)
        
        result = numerator / denominator
        val_n_minus_1 = n - 1

    # As requested, outputting each number in the final equation.
    print(f"The equation with n={n} is: (2 * pi)^({n} * ({n} - 1) / 4) / sqrt({n}!)")
    
    c1 = 2
    c2 = 4
    c3 = 2  # Denotes the '2' in 1/2 for square root
    val_pi = math.pi
    
    print("\nBreakdown of the calculation components:")
    print(f"The base of the power in the numerator is (c1 * pi) where c1 = {c1} and pi ≈ {val_pi}")
    print(f"The exponent is (n * (n - 1) / c2) where n = {n}, n - 1 = {val_n_minus_1}, and c2 = {c2}")
    print(f"The term in the denominator is (n!)^(1/c3) where n! = {val_n_factorial} and c3 = {c3}")
    
    print("\n---")
    print(f"Final Result: P({n}) = {result}")
    print("---\n")


if __name__ == '__main__':
    # Check if a command-line argument is provided
    if len(sys.argv) > 1:
        try:
            n_input = int(sys.argv[1])
            calculate_P(n_input)
        except ValueError:
            print("Invalid input. Please provide an integer.", file=sys.stderr)
    else:
        # Example usage if no command-line argument is given
        print("No integer n provided via command line. Demonstrating with examples for n = 2, 3, 4.")
        calculate_P(2)
        calculate_P(3)
        calculate_P(4)
