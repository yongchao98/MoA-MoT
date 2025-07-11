import sys

def calculate_approximation():
    """
    Calculates and prints the approximated sum for a given n.
    The value of n is taken from command-line arguments or defaults to 10.
    """
    try:
        # Check for command-line argument for n
        if len(sys.argv) > 1:
            n = int(sys.argv[1])
        else:
            # Default value for n if no argument is provided
            n = 10
        
        # n must be a positive integer
        if n <= 0:
            raise ValueError("n must be a positive integer.")

    except (ValueError, IndexError) as e:
        print(f"Error: Invalid input. {e}", file=sys.stderr)
        print("Usage: python your_script_name.py <positive_integer_n>", file=sys.stderr)
        print("Using default value n=10 instead.", file=sys.stderr)
        n = 10

    # The formula for the approximation is (n^2 / 2) + (1 / 120) + (1 / (252 * n)).
    # This is derived from the Euler-Maclaurin formula and provides an error of order O(n^-2).

    # Calculate the three terms of the formula using floating-point division.
    term1 = n * n / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # Calculate the final approximated sum.
    result = term1 + term2 + term3

    # The problem requires outputting each number in the final equation.
    # We print the equation with the calculated values of the terms.
    print(f"{term1} + {term2} + {term3} = {result}")

if __name__ == "__main__":
    calculate_approximation()