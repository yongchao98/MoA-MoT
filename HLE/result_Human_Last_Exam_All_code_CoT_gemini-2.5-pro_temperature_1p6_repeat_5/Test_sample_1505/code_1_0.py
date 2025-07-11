import sys

def calculate_approximation():
    """
    Calculates the sum approximation for a given integer n.
    The script expects n as a command-line argument.
    """
    if len(sys.argv) != 2:
        print("Usage: python your_script.name <n>", file=sys.stderr)
        sys.exit(1)

    try:
        n_str = sys.argv[1]
        n = int(n_str)
        if n <= 0:
            raise ValueError("Input must be a positive integer.")
    except ValueError as e:
        print(f"Error: Invalid input. {e}", file=sys.stderr)
        sys.exit(1)

    # The derived approximation formula is S_n ≈ n²/2 + 1/120 + 1/(252*n).
    # The numbers in the equation are 2, 120, and 252.
    num_2 = 2
    num_120 = 120
    num_252 = 252

    # Calculate the result using floating-point arithmetic.
    result = n**2 / float(num_2) + 1 / float(num_120) + 1 / (float(num_252) * n)

    # Output the final equation with the value of n, showing each number.
    print(f"{n}^2/{num_2} + 1/{num_120} + 1/({num_252}*{n}) = {result}")

if __name__ == '__main__':
    calculate_approximation()