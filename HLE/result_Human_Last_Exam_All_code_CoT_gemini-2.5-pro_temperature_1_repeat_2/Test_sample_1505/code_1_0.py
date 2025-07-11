import sys

def calculate_approximation(n_str):
    """
    Calculates the approximation for the sum for a given integer n.

    The approximation formula is S(n) ≈ n^2/2 + 1/120 + 1/(252*n).
    This formula gives an absolute error of O(n^-2).
    """
    try:
        n = int(n_str)
        if n <= 0:
            print("Error: n must be a positive integer.")
            return
    except ValueError:
        print("Error: Please provide a valid integer for n.")
        return

    # Calculate the three terms of the formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # Calculate the final result
    result = term1 + term2 + term3

    # Print the equation with the calculated values
    print(f"For n = {n}, the sum is approximated by the equation:")
    print(f"{term1} + {term2} + {term3} = {result}")


if __name__ == '__main__':
    # The user can provide n as a command-line argument.
    # If not provided, it will use a default value.
    if len(sys.argv) > 1:
        n_value = sys.argv[1]
    else:
        # Default value for n if not provided.
        n_value = "100"
        print(f"No value for n provided, using default n = {n_value}")
        print("You can also run the script with a value for n, e.g., 'python your_script_name.py 50'")
    
    calculate_approximation(n_value)
