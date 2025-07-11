import sys

def calculate_correlation_norm(n_str: str):
    """
    This function calculates the 1-norm of the correlation matrix T for a quantum state J_n,
    where n is a positive odd integer.

    The 1-norm of the correlation matrix T for the state J_n has been derived to be 2^(n+1) - 1.
    This script computes this value for a given odd integer n.

    Args:
        n_str: A string representing a positive odd integer for n.
    """
    try:
        n = int(n_str)
        if n < 1 or n % 2 == 0:
            print(f"Error: Input 'n' must be a positive odd integer. Received: {n}")
            return
    except ValueError:
        print(f"Error: Input 'n' must be an integer. Received: '{n_str}'")
        return

    # The formula for the 1-norm is 2^(n+1) - 1.
    # The numbers in the final equation are the base (2), the exponent (n+1), and the subtrahend (1).
    base = 2
    exponent_val = n + 1
    subtrahend = 1
    
    # Calculate the result
    result = base**exponent_val - subtrahend
    
    # Output the details of the calculation as requested
    print(f"For the given odd integer n = {n}:")
    print(f"The 1-norm of the correlation matrix T is given by the formula: 2^(n+1) - 1")
    print(f"Substituting the value of n into the formula gives:")
    print(f"{base}^({n} + 1) - {subtrahend} = {result}")
    print(f"\nThe final result is: {result}")


if __name__ == "__main__":
    # The problem asks for the norm for an arbitrary odd n.
    # This script can take any odd n as a command-line argument.
    # If no argument is provided, it will run with a default example value, n=5.
    if len(sys.argv) > 1:
        input_n = sys.argv[1]
    else:
        print("No command-line argument for n provided. Using default value n=5 for demonstration.")
        input_n = "5"
    
    calculate_correlation_norm(input_n)
