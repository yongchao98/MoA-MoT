import sys

def calculate_max_zeros(N):
    """
    Calculates the maximal possible number of complex zeros for the determinant
    of the matrix B(k).

    Args:
      N: The dimension of the matrix, a positive integer.
    """
    if not isinstance(N, int) or N <= 0:
        print("Error: N must be a positive integer.", file=sys.stderr)
        return

    # Based on the derivation, the maximal number of complex zeros is given by the formula N * 2^N.
    # This arises from the degree of the polynomial system describing the problem.
    # The number of solutions for k_1^2 is N * 2^(N-1).
    # Each complex solution for k_1^2 gives two complex solutions for k_1.
    
    # Calculate 2^N
    power_of_2 = 2**N
    
    # Calculate the final result
    result = N * power_of_2
    
    print(f"For a matrix of size N x N (where N = {N}):")
    print("The formula for the maximal number of complex zeros is: N * 2^N")
    print("\nCalculation steps:")
    
    # As requested, outputting each number in the final equation.
    print(f"1. The value of N is {N}.")
    print(f"2. Calculate 2 to the power of N: 2^{N} = {power_of_2}.")
    print(f"3. Multiply N by the result from step 2: {N} * {power_of_2} = {result}.")
    
    print(f"\nThe maximal possible number of complex zeros is: {result}")

if __name__ == '__main__':
    # This script can be run from the command line with N as an argument.
    # e.g., python your_script_name.py 4
    if len(sys.argv) > 1:
        try:
            N_val = int(sys.argv[1])
            calculate_max_zeros(N_val)
        except ValueError:
            print("Error: Please provide a valid integer for N.", file=sys.stderr)
            print("Usage: python script.py <N>", file=sys.stderr)
    else:
        # If no argument is given, run with an example value.
        print("No command-line argument for N found.")
        print("Running with a default example value of N = 3.")
        print("To use a different value, run as: python script.py <N>\n")
        calculate_max_zeros(3)