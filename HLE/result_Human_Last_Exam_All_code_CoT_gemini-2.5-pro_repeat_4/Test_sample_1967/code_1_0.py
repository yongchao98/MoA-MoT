import sys

def count_destabilizer_sets():
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer set {Z_1, ..., Z_n}.

    The number of such sets is equal to the number of n x n symmetric binary matrices,
    which is 2^(n*(n+1)/2). This script calculates this value.
    """
    try:
        # Prompt for user input for the number of qubits, n.
        n_str = input("Enter the number of qubits (n): ")
        n = int(n_str)
        if n < 1:
            print("Error: Number of qubits (n) must be a positive integer.", file=sys.stderr)
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.", file=sys.stderr)
        return
    except EOFError:
        print("Error: No input provided for n.", file=sys.stderr)
        return

    # Calculate the components of the formula 2**(n * (n + 1) / 2)
    n_plus_1 = n + 1
    numerator = n * n_plus_1
    exponent = numerator // 2

    # Python's ** operator handles arbitrarily large integers, which is necessary
    # as this number grows very quickly with n.
    total_sets = 2**exponent

    print(f"For n = {n}, the number of different destabilizer sets is given by the formula 2**(n * (n + 1) / 2).")
    print("Here is the breakdown of the calculation:")
    print(f"n = {n}")
    print(f"n + 1 = {n_plus_1}")
    print(f"n * (n + 1) = {n} * {n_plus_1} = {numerator}")
    print(f"The exponent is {numerator} / 2 = {exponent}")
    print(f"The final number of different destabilizer sets is 2**{exponent}, which is:")
    print(total_sets)

if __name__ == "__main__":
    count_destabilizer_sets()