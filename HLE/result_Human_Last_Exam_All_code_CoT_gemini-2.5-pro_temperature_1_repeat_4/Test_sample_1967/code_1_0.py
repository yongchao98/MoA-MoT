def count_destabilizer_sets():
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.
    """
    try:
        n_str = input("Enter the number of qubits (n): ")
        n = int(n_str)
        if n < 1:
            print("Number of qubits must be a positive integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The formula for the number of sets is 2^((n^2 + 5n) / 2)
    # We calculate the components of the exponent first.
    n_squared = n * n
    five_n = 5 * n
    
    # The exponent must be an integer, so we use integer division.
    numerator = n_squared + five_n
    exponent = numerator // 2
    
    # Calculate the final result.
    # Python's integers handle arbitrary size, so this won't overflow.
    result = 1 << exponent # Same as 2**exponent but can be faster

    print("\nFor n = {}, the number of different destabilizer sets is derived from the formula: 2^((n^2 + 5*n) / 2)".format(n))
    print("Let's calculate the value step by step:")
    print("1. Calculate the terms in the exponent's numerator:")
    print("   n^2 = {}^2 = {}".format(n, n_squared))
    print("   5*n = 5*{} = {}".format(n, five_n))
    print("2. Calculate the numerator:")
    print("   n^2 + 5*n = {} + {} = {}".format(n_squared, five_n, numerator))
    print("3. Calculate the exponent:")
    print("   ({} / 2) = {}".format(numerator, exponent))
    print("4. Calculate the final result:")
    print("   2^{} = {}".format(exponent, result))

if __name__ == '__main__':
    count_destabilizer_sets()