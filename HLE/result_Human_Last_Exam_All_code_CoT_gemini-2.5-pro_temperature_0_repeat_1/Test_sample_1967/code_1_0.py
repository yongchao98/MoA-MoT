import math

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits. Must be a non-negative integer.

    Returns:
        int: The total number of different destabilizer sets.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return None

    # The number of sets is 2^(n * (n + 1) / 2).
    # We calculate the exponent first.
    # The formula for the number of elements in the upper triangle of a
    # symmetric matrix, including the diagonal, is n * (n + 1) / 2.
    exponent = n * (n + 1) // 2
    
    # The total number is 2 to the power of this exponent.
    # Using 1 << exponent is an efficient way to calculate 2**exponent for integers.
    total_sets = 1 << exponent
    
    print(f"For n = {n} qubits:")
    # Print the formula with the value of n substituted.
    print(f"The number of different destabilizer sets is 2^({n} * ({n} + 1) / 2)")
    print(f"= 2^{exponent}")
    print(f"= {total_sets}")
    print("-" * 20)
    return total_sets

if __name__ == '__main__':
    # Demonstrate the function for a few values of n.
    count_destabilizer_sets(1)
    count_destabilizer_sets(2)
    count_destabilizer_sets(3)
    count_destabilizer_sets(4)
    # Example with a larger n
    count_destabilizer_sets(10)