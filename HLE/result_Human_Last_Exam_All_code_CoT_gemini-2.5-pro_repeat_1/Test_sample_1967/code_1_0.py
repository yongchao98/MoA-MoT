import math

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.

    Returns:
        int: The total number of different destabilizer sets.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # The formula for the number of sets is 2^((n^2 + 5n) / 2)
    print(f"For n = {n}, we calculate the number of destabilizer sets.")
    print("The formula is: 2**((n**2 + 5*n) / 2)")
    print("-" * 30)

    # Calculate each part of the exponent
    n_squared = n * n
    five_n = 5 * n
    numerator = n_squared + five_n
    
    # The exponent is guaranteed to be an integer for any integer n
    exponent = numerator // 2

    # Calculate the final result
    # Using Python's arbitrary-precision integers, this works for large n
    total_sets = pow(2, exponent)

    # Print the breakdown of the calculation
    print("Calculation of the exponent (n^2 + 5*n) / 2:")
    print(f"n^2 = {n} * {n} = {n_squared}")
    print(f"5*n = 5 * {n} = {five_n}")
    print(f"n^2 + 5*n = {n_squared} + {five_n} = {numerator}")
    print(f"Exponent = {numerator} // 2 = {exponent}")
    print("-" * 30)
    print(f"Total number of destabilizer sets = 2**{exponent}")
    print(f"Result: {total_sets}")


if __name__ == '__main__':
    # You can change the value of n here
    n_qubits = 3
    count_destabilizer_sets(n_qubits)