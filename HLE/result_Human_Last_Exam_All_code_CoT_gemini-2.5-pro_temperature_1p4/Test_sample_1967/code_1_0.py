import math

def calculate_destabilizers(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits. Must be a positive integer.

    Returns:
        None. Prints the result of the calculation.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: The number of qubits (n) must be a positive integer.")
        return

    # The number of different destabilizer sets is 2^(n*(n+1)/2).
    # First, calculate the exponent.
    exponent = n * (n + 1) // 2

    # Then, calculate the final number.
    # Use standard pow for potentially large numbers.
    result = pow(2, exponent)

    # As requested, output each number in the final equation.
    # The final equation is Result = 2 ^ Exponent.
    print(f"For n = {n} qubits:")
    print(f"The number of sets is given by the formula: 2^(n * (n + 1) / 2)")
    print(f"The calculation is: 2^({n}*({n}+1)/2) = 2^{exponent} = {result}")

# You can change the value of n to calculate for a different number of qubits.
# Example for n = 3:
calculate_destabilizers(3)

# Example for n = 5:
calculate_destabilizers(5)
