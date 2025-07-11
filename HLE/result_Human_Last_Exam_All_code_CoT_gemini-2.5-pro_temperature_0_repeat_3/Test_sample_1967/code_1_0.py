import math

def count_destabilizer_sets(n: int):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.

    Args:
        n: The number of qubits.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The number of qubits 'n' must be a positive integer.")
        return

    # The number of different destabilizer sets is given by the formula 2^(n*(n+1)/2).
    # This arises from counting the number of n x n symmetric binary matrices.

    # Step 1: Calculate the term (n + 1)
    n_plus_1 = n + 1

    # Step 2: Calculate the exponent n * (n + 1) / 2.
    # We use integer division // to ensure the result is an integer.
    exponent = n * n_plus_1 // 2

    # Step 3: Calculate the final result, 2 to the power of the exponent.
    # Python's integers handle arbitrary size, so this won't overflow for reasonable n.
    try:
        result = 2**exponent
        is_too_large = False
    except OverflowError:
        result = float('inf')
        is_too_large = True


    # Print the explanation and the breakdown of the calculation
    print(f"For an n-qubit system where n = {n}:")
    print("The stabilizer generator set is {Z_1, ..., Z_n}.")
    print("The number of corresponding destabilizer generator sets is given by the formula: 2^(n * (n + 1) / 2)")
    print("\nCalculation breakdown:")
    print(f"n = {n}")
    print(f"n + 1 = {n_plus_1}")
    print(f"exponent = n * (n + 1) / 2 = {n} * {n_plus_1} / 2 = {exponent}")

    print("\nFinal equation:")
    if is_too_large:
        print(f"2^{exponent} = infinity (The number is too large to be represented as a standard float).")
    else:
        # The final result is printed with the full equation
        print(f"2^{exponent} = {result}")

# --- Example Usage ---
# You can change the value of n here to see the result for a different number of qubits.
n_qubits = 4
count_destabilizer_sets(n_qubits)