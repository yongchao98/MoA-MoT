import math

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit 
    stabilizer generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits, must be a positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    print(f"For an n-qubit system where n = {n}, the number of possible destabilizer sets is given by the formula 2^(n*(n+1)/2).")
    print("-" * 20)
    print("Here are the steps to calculate the result:")

    # Step 1: Calculate n + 1
    n_plus_1 = n + 1
    print(f"Step 1: Calculate n + 1 = {n} + 1 = {n_plus_1}")

    # Step 2: Calculate the numerator of the exponent
    numerator = n * n_plus_1
    print(f"Step 2: Calculate n * (n + 1) = {n} * {n_plus_1} = {numerator}")

    # Step 3: Calculate the exponent
    exponent = numerator // 2
    print(f"Step 3: Calculate the exponent n * (n + 1) / 2 = {numerator} / 2 = {exponent}")
    
    # Step 4: Calculate the final number. Python's ** operator handles large integers.
    # Note: Using math.pow would return a float, which can lose precision for very large numbers.
    # The ** operator or pow(2, exponent) is preferred for integer arithmetic.
    result = 2**exponent
    print(f"Step 4: Calculate the final result 2^{exponent} = {result}")

    print("-" * 20)
    print(f"The total number of different sets of destabilizers for n = {n} is {result}.")

# You can change the value of n here to see the result for a different number of qubits.
n = 4
count_destabilizer_sets(n)