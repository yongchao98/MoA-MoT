def calculate_destabilizer_sets(n):
    """
    Calculates the number of destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The number of qubits (n) must be a positive integer.")
        return

    # The formula for the number of sets is 2^(n * (n + 1) / 2)
    # We calculate the parts of the formula step-by-step.
    
    # 1. Calculate n + 1
    n_plus_1 = n + 1
    
    # 2. Calculate the exponent n * (n + 1) / 2
    # Using integer division // to ensure the result is an integer
    exponent = n * n_plus_1 // 2
    
    # 3. Calculate the final result 2^exponent
    # Python's integers can handle arbitrarily large numbers
    result = 2**exponent

    print(f"For an n-qubit system where n = {n}:")
    print("The number of different destabilizer sets is given by the formula: 2^(n * (n + 1) / 2)")
    print("\n--- Calculation Breakdown ---")
    print(f"n = {n}")
    print(f"n + 1 = {n_plus_1}")
    print(f"The exponent is: ({n} * {n_plus_1}) / 2 = {exponent}")
    print(f"The total number of sets is: 2^{exponent} = {result}")

# --- User-runnable example ---
# You can change the value of n_qubits to see the result for different system sizes.
n_qubits = 5
calculate_destabilizer_sets(n_qubits)