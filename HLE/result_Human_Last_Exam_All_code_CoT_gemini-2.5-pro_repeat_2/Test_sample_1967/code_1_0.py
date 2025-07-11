def count_destabilizer_sets():
    """
    Calculates the number of different destabilizer sets for an n-qubit system
    stabilized by {Z_1, ..., Z_n}.
    """
    # The number of qubits in the system. You can change this value.
    n = 3

    # The formula for the number of destabilizer sets is 2^(n*(n+1)/2).
    # We calculate the exponent first.
    exponent = n * (n + 1) // 2

    # Calculate the final number, which can be very large.
    # Python's integers handle arbitrary size automatically.
    result = 2**exponent

    # Print the explanation and the step-by-step calculation.
    print(f"For an n-qubit system with n = {n}:")
    print("The number of different destabilizer sets for the stabilizer set {Z_1, ..., Z_n} is given by the formula 2^(n*(n+1)/2).")
    print("\nCalculation:")
    # Show the formula with the value of n plugged in.
    print(f"Number of sets = 2^({n}*({n}+1)/2)")
    # Show the intermediate steps of the calculation.
    print(f"                 = 2^({n * (n + 1)}/2)")
    print(f"                 = 2^{exponent}")
    # Print the final result.
    print(f"                 = {result}")

# Run the function to display the result.
count_destabilizer_sets()