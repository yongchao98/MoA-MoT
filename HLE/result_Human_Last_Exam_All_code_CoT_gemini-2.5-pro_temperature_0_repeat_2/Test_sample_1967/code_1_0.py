import sys

def count_destabilizer_sets():
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.
    The number of qubits 'n' is provided as a command-line argument.
    """
    # --- Input Handling ---
    if len(sys.argv) != 2:
        print("Usage: python <script_name>.py <n>")
        print("Error: Please provide the number of qubits 'n' as a single integer argument.")
        sys.exit(1)

    try:
        n = int(sys.argv[1])
        if n < 1:
            print("Error: The number of qubits 'n' must be a positive integer.")
            sys.exit(1)
    except ValueError:
        print("Error: The provided argument is not a valid integer.")
        sys.exit(1)

    # --- Calculation ---
    # The number of sets is given by the formula 2^(n * (n + 1) / 2).
    # We calculate the exponent first.
    n_plus_1 = n + 1
    numerator = n * n_plus_1
    # Use integer division which is equivalent to / 2 for this case
    exponent = numerator // 2
    
    # Calculate the final result. Python's integers handle arbitrary size.
    result = 2**exponent

    # --- Output ---
    print(f"For an n-qubit system where n = {n}:")
    print("The number of different destabilizer sets for the stabilizers {Z_1, ..., Z_n} is given by the formula:")
    print("  N = 2^(n * (n + 1) / 2)")
    
    print("\nCalculation steps:")
    # This part fulfills the requirement: "output each number in the final equation"
    print(f"  n = {n}")
    print(f"  n + 1 = {n_plus_1}")
    print(f"  n * (n + 1) = {n} * {n_plus_1} = {numerator}")
    print(f"  Exponent = (n * (n + 1)) / 2 = {numerator} / 2 = {exponent}")
    
    print("\nFinal Result:")
    print(f"  N = 2^{exponent}")
    
    # Only print the full number if it's not excessively large to read
    if exponent <= 64:
        print(f"  N = {result}")
    else:
        # Calculate number of digits for context
        num_digits = int(exponent * 0.30103) + 1
        print(f"  (The result is a number with {num_digits} digits, too large to display fully).")


if __name__ == "__main__":
    count_destabilizer_sets()