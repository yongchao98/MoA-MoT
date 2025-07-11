def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.

    Returns:
        int: The total number of different destabilizer sets.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The formula for the number of sets is 2^(3*n + n*(n-1)/2)
    # Calculate the terms in the exponent
    phase_and_diag_part = 3 * n
    off_diag_part_numerator = n * (n - 1)
    
    # This division will always result in an integer because n or n-1 is even
    off_diag_part = off_diag_part_numerator // 2
    
    exponent = phase_and_diag_part + off_diag_part
    
    # The result can be very large, Python handles large integers automatically
    result = 2**exponent
    
    print(f"For n = {n}, the number of different sets of destabilizers is given by the formula: 2^(3*n + n*(n-1)/2)")
    print(f"Substituting n = {n} into the formula:")
    print(f"Number of sets = 2^(3*{n} + ({n}*({n}-1))/2)")
    print(f"               = 2^({phase_and_diag_part} + {off_diag_part})")
    print(f"               = 2^{exponent}")
    print(f"               = {result}")

# Example usage for n=4
count_destabilizer_sets(4)