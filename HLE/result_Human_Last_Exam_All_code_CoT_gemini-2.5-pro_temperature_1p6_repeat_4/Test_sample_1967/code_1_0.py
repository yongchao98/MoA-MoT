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
    if not isinstance(n, int) or n <= 0:
        print("Error: The number of qubits 'n' must be a positive integer.")
        return

    # Number of ways to choose the Pauli strings.
    # This corresponds to the number of n x n symmetric binary matrices.
    pauli_string_exponent = n * (n + 1) // 2
    num_pauli_choices = 2**pauli_string_exponent

    # Number of ways to choose the phases for the n destabilizers.
    # Each of the n destabilizers can have a phase from {1, -1, i, -i}.
    num_phase_choices = 4**n

    # Total number of different destabilizer sets
    total_count = num_phase_choices * num_pauli_choices

    # Output the result with the calculation steps
    print(f"For n = {n} qubits:")
    print(f"The number of choices for the Pauli strings is 2^({n}*({n}+1)/2) = 2^{pauli_string_exponent} = {num_pauli_choices}.")
    print(f"The number of choices for the phases is 4^{n} = {num_phase_choices}.")
    print("The total number of different destabilizer sets is:")
    print(f"({num_phase_choices}) * ({num_pauli_choices}) = {total_count}")
    
    # An alternative simplified calculation
    # total_exponent = (n**2 + 5*n) / 2
    # print(f"\nAlternatively, using the simplified formula 2^((n^2+5*n)/2):")
    # print(f"2^(({n}^2+5*{n})/2) = 2^({int(total_exponent)}) = {total_count}")


# --- User Input ---
# Set the number of qubits 'n' here.
# For example, for a 3-qubit system:
n = 3
# ------------------

count_destabilizer_sets(n)
