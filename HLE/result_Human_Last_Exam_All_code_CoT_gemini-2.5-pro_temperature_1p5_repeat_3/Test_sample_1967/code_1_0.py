def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.

    The number of choices for the Pauli operators of the destabilizers
    corresponds to the number of n x n symmetric binary matrices.
    The number of such matrices is 2^(n*(n+1)/2).

    For each of the n destabilizer generators, there are 4 choices for the
    global phase (+1, -1, +i, -i), giving 4^n choices for the phases.

    Args:
        n (int): The number of qubits.

    Returns:
        None. Prints the calculation steps and the result.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: The number of qubits (n) must be a positive integer.")
        return

    # Calculate the exponent for the number of symmetric matrices
    # This is the number of elements in the upper triangle plus the diagonal:
    # n*(n-1)/2 + n = (n^2 - n + 2n)/2 = n*(n+1)/2
    num_pauli_choices_exp = n * (n + 1) // 2

    # Calculate the number of choices for the Pauli part (the symmetric matrix)
    # Using integer arithmetic to handle large numbers
    num_pauli_choices = 1 << num_pauli_choices_exp

    # Calculate the number of phase choices (4^n)
    num_phase_choices_base = 4
    num_phase_choices = num_phase_choices_base**n

    # Total number is the product of the two
    total_choices = num_phase_choices * num_pauli_choices

    print(f"For an n={n} qubit system with stabilizer generators {{Z_1, ..., Z_n}}:")
    print("-" * 30)
    print("1. Counting choices for the Pauli operators (up to phase):")
    print(f"   This is equivalent to counting n x n symmetric binary matrices.")
    print(f"   The number of choices is 2^(n*(n+1)/2) = 2^({n}*({n}+1)/2) = 2^{num_pauli_choices_exp} = {num_pauli_choices}.")
    print("\n2. Counting choices for the global phases:")
    print(f"   Each of the n generators can have a phase from {{+1, -1, +i, -i}}.")
    print(f"   The number of choices is {num_phase_choices_base}^n = {num_phase_choices}.")
    print("\n3. Total number of different destabilizer sets:")
    print("   Final Equation:")
    print(f"   {num_phase_choices_base}^{n} * 2^({n}*({n}+1)/2) = {num_phase_choices} * {num_pauli_choices} = {total_choices}")

# --- Main execution ---
# Set the number of qubits, n.
# You can change this value to calculate the result for a different number of qubits.
n = 3
count_destabilizer_sets(n)
