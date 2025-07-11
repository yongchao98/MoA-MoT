def count_destabilizer_sets(n: int):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.

    Args:
        n: The number of qubits.

    Returns:
        The total number of different destabilizer sets.
    """
    # The number of choices for the Pauli string part of the destabilizers corresponds
    # to the number of n x n symmetric binary matrices.
    # The formula for this is 2^(n * (n+1) / 2).
    # Using integer division // for safety.
    pauli_choices_exponent = n * (n + 1) // 2
    num_pauli_choices = 2**pauli_choices_exponent

    # For each of the n destabilizers, there are 4 independent phase choices
    # from {+1, -1, +i, -i}. This gives 4^n combinations.
    num_phase_choices = 4**n

    # The total number of sets is the product of the two counts.
    total_count = num_pauli_choices * num_phase_choices

    # Print the breakdown of the calculation as requested.
    print(f"For an n={n} qubit system:")
    print(f"The number of valid sets of Pauli strings is 2^({n}*({n}+1)/2), which is {num_pauli_choices}.")
    print(f"The number of ways to assign phases is 4^{n}, which is {num_phase_choices}.")
    print("\nThe total number of different destabilizer sets is given by the final equation:")
    print(f"{num_pauli_choices} * {num_phase_choices} = {total_count}")
    
    return total_count

if __name__ == '__main__':
    # We choose a specific n for a concrete example, e.g., n=4.
    n_qubits = 4
    count_destabilizer_sets(n_qubits)
