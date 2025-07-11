def solve_destabilizers(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.
    """
    if not isinstance(n, int) or n < 1:
        print("Please provide a positive integer for n.")
        return

    # Part 1: Calculate the number of possible Pauli operator structures.
    # The exponent for the number of choices for the diagonal part (X or Y) is n.
    # The exponent for the number of choices for the off-diagonal part (I or Z pairs)
    # corresponds to the upper triangle of an nxn matrix, which is n*(n-1)/2.
    # Total exponent for Pauli choices = n + n*(n-1)/2 = n*(n+1)/2
    pauli_exponent = n * (n + 1) // 2
    
    # Using python's arbitrary-precision integers to handle large numbers
    num_pauli_choices = 1 << pauli_exponent # Same as 2**pauli_exponent

    # Part 2: Calculate the number of phase choices.
    # For each of the n destabilizers, there are 4 phase choices (+1, -1, +i, -i).
    num_phase_choices = 4**n
    
    # Part 3: The total number of destabilizer sets is the product of the two parts.
    total_choices = num_pauli_choices * num_phase_choices

    print(f"For n = {n} qubits:")
    print(f"The number of ways to choose the set of Pauli operators is 2^({n}({n}+1)/2) = {num_pauli_choices}")
    print(f"The number of ways to choose the set of phases is 4^{n} = {num_phase_choices}")
    print("\nThe final equation for the total number of sets is:")
    print(f"{num_pauli_choices} * {num_phase_choices} = {total_choices}")

    # For verification, the simplified total formula is 2^(n*(n+5)/2)
    final_exponent = n * (n + 5) // 2
    verification_total = 1 << final_exponent
    print(f"\nThis result is consistent with the simplified formula 2^({n}({n}+5)/2) = 2^{final_exponent} = {verification_total}")


# You can change the number of qubits 'n' here to see the result.
n_qubits = 4
solve_destabilizers(n_qubits)