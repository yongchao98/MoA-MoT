import math

def calculate_destabilizer_sets(n: int):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer set {Z_1, ..., Z_n}.

    Args:
        n: The number of qubits.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Number of qubits (n) must be a positive integer.")
        return

    print(f"Calculating the number of destabilizer sets for n = {n} qubits.")
    print("-" * 50)
    print(f"The stabilizer generator set is {{S_i}} = {{Z_1, ..., Z_{n}}}.")
    print("A destabilizer set {D_1, ..., D_n} must satisfy:")
    print("1. [D_i, D_j] = 0 for all i, j")
    print("2. [S_i, D_j] = 0 for i != j")
    print("3. {S_i, D_i} = 0 for all i")
    print("")

    # Step 1: Count the number of possible Pauli string configurations.
    # This is equivalent to counting the number of n x n symmetric binary matrices.
    num_free_elements = n * (n + 1) // 2
    num_pauli_sets = 2**num_free_elements
    
    print("The commutation relations dictate the structure of the destabilizers.")
    print("Each destabilizer D_i must contain an X or Y on qubit i and only I or Z on other qubits.")
    print("The condition [D_i, D_j]=0 means the matrix of Z parts must be symmetric.")
    print(f"The number of such sets of Pauli strings is 2^(n*(n+1)/2).")
    print(f"For n = {n}, this is 2^({n}*({n}+1)/2) = 2^{num_free_elements} = {num_pauli_sets}")
    print("")

    # Step 2: Count the number of phase choices.
    # Each of the n destabilizers can have a phase from {1, -1, i, -i}.
    num_phase_choices = 4**n

    print("Each of the n destabilizers can have one of 4 global phases.")
    print(f"The number of ways to choose these phases is 4^n.")
    print(f"For n = {n}, this is 4^{n} = {num_phase_choices}")
    print("")

    # Step 3: Calculate the total number of sets.
    total_sets = num_pauli_sets * num_phase_choices

    print("The total number of different destabilizer sets is the product of these two numbers.")
    print("Total = (Number of Pauli sets) * (Number of phase choices)")
    print(f"Total = {num_pauli_sets} * {num_phase_choices} = {total_sets}")
    print("-" * 50)


if __name__ == '__main__':
    # Example: Calculate for n=4
    n_qubits = 4
    calculate_destabilizer_sets(n_qubits)
