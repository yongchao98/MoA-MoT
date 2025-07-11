import math

def calculate_destabilizers(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The number of qubits 'n' must be a positive integer.")
        return

    # The number of destabilizer sets is given by the formula:
    # 4^n * 2^(n*(n+1)/2)

    # 1. Calculate the number of choices for global phases.
    # There are n destabilizers, and each can have a phase from {1, -1, i, -i}.
    num_phase_choices = 4**n

    # 2. Calculate the number of choices for the Pauli Z parts.
    # This is the number of n x n symmetric binary matrices.
    # The exponent is n*(n+1)/2.
    exponent_pauli = n * (n + 1) // 2
    num_pauli_choices = 2**exponent_pauli

    # 3. The total number is the product of the two.
    total_destabilizers = num_phase_choices * num_pauli_choices

    print(f"For n = {n} qubits:")
    print("The total number of destabilizer sets is the product of phase choices and Pauli operator choices.")
    print("\nEquation Breakdown:")
    print(f"Number of phase choices = 4^n = 4^{n} = {num_phase_choices}")
    print(f"Number of Pauli operator choices = 2^(n*(n+1)/2) = 2^({n}*({n}+1)/2) = 2^{exponent_pauli} = {num_pauli_choices}")
    
    print(f"\nFinal Calculation:")
    print(f"Total number of sets = {num_phase_choices} * {num_pauli_choices} = {total_destabilizers}")

# Example: Calculate for n=3 qubits. You can change this value.
n_qubits = 3
calculate_destabilizers(n_qubits)
