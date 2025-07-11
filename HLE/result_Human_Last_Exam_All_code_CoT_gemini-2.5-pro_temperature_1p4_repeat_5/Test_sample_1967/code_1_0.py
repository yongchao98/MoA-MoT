import math

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.

    Returns:
        int: The total number of different destabilizer sets.
             Returns -1 if n is not a positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The number of qubits 'n' must be a positive integer.")
        return -1

    # The number of ways to choose the Pauli strings corresponds to the number
    # of n x n symmetric binary matrices.
    # The exponent is n for the diagonal elements plus n*(n-1)/2 for the
    # upper-triangular elements, which simplifies to n*(n+1)/2.
    exponent_pauli = n * (n + 1) // 2
    num_pauli_strings = 2**exponent_pauli

    # For each of the n generators, there are 4 choices for the global phase.
    num_phases = 4**n

    # The total number of sets is the product of the two.
    total_count = num_pauli_strings * num_phases
    
    print(f"For n = {n} qubits:")
    print("-" * 30)
    print(f"The number of possible sets of Pauli strings is:")
    print(f"  2^(n*(n+1)/2) = 2^({n}*({n}+1)/2) = 2^{exponent_pauli} = {num_pauli_strings}")
    print(f"The number of choices for global phases is:")
    print(f"  4^n = 4^{n} = {num_phases}")
    print("-" * 30)
    print(f"Total number of destabilizer sets = {num_pauli_strings} * {num_phases}")
    print(f"Final Answer: {total_count}")
    return total_count

if __name__ == '__main__':
    # Let's calculate the result for a 4-qubit system as an example.
    num_qubits = 4
    count_destabilizer_sets(num_qubits)