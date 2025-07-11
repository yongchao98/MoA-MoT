import math

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    print(f"Calculating the number of destabilizer sets for n = {n} qubits.")
    print("-" * 50)

    # 1. Calculate the number of possible Pauli operator structures.
    # This is equal to the number of n x n symmetric binary matrices.
    # The exponent is n (for the diagonal) + n*(n-1)/2 (for the upper triangle).
    # This simplifies to n*(n+1)/2.
    num_pauli_exponent = n * (n + 1) // 2
    
    # Use standard integer exponentiation for potentially large numbers
    num_pauli_structures = 2**num_pauli_exponent

    print(f"Step 1: Find the number of valid structures for the Pauli operators.")
    print(f"This is given by the formula 2^(n*(n+1)/2).")
    print(f"For n={n}, the calculation is: 2^({n}*({n}+1)/2) = 2^{num_pauli_exponent}")
    # Printing large numbers can be slow/long, let's only do it if it's manageable
    if num_pauli_exponent < 100:
        print(f"Result: {num_pauli_structures}")
    else:
        print(f"Result is 2^{num_pauli_exponent}, a very large number.")
    print("-" * 50)

    # 2. Calculate the number of phase combinations.
    # Each of the n destabilizers can have a phase from {+1, -1, +i, -i}.
    # This gives 4^n choices.
    num_phase_choices = 4**n

    print(f"Step 2: Find the number of phase combinations.")
    print(f"This is given by the formula 4^n.")
    print(f"For n={n}, the calculation is: 4^{n} = {num_phase_choices}")
    print("-" * 50)

    # 3. Calculate the total number of destabilizer sets.
    # This is the product of the results from Step 1 and Step 2.
    # Total = 2^(n(n+1)/2) * 4^n = 2^(n(n+1)/2) * 2^(2n) = 2^((n^2+n+4n)/2)
    total_exponent = num_pauli_exponent + 2 * n
    total_count = 2**total_exponent

    print(f"Step 3: Calculate the total number of sets.")
    print(f"The final equation is: (Number of Pauli Structures) * (Number of Phase Choices)")
    print(f"Calculation: {num_pauli_structures} * {num_phase_choices} = 2^{num_pauli_exponent} * 4^{n} = 2^{num_pauli_exponent} * 2^{2*n} = 2^{total_exponent}")
    if total_exponent < 100:
      print(f"Final Answer: {total_count}")
    else:
      print(f"Final Answer is 2^{total_exponent}, a very large number.")

if __name__ == '__main__':
    # You can change this value to calculate for a different number of qubits.
    n_qubits = 3
    count_destabilizer_sets(n_qubits)