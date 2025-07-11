import math

def solve_destabilizers(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The number of qubits (n) must be a positive integer.")
        return

    # 1. Calculate the number of possible Pauli string structures.
    # This is equivalent to the number of n x n symmetric binary matrices.
    # The number of elements to choose freely is n (diagonal) + n*(n-1)/2 (upper triangle).
    # This sums to n*(n+1)/2.
    num_pauli_string_configs = 2**(n * (n + 1) // 2)

    # 2. Calculate the number of phase configurations.
    # Each of the n destabilizers can have a phase from {1, -1, i, -i}.
    num_phase_configs = 4**n

    # 3. The total number of sets is the product of the two.
    total_sets = num_pauli_string_configs * num_phase_configs

    # Output the results, including the numbers in the final equation as requested.
    print(f"For n = {n} qubits:")
    print("-" * 30)

    # Explanation for the Pauli part
    pauli_exponent = n * (n + 1) // 2
    print("Step 1: Count Pauli String Structures")
    print(f"The structure is determined by an n x n symmetric binary matrix.")
    print(f"The number of such matrices is 2^(n*(n+1)/2) = 2^{pauli_exponent}.")
    # To handle potentially large numbers, we don't print the expanded number if it's too big
    if pauli_exponent < 100:
        print(f"Number of Pauli string configurations = {num_pauli_string_configs}")
    else:
        print(f"Number of Pauli string configurations = 2^{pauli_exponent}")
    print()

    # Explanation for the Phase part
    print("Step 2: Count Phase Configurations")
    print(f"Each of the {n} destabilizers can have one of 4 phases.")
    print(f"The number of phase configurations is 4^n = 4^{n}.")
    if n < 50:
         print(f"Number of phase configurations = {num_phase_configs}")
    else:
        print(f"Number of phase configurations = 4^{n} (which is 2^{2*n})")
    print()
    
    # Final Equation and Result
    print("Step 3: Calculate Total")
    print("The final equation is: Total Sets = (Number of Pauli String Configs) * (Number of Phase Configs)")
    if pauli_exponent + 2*n < 100:
        print(f"Total number of destabilizer sets = {num_pauli_string_configs} * {num_phase_configs} = {total_sets}")
    else:
        total_exponent = pauli_exponent + 2 * n
        print(f"Total number of destabilizer sets = 2^{pauli_exponent} * 4^{n} = 2^{{{total_exponent}}}")


# --- Main Execution ---
# You can change the value of n here
n_qubits = 3
solve_destabilizers(n_qubits)
