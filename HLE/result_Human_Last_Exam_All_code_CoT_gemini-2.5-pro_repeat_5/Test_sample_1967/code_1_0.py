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

    print(f"For an n-qubit system where n = {n}:")
    print("-" * 40)

    # 1. Counting choices for the Pauli strings
    print("Step 1: Count the number of valid sets of Pauli strings {P_1, ..., P_n}.")

    # For each D_i, its i-th Pauli operator must anti-commute with Z_i.
    # The choices are X_i or Y_i. There are 2 choices for each of the n operators.
    num_diag_choices = 2**n
    print(f"  - The 'diagonal' condition {{S_i, D_i}} = 0 requires the i-th Pauli of D_i to be X or Y.")
    print(f"    Number of choices for these n operators: 2^{n} = {num_diag_choices}")

    # For each D_i, its j-th Pauli operator (where j!=i) must commute with Z_j.
    # The choices are I_j or Z_j.
    # The condition [D_i, D_j] = 0 creates a symmetry constraint on these choices.
    # We only need to choose for the upper triangle of a conceptual n x n matrix.
    num_off_diagonal_pairs = n * (n - 1) // 2
    num_off_diag_choices = 2**num_off_diagonal_pairs
    print(f"  - The 'off-diagonal' conditions [S_i, D_j] = 0 and [D_i, D_j] = 0 constrain the other Pauli operators.")
    print(f"    This leads to a symmetric choice of I or Z for {num_off_diagonal_pairs} pairs of positions.")
    print(f"    Number of choices for these operators: 2^({n}*({n}-1)/2) = 2^{num_off_diagonal_pairs} = {num_off_diag_choices}")
    
    total_pauli_choices = num_diag_choices * num_off_diag_choices
    print(f"  - Total number of Pauli string sets = {num_diag_choices} * {num_off_diag_choices} = {total_pauli_choices}")
    print("-" * 40)

    # 2. Counting choices for the phases
    print("Step 2: Count the number of phase choices.")
    # Each of the n operators D_i can have a global phase from {+1, -1, +i, -i}.
    num_phase_choices = 4**n
    print(f"  - Each of the {n} destabilizers can have one of 4 phases (1, -1, i, -i).")
    print(f"  - Total number of phase combinations: 4^{n} = {num_phase_choices}")
    print("-" * 40)

    # 3. Final calculation
    print("Step 3: Calculate the total number of destabilizer sets.")
    total_choices = total_pauli_choices * num_phase_choices
    print(f"  - Total = (Pauli choices) * (Phase choices)")
    print(f"  - Final Answer = {total_pauli_choices} * {num_phase_choices} = {total_choices}")
    print("-" * 40)
    
    # Verify with the final formula 2**(n*(n+5)/2)
    final_exponent = n * (n + 5) // 2
    final_formula_result = 2**final_exponent
    print(f"The general formula is 2^(n*(n+5)/2).")
    print(f"For n={n}, this is 2^({n}*({n}+5)/2) = 2^{final_exponent} = {final_formula_result}")
    
    # Check that the step-by-step calculation matches the formula
    if total_choices == final_formula_result:
        print("The step-by-step result matches the final formula.")
    else:
        print("Error: Mismatch between step-by-step result and final formula.")


# --- Main execution ---
# Let's use n=3 as an example. You can change this value.
n_qubits = 3
count_destabilizer_sets(n_qubits)
