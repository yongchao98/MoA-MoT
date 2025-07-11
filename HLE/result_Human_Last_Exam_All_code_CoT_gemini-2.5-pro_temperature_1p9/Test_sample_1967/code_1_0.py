import math

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer 
    generator set {Z_1, ..., Z_n}.

    The formula is derived as follows:
    1. A destabilizer D_j must be of the form c_j * P_j, where P_j is a Pauli string
       and c_j is a phase.
    2. The commutation relations {S_i, D_i}=0 and [S_i, D_j]=0 for i!=j (with S_i=Z_i)
       fix the X-part of the Pauli strings P_j. Specifically, P_j must contain X_j
       and must not contain any other X_k or Y_k terms (for k!=j).
    3. The commutation relation [D_i, D_j]=0 constrains the Z-part of the strings P_j.
       The matrix of coefficients for the Z parts must be symmetric. The number of such
       n x n binary symmetric matrices is 2^(n*(n+1)/2). This is the number of possible
       sets of Pauli strings {P_1, ..., P_n}.
    4. For each of the n operators in a set, we can choose a global phase from {1, i, -1, -i}.
       This gives 4^n choices for the phases.
    5. The total number of different destabilizer sets is the product of these two numbers:
       (4^n) * (2^(n*(n+1)/2)).
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # Calculate the number of choices for the Pauli strings
    # This corresponds to the number of n x n symmetric binary matrices
    exponent_pauli = n * (n + 1) / 2
    
    # In Python, using integer division // for safety, though n*(n+1) is always even
    if n*(n+1) % 2 != 0:
        print("Error: Exponent for Pauli strings is not an integer. This should not happen.")
        return
    exponent_pauli_int = n * (n + 1) // 2
    string_choices = 2**exponent_pauli_int

    # Calculate the number of choices for the phases
    # Each of the n destabilizers can have one of 4 phases
    exponent_phase = n
    phase_choices = 4**exponent_phase
    
    # Calculate the total number of sets
    total_count = phase_choices * string_choices
    
    # Print the breakdown of the final equation as requested
    print(f"For n = {n}:")
    print("The number of possible sets of destabilizer Pauli strings is:")
    print(f"2^({n}*({n}+1)/2) = 2^{exponent_pauli_int} = {string_choices}")
    print("\nThe number of ways to assign phases to this set is:")
    print(f"4^{exponent_phase} = {phase_choices}")
    print("\nThe total number of different destabilizer sets is the product of these two numbers.")
    print("Final Equation:")
    print(f"({4}^{exponent_phase}) * (2^({n}*({n}+1)/2)) = {total_count}")


# --- Main execution ---
# Set the number of qubits 'n' here.
n_qubits = 3
count_destabilizer_sets(n_qubits)
