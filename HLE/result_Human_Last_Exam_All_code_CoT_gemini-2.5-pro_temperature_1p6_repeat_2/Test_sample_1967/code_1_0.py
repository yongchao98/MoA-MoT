import math

def count_destabilizer_sets_for_zi_basis(n):
    """
    Calculates and explains the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.
    """
    if not isinstance(n, int) or n < 1:
        print("Please provide a positive integer for n.")
        return

    print(f"Derivation for an n = {n} qubit system.")
    print("==========================================")
    print("Given the stabilizer generators S_i = Z_i for i = 1, ..., n.")
    print("A set of destabilizer generators {D_1, ..., D_n} must satisfy:")
    print("1. [D_i, D_j] = 0 (Commute with each other)")
    print("2. [S_i, D_j] = 0 for i != j (Commute with other stabilizers)")
    print("3. {S_i, D_i} = 0 (Anti-commute with its corresponding stabilizer)")

    print("\nStep 1: Determine the structure of each destabilizer generator D_j.")
    print("From {S_j, D_j} = 0, the part of D_j on qubit j must be X_j or Y_j (2 choices).")
    print("From [S_i, D_j] = 0 for i!=j, the part of D_j on qubit i must be I_i or Z_i.")
    
    print("\nStep 2: Count the choices for the anti-commuting parts.")
    print("For each of the n destabilizers D_j, we must choose between X_j or Y_j.")
    num_A_choices = 2**n
    print(f"This gives 2^{n} choices in total.")
    print(f"For n = {n}, this is 2**{n} = {num_A_choices} choices.")

    print("\nStep 3: Count the choices for the commuting parts.")
    print("The choices of I_k or Z_k can be represented by a matrix of binary exponents 'a_ij'.")
    print("The condition [D_i, D_j] = 0 imposes a symmetry constraint: a_ij = a_ji.")
    print("This means we only need to choose the exponents for the upper triangle of the matrix.")
    
    num_independent_exponents = n * (n - 1) // 2
    print(f"The number of independent exponents is n*(n-1)/2.")
    print(f"For n = {n}, this is {n}*({n}-1)/2 = {num_independent_exponents}.")
    
    num_exponent_choices = 2**num_independent_exponents
    print("Since each exponent can be 0 or 1, the number of choices for these exponents is 2**(n*(n-1)/2).")
    print(f"For n = {n}, this is 2**{num_independent_exponents} = {num_exponent_choices}.")
    
    print("\nStep 4: Combine the counts to find the total number of sets.")
    total_sets = num_A_choices * num_exponent_choices
    print("The total number of sets is the product of the choices from Step 2 and Step 3.")
    print(f"Total = (choices for X/Y parts) * (choices for I/Z parts)")
    print(f"      = {num_A_choices} * {num_exponent_choices} = {total_sets}")
    
    print("\nFinal formula check:")
    final_exponent = n * (n + 1) // 2
    final_result_from_formula = 2**final_exponent
    print("The general formula is 2**(n*(n+1)/2).")
    print(f"For n = {n}, the final equation is 2**({n}*({n}+1)/2) = 2**{final_exponent} = {final_result_from_formula}.")
    print(f"Our calculated total {total_sets} matches the formula result {final_result_from_formula}.")


# We will demonstrate the calculation for a 4-qubit system.
n_qubits = 4
count_destabilizer_sets_for_zi_basis(n_qubits)