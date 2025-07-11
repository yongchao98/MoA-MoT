import numpy as np

def solve_quantum_fidelity():
    """
    Solves the quantum fidelity problem based on the rules of universe U.
    """
    # Define the observable operator O
    O = np.array([[3, 1],
                  [1, 2]])

    # Step 1: Find the eigenvalues and eigenvectors of O
    # eigvals contains the eigenvalues, eigvecs contains the corresponding
    # eigenvectors as columns.
    eigvals, _ = np.linalg.eig(O)

    # Step 2: Sort eigenvalues in descending order to identify the largest (lambda_1)
    # and second-largest (lambda_2).
    eigvals_sorted = np.sort(eigvals)[::-1]
    lambda_1 = eigvals_sorted[0]
    lambda_2 = eigvals_sorted[1]

    # Step 3: Calculate the terms needed for the fidelity equation.
    # The fidelity F is given by lambda_2^6 / (lambda_1^6 + lambda_2^6).
    # Since lambda_1 and lambda_2 are positive, we don't need absolute values.
    lambda_1_pow_6 = lambda_1**6
    lambda_2_pow_6 = lambda_2**6
    
    sum_of_powers = lambda_1_pow_6 + lambda_2_pow_6

    # Calculate the final fidelity
    fidelity = lambda_2_pow_6 / sum_of_powers

    # Step 4: Print the results as requested
    print("This problem is solved by interpreting the rules of universe U to mean that the final state")
    print("is determined solely by the observable O, making the initial state irrelevant. The fidelity F is")
    print("calculated using the formula derived from the rules:")
    print("\nF = (lambda_2^6) / (lambda_1^6 + lambda_2^6)\n")

    print(f"The second-largest eigenvalue is lambda_2 = {lambda_2:.6f}")
    print(f"The largest eigenvalue is lambda_1 = {lambda_1:.6f}\n")
    
    print("Values for the final equation:")
    print(f"Numerator (lambda_2^6): {lambda_2_pow_6}")
    print(f"Denominator (lambda_1^6 + lambda_2^6): {sum_of_powers}")
    
    print("\nFinal calculated fidelity:")
    print(f"F = {fidelity}")

solve_quantum_fidelity()