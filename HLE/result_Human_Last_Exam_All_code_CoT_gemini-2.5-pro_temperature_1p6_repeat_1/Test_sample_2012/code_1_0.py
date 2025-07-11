import numpy as np

def solve_quantum_fidelity():
    """
    Calculates the fidelity of the post-measurement state with respect to the
    eigenstate of the second-largest eigenvalue based on the problem's rules.
    """
    # Define the observable operator matrix
    O = np.array([[3, 1], [1, 2]])

    # Step 1: Find the eigenvalues of the operator O
    # numpy.linalg.eig returns eigenvalues and normalized eigenvectors
    eigenvalues, _ = np.linalg.eig(O)

    # Step 2: Identify the largest and second-largest eigenvalues
    # Sort eigenvalues in descending order
    eigenvalues.sort()
    # The largest is the last element, second-largest is the first
    lambda_1 = eigenvalues[1]  # Largest eigenvalue
    lambda_2 = eigenvalues[0]  # Second-largest eigenvalue

    print(f"The eigenvalues of the operator are:")
    print(f"Largest eigenvalue (lambda_1): {lambda_1}")
    print(f"Second-largest eigenvalue (lambda_2): {lambda_2}\n")
    
    # Step 3: Calculate the terms for the fidelity formula
    # The fidelity F is given by lambda_2^6 / (lambda_1^6 + lambda_2^6)
    l1_pow_6 = lambda_1**6
    l2_pow_6 = lambda_2**6
    
    # Step 4: Calculate the final fidelity and print the equation
    fidelity = l2_pow_6 / (l1_pow_6 + l2_pow_6)

    print("The fidelity F is calculated using the formula: F = lambda_2^6 / (lambda_1^6 + lambda_2^6)")
    print("Let's calculate the components:")
    print(f"lambda_1^6 = {lambda_1**6:.4f}")
    print(f"lambda_2^6 = {lambda_2**6:.4f}\n")
    print(f"So, the final equation is:")
    print(f"F = {l2_pow_6:.4f} / ({l1_pow_6:.4f} + {l2_pow_6:.4f})")
    print(f"F = {fidelity:.8f}")

    print("\n<<<{}>>>".format(fidelity))


solve_quantum_fidelity()
