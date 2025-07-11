import numpy as np

def solve_quantum_fidelity():
    """
    Calculates the fidelity based on the problem's unique physical laws.
    """
    # Define the observable operator O
    O = np.array([[3, 1],
                  [1, 2]])

    # Step 1: Find the eigenvalues of the operator O.
    # The eigenvectors are not explicitly needed for the final formula.
    eigenvalues = np.linalg.eigvals(O)

    # Step 2: Sort the eigenvalues to identify the largest and second-largest.
    # np.sort provides ascending order, so we reverse it for descending.
    sorted_eigenvalues = np.sort(eigenvalues)[::-1]
    lambda_1 = sorted_eigenvalues[0]  # Largest eigenvalue
    lambda_2 = sorted_eigenvalues[1]  # Second-largest eigenvalue

    # Step 3: Calculate the terms needed for the fidelity formula F = lambda_2^6 / (lambda_1^6 + lambda_2^6).
    lambda_1_pow_6 = lambda_1**6
    lambda_2_pow_6 = lambda_2**6
    denominator = lambda_1_pow_6 + lambda_2_pow_6

    # Step 4: Calculate the final fidelity.
    fidelity = lambda_2_pow_6 / denominator

    # Step 5: Print the results and the final equation with numerical values.
    print("Step 1: Find the eigenvalues of the operator.")
    print(f"The eigenvalues are: {lambda_1:.6f} and {lambda_2:.6f}")
    print("\nStep 2: The fidelity F is calculated using the formula F = (lambda_2^6) / (lambda_1^6 + lambda_2^6).")
    
    print("\nStep 3: Calculate the numerical values for the equation.")
    print(f"Value of lambda_1^6 = {lambda_1_pow_6:.6f}")
    print(f"Value of lambda_2^6 = {lambda_2_pow_6:.6f}")
    print(f"Value of the denominator (lambda_1^6 + lambda_2^6) = {denominator:.6f}")

    print("\nStep 4: The final equation with these numbers is:")
    print(f"Fidelity = {lambda_2_pow_6:.6f} / ({lambda_1_pow_6:.6f} + {lambda_2_pow_6:.6f})")
    
    print(f"\nFinal calculated fidelity is: {fidelity:.6f}")

    # Output the final answer in the required format
    print(f"\n<<<{fidelity}>>>")

solve_quantum_fidelity()