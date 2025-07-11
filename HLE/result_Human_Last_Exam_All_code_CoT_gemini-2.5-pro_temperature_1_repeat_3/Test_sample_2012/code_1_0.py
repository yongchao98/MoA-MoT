import numpy as np

def solve_quantum_fidelity():
    """
    Calculates the fidelity of the post-measurement state in a hypothetical universe.
    """
    # Define the observable operator O
    O = np.array([[3, 1],
                  [1, 2]])

    # Step 1: Find the eigenvalues of the operator O.
    # np.linalg.eigvalsh returns sorted eigenvalues for a Hermitian matrix.
    eigenvalues = np.linalg.eigvalsh(O)

    # In a 2-dimensional space, the second-largest eigenvalue is the smallest one.
    lambda_2 = eigenvalues[0] # The second-largest (smallest) eigenvalue
    lambda_1 = eigenvalues[1] # The largest eigenvalue

    # Step 2: Calculate the terms for the fidelity formula.
    # The eigenvalues of O are real and positive, so |lambda| = lambda.
    # The formula for fidelity is F = (lambda_2^6) / (lambda_1^6 + lambda_2^6).
    lambda_2_pow_6 = lambda_2**6
    lambda_1_pow_6 = lambda_1**6

    # Step 3: Calculate the final fidelity.
    fidelity = lambda_2_pow_6 / (lambda_1_pow_6 + lambda_2_pow_6)

    # Print the numbers used in the final equation as requested.
    print(f"The final fidelity is calculated using the formula: F = Numerator / (Denominator + Numerator)")
    print(f"The second-largest eigenvalue (lambda_2) is: {lambda_2:.4f}")
    print(f"The largest eigenvalue (lambda_1) is: {lambda_1:.4f}")
    print(f"Numerator (lambda_2^6): {lambda_2_pow_6:.4f}")
    print(f"Part of Denominator (lambda_1^6): {lambda_1_pow_6:.4f}")
    print(f"Final Fidelity: {fidelity:.8f}")

    # Return the final answer for the wrapper
    return fidelity

# Execute the function
final_answer = solve_quantum_fidelity()
# The final answer is printed in the required format below.
# This is just for capturing the final value.
# The user sees the output from the print statements above.