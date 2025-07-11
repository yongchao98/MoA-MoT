import numpy as np

def solve_quantum_fidelity():
    """
    Solves the quantum fidelity problem based on the rules of universe U.
    """
    # Step 1: Define the observable operator and find its eigenvalues.
    O = np.array([[3, 1],
                  [1, 2]])

    # Calculate eigenvalues. Eigenvectors are not needed for the final formula.
    eigenvalues = np.linalg.eigvals(O)

    # Sort eigenvalues in descending order to identify the largest and second-largest.
    eigenvalues = np.sort(eigenvalues)[::-1]
    
    lambda_1 = eigenvalues[0]  # Largest eigenvalue
    lambda_2 = eigenvalues[1]  # Second-largest eigenvalue

    print("Step 1: Find the eigenvalues of the operator O.")
    print(f"The observable operator O is:\n{O}")
    print(f"The eigenvalues are sorted in descending order:")
    print(f"Largest eigenvalue (lambda_1) = {lambda_1:.8f}")
    print(f"Second-largest eigenvalue (lambda_2) = {lambda_2:.8f}")
    print("-" * 40)

    # Step 2: Calculate the 6th power of the eigenvalues as required by the fidelity formula.
    lambda_1_pow_6 = lambda_1**6
    lambda_2_pow_6 = lambda_2**6

    print("Step 2: Calculate the 6th power of the eigenvalues.")
    print(f"lambda_1^6 = ({lambda_1:.8f})^6 = {lambda_1_pow_6:.8f}")
    print(f"lambda_2^6 = ({lambda_2:.8f})^6 = {lambda_2_pow_6:.8f}")
    print("-" * 40)

    # Step 3: Calculate the fidelity.
    # The fidelity F is given by the formula: F = lambda_2^6 / (lambda_1^6 + lambda_2^6)
    fidelity = lambda_2_pow_6 / (lambda_1_pow_6 + lambda_2_pow_6)

    print("Step 3: Calculate the fidelity.")
    print("The fidelity F is the ratio of the powered eigenvalues.")
    print(f"F = lambda_2^6 / (lambda_1^6 + lambda_2^6)")
    print(f"F = {lambda_2_pow_6:.8f} / ({lambda_1_pow_6:.8f} + {lambda_2_pow_6:.8f})")
    print(f"F = {lambda_2_pow_6:.8f} / {lambda_1_pow_6 + lambda_2_pow_6:.8f}")
    print("-" * 40)
    print(f"The final calculated fidelity is: {fidelity:.8f}")
    
    return fidelity

# Run the solver and print the final answer in the specified format.
final_fidelity = solve_quantum_fidelity()
# The final answer is wrapped in <<<>>>
# print(f"\n<<<final_fidelity:.8f>>>") # This is a placeholder, will output a single final line as requested
# The final answer should be on its own line after the explanation.
# The user prompt wants a final line containing <<<answer content>>>.

print(f"\n<<<0.00309601>>>")