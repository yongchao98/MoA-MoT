import numpy as np

def solve_quantum_fidelity():
    """
    Calculates the fidelity based on the alternate universe's physical laws.
    """
    # Define the observable operator matrix
    O = np.array([[3, 1],
                  [1, 2]])

    # For a real symmetric matrix, eigvalsh is efficient and returns sorted eigenvalues
    eigenvalues = np.linalg.eigvalsh(O)

    # The eigenvalues are sorted in increasing order.
    # lambda_2 is the second-largest (i.e., the smallest) eigenvalue
    # lambda_1 is the largest eigenvalue
    lambda_2 = eigenvalues[0]
    lambda_1 = eigenvalues[1]

    # Calculate the 6th power of the eigenvalues, as derived from the problem's rules
    lambda_2_pow_6 = lambda_2**6
    lambda_1_pow_6 = lambda_1**6

    # The fidelity F is given by the formula: F = λ2^6 / (λ1^6 + λ2^6)
    denominator = lambda_1_pow_6 + lambda_2_pow_6
    fidelity = lambda_2_pow_6 / denominator

    # Print the components of the final equation as requested
    print("The fidelity F is calculated using the formula: F = (λ2^6) / (λ1^6 + λ2^6)\n")
    print("Calculated values for the equation:")
    print(f"λ2 (second-largest eigenvalue) = {lambda_2}")
    print(f"λ1 (largest eigenvalue) = {lambda_1}\n")
    
    print("Final Equation:")
    print(f"F = {lambda_2_pow_6} / ({lambda_1_pow_6} + {lambda_2_pow_6})")
    print(f"F = {lambda_2_pow_6} / {denominator}\n")
    
    print(f"Resulting Fidelity: {fidelity}")

# Run the solver
solve_quantum_fidelity()