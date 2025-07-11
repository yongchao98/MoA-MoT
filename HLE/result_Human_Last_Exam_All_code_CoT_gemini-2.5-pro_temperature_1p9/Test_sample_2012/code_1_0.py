import numpy as np

def solve_quantum_fidelity():
    """
    Solves the fidelity problem based on the rules of the hypothetical universe U.
    """
    # Define the observable operator O
    O = np.array([
        [3, 1],
        [1, 2]
    ])

    # 1. Calculate the eigenvalues of the operator O.
    # Since the matrix is Hermitian (real and symmetric), we can use eigvalsh
    # which is efficient and guarantees real eigenvalues.
    eigenvalues = np.linalg.eigvalsh(O)

    # 2. Sort eigenvalues in descending order to identify the largest (lambda_1)
    # and second-largest (lambda_2).
    eigenvalues.sort()
    lambda_2 = eigenvalues[0] # The smaller eigenvalue
    lambda_1 = eigenvalues[1] # The larger eigenvalue

    print(f"The eigenvalues are: {lambda_1} and {lambda_2}")
    
    # 3. According to the problem, the fidelity F is given by the formula:
    # F = lambda_2^6 / (lambda_1^6 + lambda_2^6)
    
    # Calculate the components of the formula
    lambda_1_pow_6 = lambda_1**6
    lambda_2_pow_6 = lambda_2**6
    denominator = lambda_1_pow_6 + lambda_2_pow_6
    
    # 4. Calculate the final fidelity
    fidelity = lambda_2_pow_6 / denominator

    # 5. Print the final equation with all numbers and the result
    print("\nThe fidelity F is calculated as lambda_2^6 / (lambda_1^6 + lambda_2^6).")
    print("\nCalculating the components:")
    print(f"Second-largest eigenvalue (lambda_2) to the power of 6 is: {lambda_2_pow_6}")
    print(f"Largest eigenvalue (lambda_1) to the power of 6 is: {lambda_1_pow_6}")
    print(f"The sum (lambda_1^6 + lambda_2^6) is: {denominator}")

    print("\nFinal equation for fidelity:")
    print(f"F = {lambda_2_pow_6} / {denominator}")
    
    print("\nThe final fidelity is:")
    print(fidelity)

solve_quantum_fidelity()