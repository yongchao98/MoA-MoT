import numpy as np

def solve_quantum_problem():
    """
    Solves the quantum mechanics problem based on the rules of universe U.
    """
    # Define the observable operator and the initial state
    O = np.array([[3, 1],
                  [1, 2]], dtype=float)

    psi = np.array([np.sqrt(3)/2,
                    1/2], dtype=float)

    # Step 1: Find the eigenvalues and eigenvectors of O
    eigenvalues, eigenvectors = np.linalg.eigh(O)

    # np.linalg.eigh returns eigenvalues in ascending order.
    # To match the problem's "largest" and "second-largest" terminology, we sort descending.
    sort_indices = np.argsort(eigenvalues)[::-1]
    
    # Eigenvalue associated with the largest eigenvalue
    lambda_1 = eigenvalues[sort_indices[0]] 
    # Eigenvector for the largest eigenvalue
    v1 = eigenvectors[:, sort_indices[0]]

    # Eigenvalue associated with the second-largest eigenvalue
    lambda_2 = eigenvalues[sort_indices[1]]
    # Eigenvector for the second-largest eigenvalue (our target state)
    v2 = eigenvectors[:, sort_indices[1]]
    
    # Step 2: Express the initial state in the eigenbasis
    # c_i = <v_i | psi>
    c1 = np.dot(v1, psi)
    c2 = np.dot(v2, psi)

    # Step 3: Apply the measurement rule of universe U to find the new coefficients
    # c'_i = c_i * lambda_i^3
    c_prime_1 = c1 * (lambda_1 ** 3)
    c_prime_2 = c2 * (lambda_2 ** 3)
    
    # The components for the fidelity calculation are the squared magnitudes
    prob_prime_1 = np.abs(c_prime_1)**2
    prob_prime_2 = np.abs(c_prime_2)**2

    # Step 4 & 5: Calculate the fidelity
    # F = |<v2 | psi'>|^2 = |c'_2|^2 / (|c'_1|^2 + |c'_2|^2)
    fidelity = prob_prime_2 / (prob_prime_1 + prob_prime_2)
    
    print("The final fidelity is calculated using the formula:")
    print("F = |c'_2|^2 / (|c'_1|^2 + |c'_2|^2)\n")
    print(f"Where |c'_1|^2 is the squared magnitude of the new amplitude for the largest eigenstate, and |c'_2|^2 is for the second-largest.")
    print("\nCalculated numeric values for the equation:")
    print(f"|c'_1|^2 = {prob_prime_1}")
    print(f"|c'_2|^2 = {prob_prime_2}")
    
    print(f"\nFinal Equation:")
    print(f"Fidelity = {prob_prime_2} / ({prob_prime_1} + {prob_prime_2})")
    
    print(f"\nResulting Fidelity:")
    print(fidelity)
    
    # Return final answer in the specified format
    print(f"\n<<<{fidelity:.10f}>>>")

solve_quantum_problem()