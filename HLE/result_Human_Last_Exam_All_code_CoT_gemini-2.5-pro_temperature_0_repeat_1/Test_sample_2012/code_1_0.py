import numpy as np

def solve_quantum_fidelity():
    """
    Solves the quantum mechanics problem as described for the hypothetical universe U.
    """
    # Define the observable operator O and the initial state psi
    O = np.array([[3, 1],
                  [1, 2]])

    psi_initial = np.array([np.sqrt(3)/2, 1/2])

    # Step 1: Find eigenvalues and eigenvectors of O
    eigenvalues, eigenvectors = np.linalg.eig(O)

    # Step 2: Sort eigenvalues and eigenvectors in descending order
    # The eigenvalues are real, so we can sort them directly.
    sorted_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]

    # Largest and second-largest eigenvalues and their eigenvectors
    lambda1 = sorted_eigenvalues[0]
    lambda2 = sorted_eigenvalues[1]
    v1 = sorted_eigenvectors[:, 0]
    v2 = sorted_eigenvectors[:, 1]

    # Step 3: The target state is the eigenstate of the second-largest eigenvalue
    psi_target = v2

    # Step 4: Decompose the initial state in the eigenbasis to find coefficients c1 and c2
    # c_i = <v_i | psi>
    c1 = np.vdot(v1, psi_initial)
    c2 = np.vdot(v2, psi_initial)

    # Step 5: Calculate the terms for the fidelity formula
    # F = |<psi_target | psi_post>|^2
    # psi_post is proportional to c1*lambda1^3*v1 + c2*lambda2^3*v2
    # F = (|c2|^2 * lambda2^6) / (|c1|^2 * lambda1^6 + |c2|^2 * lambda2^6)
    
    c1_sq_mag = np.abs(c1)**2
    c2_sq_mag = np.abs(c2)**2
    
    lambda1_pow6 = lambda1**6
    lambda2_pow6 = lambda2**6
    
    numerator = c2_sq_mag * lambda2_pow6
    denominator = (c1_sq_mag * lambda1_pow6) + (c2_sq_mag * lambda2_pow6)
    
    fidelity = numerator / denominator

    # Step 6: Print the results and the final equation
    print("This script calculates the fidelity based on the rules of universe U.")
    print("\n--- Calculated Values ---")
    print(f"Largest eigenvalue (lambda1): {lambda1}")
    print(f"Second-largest eigenvalue (lambda2): {lambda2}")
    print(f"Coefficient c1 = <v1|psi>: {c1}")
    print(f"Coefficient c2 = <v2|psi>: {c2}")
    
    print("\n--- Fidelity Calculation ---")
    print("The formula for fidelity (F) is: (|c2|^2 * lambda2^6) / (|c1|^2 * lambda1^6 + |c2|^2 * lambda2^6)")
    print("\nSubstituting the values:")
    print(f"F = ({c2_sq_mag:.6f} * {lambda2_pow6:.6f}) / (({c1_sq_mag:.6f} * {lambda1_pow6:.6f}) + ({c2_sq_mag:.6f} * {lambda2_pow6:.6f}))")
    print(f"F = {numerator:.6f} / ({c1_sq_mag * lambda1_pow6:.6f} + {numerator:.6f})")
    print(f"F = {numerator:.6f} / {denominator:.6f}")
    
    print("\n--- Final Answer ---")
    print(f"The final fidelity is: {fidelity}")
    
    # Output the final answer in the required format
    print(f"\n<<<{fidelity}>>>")

solve_quantum_fidelity()