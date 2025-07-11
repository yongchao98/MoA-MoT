import numpy as np

def solve_quantum_fidelity():
    """
    Solves the quantum fidelity problem based on the rules of universe U.
    """
    # Step 1: Define the operator and the initial state
    O = np.array([[3, 1],
                  [1, 2]], dtype=float)
    psi = np.array([np.sqrt(3)/2, 1/2], dtype=float)

    # Step 2: Find the eigenvalues and eigenvectors of O
    # numpy.linalg.eigh returns eigenvalues in ascending order and corresponding eigenvectors as columns.
    eigenvalues, eigenvectors = np.linalg.eigh(O)

    # For a 2D system, the second-largest eigenvalue is the smaller one.
    lambda_2 = eigenvalues[0]
    lambda_1 = eigenvalues[1]

    # The eigenvectors are the columns of the returned matrix
    # |v2> corresponds to lambda_2, and |v1> corresponds to lambda_1
    v_2 = eigenvectors[:, 0]
    v_1 = eigenvectors[:, 1]

    # Step 3: Project the initial state onto the eigenbasis to find coefficients c1 and c2
    # c = <v|psi>. Since eigenvectors are real, this is a simple dot product.
    c1 = np.dot(v_1, psi)
    c2 = np.dot(v_2, psi)

    # Step 4 & 5: Calculate the final fidelity using the derived formula:
    # F = (|c2|^2 * lambda_2^6) / (|c1|^2 * lambda_1^6 + |c2|^2 * lambda_2^6)
    
    # Calculate the components of the formula
    c1_sq = c1**2
    c2_sq = c2**2
    lambda1_pow6 = lambda_1**6
    lambda2_pow6 = lambda_2**6
    
    numerator = c2_sq * lambda2_pow6
    denominator = c1_sq * lambda1_pow6 + c2_sq * lambda2_pow6
    
    fidelity = numerator / denominator

    # Print all the numbers involved in the final equation
    print("This script calculates the fidelity of the post-measurement state with the eigenstate of the second-largest eigenvalue.")
    print("\nThe final fidelity F is calculated using the formula:")
    print("F = (c2^2 * λ2^6) / (c1^2 * λ1^6 + c2^2 * λ2^6)\n")
    print("Here are the values for each component of the equation:")
    print(f"Second-largest eigenvalue, λ2 = {lambda_2}")
    print(f"Largest eigenvalue, λ1 = {lambda_1}")
    print(f"Coefficient c1 = <v1|ψ> = {c1}")
    print(f"Coefficient c2 = <v2|ψ> = {c2}")
    print("-" * 30)
    print("Values for the final equation:")
    print(f"c1^2 = {c1_sq}")
    print(f"c2^2 = {c2_sq}")
    print(f"λ1^6 = {lambda1_pow6}")
    print(f"λ2^6 = {lambda2_pow6}")
    print("-" * 30)
    print(f"Numerator = ({c2_sq}) * ({lambda2_pow6}) = {numerator}")
    print(f"Denominator = ({c1_sq}) * ({lambda1_pow6}) + ({c2_sq}) * ({lambda2_pow6}) = {denominator}")
    print("-" * 30)
    print(f"Final Fidelity F = {numerator} / {denominator} = {fidelity}")
    
    # Final answer in the required format
    print(f"\n<<<{fidelity}>>>")

# Execute the function
solve_quantum_fidelity()