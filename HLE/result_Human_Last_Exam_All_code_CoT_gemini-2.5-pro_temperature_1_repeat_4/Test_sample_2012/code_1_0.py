import numpy as np

def solve_quantum_fidelity():
    """
    Solves the quantum mechanics problem based on the specified rules of universe U.
    """
    # Step 1: Define the system properties in universe U
    # The observable operator
    O = np.array([[3.0, 1.0],
                  [1.0, 2.0]])

    # The initial state of the system
    psi = np.array([np.sqrt(3)/2.0, 1.0/2.0])

    # Step 2: Find the eigenvalues and eigenvectors of the observable O.
    # numpy.linalg.eigh is ideal for Hermitian matrices. It returns eigenvalues
    # sorted in ascending order and their corresponding eigenvectors as columns.
    eigenvalues, eigenvectors = np.linalg.eigh(O)

    # Assign eigenvalues and eigenvectors based on the problem's terminology.
    # lambda_2 is the second-largest (i.e., the smallest) eigenvalue.
    lambda_2 = eigenvalues[0]
    # lambda_1 is the largest eigenvalue.
    lambda_1 = eigenvalues[1]

    # v2 is the eigenvector for the second-largest eigenvalue lambda_2.
    v2 = eigenvectors[:, 0]
    # v1 is the eigenvector for the largest eigenvalue lambda_1.
    v1 = eigenvectors[:, 1]

    # Step 3: Decompose the initial state into the eigenbasis of O.
    # The coefficients are found by taking the inner product of the eigenvectors
    # with the initial state vector. For real vectors, this is the dot product.
    c1 = np.dot(v1, psi)
    c2 = np.dot(v2, psi)

    # The squared coefficients represent the probabilities in a standard measurement.
    # Let's call them prob1 and prob2.
    prob1 = c1**2
    prob2 = c2**2

    # Step 4: Calculate the components needed for the fidelity formula.
    # The fidelity F of the post-measurement state with respect to |v2> is given by:
    # F = (c2^2 * |lambda_2|^6) / (c1^2 * |lambda_1|^6 + c2^2 * |lambda_2|^6)
    # Since the eigenvalues are positive, the absolute value is redundant.
    
    lambda_1_pow6 = lambda_1**6
    lambda_2_pow6 = lambda_2**6

    # Calculate the numerator of the fidelity equation
    numerator = prob2 * lambda_2_pow6

    # Calculate the denominator of the fidelity equation
    denominator = prob1 * lambda_1_pow6 + prob2 * lambda_2_pow6

    # Calculate the final fidelity
    fidelity = numerator / denominator

    # Step 5: Print the results, showing the numbers in the final equation as requested.
    print("The fidelity F is calculated using the formula:")
    print("F = (c2^2 * lambda_2^6) / (c1^2 * lambda_1^6 + c2^2 * lambda_2^6)\n")
    print("Calculated values for the equation's components:")
    print(f"c1^2 (|<v1|psi>|^2) = {prob1}")
    print(f"lambda_1^6 = {lambda_1_pow6}")
    print(f"c2^2 (|<v2|psi>|^2) = {prob2}")
    print(f"lambda_2^6 = {lambda_2_pow6}\n")
    print("Substituting these values into the formula:")
    print(f"F = ({prob2} * {lambda_2_pow6}) / (({prob1} * {lambda_1_pow6}) + ({prob2} * {lambda_2_pow6}))")
    print(f"F = {numerator} / {denominator}")
    print(f"\nThe final fidelity is: {fidelity}")
    print(f"<<<{fidelity}>>>")

if __name__ == '__main__':
    solve_quantum_fidelity()