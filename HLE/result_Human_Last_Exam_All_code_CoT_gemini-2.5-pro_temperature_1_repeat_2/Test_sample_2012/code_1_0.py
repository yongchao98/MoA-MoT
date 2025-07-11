import numpy as np

def solve_quantum_problem():
    """
    Calculates the fidelity based on the quantum rules of universe U.
    """
    # Define the observable operator
    O = np.array([[3, 1],
                  [1, 2]])

    # Step 1: Calculate eigenvalues and eigenvectors of the operator O.
    # The numpy linalg.eig function computes the eigenvalues and right eigenvectors.
    eigenvalues, eigenvectors = np.linalg.eig(O)

    # The eigenvalues are not guaranteed to be sorted, so we sort them in descending order.
    # In this 2D case, this separates the largest from the second-largest eigenvalue.
    sort_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[sort_indices]

    # Assign eigenvalues based on their size.
    # lambda_1 is the largest, lambda_2 is the second-largest.
    lambda_1 = sorted_eigenvalues[0]
    lambda_2 = sorted_eigenvalues[1]

    print("Step 1: The eigenvalues of the operator have been found.")
    print(f"Largest eigenvalue (lambda_1): {lambda_1:.6f}")
    print(f"Second-largest eigenvalue (lambda_2): {lambda_2:.6f}\n")

    # Step 2: Determine the formula for fidelity.
    # According to Rule 2, the post-measurement state |psi_new> is a superposition
    # of the eigenstates, with coefficients proportional to the cube of the eigenvalues.
    # |psi_new> = N * (lambda_1^3 * |v1> + lambda_2^3 * |v2>),
    # where N is the normalization constant, N = 1 / sqrt(lambda_1^6 + lambda_2^6).
    #
    # The fidelity F is calculated with respect to the eigenstate of the second-largest
    # eigenvalue, |v2>.
    # F = |<v2 | psi_new>|^2
    #
    # Due to the orthonormality of eigenvectors (<v2|v1>=0, <v2|v2>=1), this simplifies to:
    # F = lambda_2^6 / (lambda_1^6 + lambda_2^6)
    print("Step 2: The fidelity is given by the formula F = lambda_2^6 / (lambda_1^6 + lambda_2^6).\n")


    # Step 3: Calculate the numerical values for the formula.
    lambda_1_pow_6 = lambda_1**6
    lambda_2_pow_6 = lambda_2**6
    sum_lambda_pow_6 = lambda_1_pow_6 + lambda_2_pow_6

    fidelity = lambda_2_pow_6 / sum_lambda_pow_6
    
    print("Step 3: Plugging the values into the formula:")
    print(f"F = {lambda_2_pow_6:.6f} / ({lambda_1_pow_6:.6f} + {lambda_2_pow_6:.6f})")
    print(f"F = {lambda_2_pow_6:.6f} / {sum_lambda_pow_6:.6f}\n")

    print("Final Result:")
    print(f"The fidelity of the post-measurement state is: {fidelity}")

solve_quantum_problem()
<<<0.0030958611>>>