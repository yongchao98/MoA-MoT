import numpy as np

def calculate_fidelity_in_universe_u():
    """
    Solves for the fidelity based on the rules of the fictional universe U.
    """
    # Step 1: Define the observable operator O
    O = np.array([[3, 1],
                  [1, 2]])

    # Step 2: Find the eigenvalues and eigenvectors of O.
    # np.linalg.eigh is used for Hermitian matrices and returns eigenvalues
    # in ascending order.
    eigenvalues, eigenvectors = np.linalg.eigh(O)

    # The second-largest eigenvalue is the first element of the sorted list.
    lambda_2 = eigenvalues[0]
    # The largest eigenvalue is the second element.
    lambda_1 = eigenvalues[1]

    # The target eigenstate is the one associated with lambda_2.
    # v_2 = eigenvectors[:, 0]

    # Step 3 & 4: Determine the fidelity.
    # According to Law 2, the post-measurement state |psi_final> is a superposition
    # of eigenvectors |v_i> with coefficients proportional to lambda_i^3.
    # |psi_final> = C * (lambda_1^3 * |v_1> + lambda_2^3 * |v_2>)
    # After normalization, the state is:
    # |psi_final> = (lambda_1^3 * |v_1> + lambda_2^3 * |v_2>) / sqrt(lambda_1^6 + lambda_2^6)
    #
    # The fidelity F is |<v_2 | psi_final>|^2.
    # <v_2 | psi_final> = <v_2 | (lambda_1^3*|v_1> + lambda_2^3*|v_2>) / sqrt(...)
    # Using orthonormality (<v_2|v_1>=0, <v_2|v_2>=1), this simplifies to:
    # <v_2 | psi_final> = lambda_2^3 / sqrt(lambda_1^6 + lambda_2^6)
    #
    # Squaring this gives the final formula for fidelity:
    # F = lambda_2^6 / (lambda_1^6 + lambda_2^6)

    numerator = lambda_2**6
    denominator = lambda_1**6 + lambda_2**6
    fidelity = numerator / denominator

    # Output the steps of the calculation as requested
    print("The observable operator O is:")
    print(O)
    print(f"\nThe largest eigenvalue (lambda_1) is: {lambda_1:.6f}")
    print(f"The second-largest eigenvalue (lambda_2) is: {lambda_2:.6f}")
    print("\nThe fidelity F is calculated using the formula: (lambda_2^6) / (lambda_1^6 + lambda_2^6)")
    print("\nPlugging in the numbers for the final equation:")
    print(f"F = {numerator:.6f} / ({lambda_1**6:.6f} + {numerator:.6f})")
    print(f"F = {numerator:.6f} / {denominator:.6f}")
    print(f"\nThe calculated fidelity is: {fidelity}")

calculate_fidelity_in_universe_u()
<<<0.0031619586326753177>>>