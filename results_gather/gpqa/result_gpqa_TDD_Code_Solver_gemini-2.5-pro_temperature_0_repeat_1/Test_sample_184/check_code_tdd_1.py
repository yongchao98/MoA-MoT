import numpy as np

def check_correctness():
    """
    Checks the correctness of the answer for the eigenvalues of the Hamiltonian
    H = epsilon * sigma . n.

    The question states H = epsilon * vec(sigma) . vec(n), where vec(n) is an
    arbitrary unit vector. The proposed answer is that the eigenvalues are
    +epsilon and -epsilon.

    This function will:
    1. Define the Pauli matrices (sigma_x, sigma_y, sigma_z).
    2. Choose an arbitrary value for the energy constant epsilon.
    3. Define the expected eigenvalues based on the proposed answer.
    4. Test the claim for several different unit vectors vec(n), including
       both cardinal directions and a randomly generated direction.
    5. For each test case, it constructs the Hamiltonian matrix H.
    6. It then computes the eigenvalues of H using numpy.linalg.eigvals.
    7. It compares the computed eigenvalues with the expected eigenvalues.
    8. If all tests pass, it returns "Correct". Otherwise, it returns the
       reason for the failure.
    """
    # 1. Define Pauli matrices
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # 2. Choose an arbitrary value for epsilon
    epsilon = 4.2  # An arbitrary non-integer value

    # 3. Define the expected eigenvalues from answer A
    # The eigenvalues should be [-epsilon, +epsilon]
    expected_eigenvalues = np.array([-epsilon, epsilon])

    # 4. Define a list of test unit vectors
    test_vectors = [
        np.array([1., 0., 0.]),  # n along x-axis
        np.array([0., 1., 0.]),  # n along y-axis
        np.array([0., 0., 1.]),  # n along z-axis
        np.array([1., 1., 1.]) / np.sqrt(3),  # n along a diagonal
    ]
    # Add a randomly generated unit vector for a more general test
    np.random.seed(42) # for reproducibility
    random_vec = np.random.randn(3)
    random_unit_vec = random_vec / np.linalg.norm(random_vec)
    test_vectors.append(random_unit_vec)

    # 5. Loop through test vectors
    for n_vec in test_vectors:
        # 5a. Construct the Hamiltonian matrix H
        H = epsilon * (n_vec[0] * sigma_x + n_vec[1] * sigma_y + n_vec[2] * sigma_z)

        # 6. Compute the eigenvalues
        calculated_eigenvalues = np.linalg.eigvals(H)

        # Sort eigenvalues for consistent comparison
        calculated_eigenvalues.sort()

        # 7. Compare computed eigenvalues with expected ones
        if not np.allclose(calculated_eigenvalues, expected_eigenvalues):
            return (f"Incorrect. The answer claims the eigenvalues are always "
                    f"[{-epsilon}, {epsilon}]. However, for the unit vector "
                    f"n = {np.round(n_vec, 3)}, the calculated eigenvalues are "
                    f"{np.round(calculated_eigenvalues, 5)}.")

    # 8. If all tests pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)