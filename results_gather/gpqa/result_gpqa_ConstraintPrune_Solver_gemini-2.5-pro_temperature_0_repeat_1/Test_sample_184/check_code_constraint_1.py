import numpy as np

def check_answer():
    """
    This function checks the correctness of the given answer by numerically
    calculating the eigenvalues of the Hamiltonian H = epsilon * (sigma . n)
    for various unit vectors n.
    """
    # Let's choose an arbitrary value for the energy constant epsilon.
    # The result should be independent of this choice.
    epsilon = 5.0

    # The answer claims the eigenvalues are +epsilon and -epsilon.
    # We sort them for consistent comparison later.
    expected_eigenvalues = np.array([-epsilon, epsilon])

    # Define the Pauli matrices
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # We will test for a variety of unit vectors n to ensure the
    # result is general.
    test_vectors_n = [
        np.array([1, 0, 0]),  # n along x-axis
        np.array([0, 1, 0]),  # n along y-axis
        np.array([0, 0, 1]),  # n along z-axis
        np.array([1/np.sqrt(2), 1/np.sqrt(2), 0]), # n in xy-plane
        np.array([1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]), # n in a general direction
        np.array([0.6, -0.8, 0]), # Another unit vector
    ]
    
    # Add a randomly generated unit vector for a non-trivial test case
    random_vec = np.random.rand(3)
    random_unit_vec = random_vec / np.linalg.norm(random_vec)
    test_vectors_n.append(random_unit_vec)

    for n_vec in test_vectors_n:
        # Unpack the components of the unit vector
        nx, ny, nz = n_vec
        
        # Construct the Hamiltonian matrix
        H = epsilon * (nx * sigma_x + ny * sigma_y + nz * sigma_z)
        
        # Calculate the eigenvalues of the Hamiltonian matrix
        calculated_eigenvalues = np.linalg.eigvals(H)
        
        # Sort the calculated eigenvalues to ensure a consistent order for comparison
        calculated_eigenvalues.sort()
        
        # Check if the calculated eigenvalues are close to the expected eigenvalues.
        # We use np.allclose to handle potential floating-point inaccuracies.
        if not np.allclose(calculated_eigenvalues, expected_eigenvalues):
            return (f"Incorrect. For the unit vector n = {n_vec.round(3)}, the calculated "
                    f"eigenvalues are {calculated_eigenvalues.round(5)}, but they "
                    f"should be {expected_eigenvalues}.")

    # If all test cases pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)