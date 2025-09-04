import numpy as np

def check_correctness():
    """
    This function checks the correctness of the given answer by numerically calculating
    the eigenvalues of the Hamiltonian H = epsilon * (sigma . n) for random unit vectors n.
    
    The question asks for the eigenvalues of H.
    The provided answer is that the eigenvalues are +epsilon and -epsilon.
    """
    
    # Define the Pauli matrices
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Choose an arbitrary non-zero value for epsilon (a constant with dimensions of energy)
    epsilon = 5.0
    
    # The answer to check states the eigenvalues are +epsilon and -epsilon.
    # We sort them for consistent comparison later.
    expected_eigenvalues = np.array([-epsilon, epsilon])

    # The property must hold for any unit vector n.
    # We test this for several randomly generated unit vectors to ensure robustness.
    num_tests = 25
    for _ in range(num_tests):
        # Generate a random 3D vector and normalize it to create a unit vector n.
        # This ensures we test arbitrary directions in space.
        n_random = np.random.randn(3)
        n_unit = n_random / np.linalg.norm(n_random)
        nx, ny, nz = n_unit

        # Verify the constraint that n is a unit vector.
        if not np.isclose(nx**2 + ny**2 + nz**2, 1.0):
            return (f"Constraint check failed: The generated vector n={n_unit} "
                    f"is not a unit vector. Its squared norm is {nx**2 + ny**2 + nz**2}.")

        # Construct the Hamiltonian matrix H = epsilon * (nx*sigma_x + ny*sigma_y + nz*sigma_z)
        H = epsilon * (nx * sigma_x + ny * sigma_y + nz * sigma_z)

        # Calculate the eigenvalues of the constructed Hamiltonian matrix
        calculated_eigenvalues = np.linalg.eigvals(H)
        
        # Sort the calculated eigenvalues for a consistent comparison
        calculated_eigenvalues.sort()

        # Check if the calculated eigenvalues are numerically close to the expected ones.
        if not np.allclose(calculated_eigenvalues, expected_eigenvalues):
            return (f"Incorrect. For the unit vector n = {n_unit}, the calculated "
                    f"eigenvalues are {np.round(calculated_eigenvalues, 5)}. "
                    f"The expected eigenvalues based on the answer are {expected_eigenvalues}. "
                    "The calculated eigenvalues do not match the expected ones.")

    # If all tests with random unit vectors pass, the answer is correct.
    # This numerical check confirms the analytical result that lambda = +/- epsilon.
    return "Correct"

# Run the check
result = check_correctness()
if result == "Correct":
    # The answer is correct, so we return the option letter.
    print("Correct")
else:
    # The answer is incorrect, so we return the reason.
    print(result)