import sympy

def check_hamiltonian_eigenvalues():
    """
    Symbolically calculates the eigenvalues of the Hamiltonian H = ε * σ.n
    to verify the correctness of the given answer.
    """
    try:
        # Define the symbols used in the Hamiltonian
        # ε is a real, positive energy constant
        epsilon = sympy.Symbol('varepsilon', real=True, positive=True)
        # n_x, n_y, n_z are the real components of the unit vector n
        nx = sympy.Symbol('n_x', real=True)
        ny = sympy.Symbol('n_y', real=True)
        nz = sympy.Symbol('n_z', real=True)
        # lambda will represent the eigenvalues
        lam = sympy.Symbol('lambda')

        # Define the Pauli matrices as 2x2 sympy Matrix objects
        sigma_x = sympy.Matrix([[0, 1], [1, 0]])
        sigma_y = sympy.Matrix([[0, -sympy.I], [sympy.I, 0]])
        sigma_z = sympy.Matrix([[1, 0], [0, -1]])

        # Construct the Hamiltonian matrix H = ε * (n_x*σ_x + n_y*σ_y + n_z*σ_z)
        H = epsilon * (nx * sigma_x + ny * sigma_y + nz * sigma_z)

        # The eigenvalues are the roots of the characteristic polynomial, det(H - λI) = 0
        # First, construct the matrix (H - λI)
        identity_matrix = sympy.eye(2)
        char_matrix = H - lam * identity_matrix

        # Calculate the determinant, which is the characteristic polynomial
        char_poly = sympy.det(char_matrix)
        
        # The polynomial will contain terms like nx**2, ny**2, nz**2.
        # We must use the constraint that n is a unit vector: nx**2 + ny**2 + nz**2 = 1.
        # Sympy's expand() function helps to group these terms.
        expanded_poly = sympy.expand(char_poly)
        
        # Substitute the unit vector constraint into the polynomial
        unit_vector_sum_of_squares = nx**2 + ny**2 + nz**2
        final_poly = expanded_poly.subs(unit_vector_sum_of_squares, 1)

        # Solve the final polynomial for lambda to find the eigenvalues
        eigenvalues = sympy.solve(final_poly, lam)

        # The expected answer is A, which corresponds to eigenvalues {+ε, -ε}
        expected_eigenvalues = {epsilon, -epsilon}

        # Check if the calculated eigenvalues match the expected ones
        if set(eigenvalues) == expected_eigenvalues:
            return "Correct"
        else:
            return (f"Incorrect. The symbolic calculation yielded the wrong eigenvalues.\n"
                    f"Expected: {expected_eigenvalues}\n"
                    f"Got: {set(eigenvalues)}")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_hamiltonian_eigenvalues()
print(result)