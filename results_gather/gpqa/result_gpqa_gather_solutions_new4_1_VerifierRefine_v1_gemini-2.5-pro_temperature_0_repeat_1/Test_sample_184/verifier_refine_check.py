import sympy

def check_correctness():
    """
    This function verifies the eigenvalues of the Hamiltonian operator H = ε * σ.n
    using symbolic mathematics with the sympy library.

    The process is as follows:
    1. Define all components (ε, n_x, n_y, n_z, σ_x, σ_y, σ_z) as symbolic objects.
    2. Construct the Hamiltonian matrix H symbolically.
    3. Calculate the eigenvalues of the symbolic matrix H.
    4. Apply the constraint that n is a unit vector (n_x^2 + n_y^2 + n_z^2 = 1) to simplify the eigenvalue expressions.
    5. Compare the simplified, calculated eigenvalues with the expected answer {+ε, -ε} from option C.
    6. Return "Correct" if they match, or a detailed reason for the mismatch.
    """
    try:
        # 1. Define symbolic variables
        # ε (epsilon) is the energy constant, defined as a real, positive symbol.
        epsilon = sympy.Symbol('varepsilon', real=True, positive=True)
        # n_x, n_y, n_z are the real components of the unit vector n.
        nx, ny, nz = sympy.symbols('n_x n_y n_z', real=True)

        # 2. Define the Pauli matrices as sympy Matrix objects.
        sigma_x = sympy.Matrix([[0, 1], [1, 0]])
        sigma_y = sympy.Matrix([[0, -sympy.I], [sympy.I, 0]])
        sigma_z = sympy.Matrix([[1, 0], [0, -1]])

        # 3. Construct the Hamiltonian matrix H = ε * (σ_x*n_x + σ_y*n_y + σ_z*n_z)
        H = epsilon * (nx * sigma_x + ny * sigma_y + nz * sigma_z)

        # 4. Calculate the eigenvalues of the symbolic matrix H.
        # The .eigenvals() method returns a dictionary of {eigenvalue: multiplicity}.
        eigenvals_expressions = H.eigenvals()
        
        # We are interested in the unique eigenvalues, which are the keys of the dictionary.
        calculated_eigenvals = set(eigenvals_expressions.keys())

        # 5. Apply the unit vector constraint: n_x^2 + n_y^2 + n_z^2 = 1.
        # The eigenvalue expressions from sympy will contain sqrt(nx**2 + ny**2 + nz**2).
        # We substitute this term with 1 to simplify.
        norm_squared = nx**2 + ny**2 + nz**2
        simplified_eigenvals = {ev.subs(norm_squared, 1).simplify() for ev in calculated_eigenvals}

        # 6. Define the expected eigenvalues based on the provided answer C: {+ε, -ε}.
        expected_eigenvals = {epsilon, -epsilon}

        # 7. Compare the calculated set of eigenvalues with the expected set.
        if simplified_eigenvals == expected_eigenvals:
            return "Correct"
        else:
            # If they don't match, provide a detailed failure reason.
            reason = (
                f"The calculated eigenvalues do not match the expected answer.\n"
                f"The question asks for the eigenvalues of H = ε * σ.n.\n"
                f"The provided answer 'C' states the eigenvalues are {{+ε, -ε}}.\n"
                f"Symbolic calculation resulted in eigenvalues: {simplified_eigenvals}.\n"
                f"Expected eigenvalues: {expected_eigenvals}.\n"
                f"The sets do not match."
            )
            return reason

    except Exception as e:
        # Catch any potential errors during the symbolic calculation.
        return f"An error occurred during the symbolic verification: {e}"

# Execute the check and print the result.
result = check_correctness()
print(result)