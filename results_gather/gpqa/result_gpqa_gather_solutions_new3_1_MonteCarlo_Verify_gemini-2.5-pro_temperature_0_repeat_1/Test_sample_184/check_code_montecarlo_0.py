import sympy

def check_hamiltonian_eigenvalues():
    """
    Checks the eigenvalues of the Hamiltonian H = ε * σ.n using symbolic math.
    """
    try:
        # 1. Define symbolic variables
        # ε is a positive real energy constant
        epsilon = sympy.symbols('varepsilon', real=True, positive=True)
        # n_x, n_y, n_z are real components of the unit vector n
        nx, ny, nz = sympy.symbols('n_x n_y n_z', real=True)

        # 2. Define the Pauli matrices
        sigma_x = sympy.Matrix([[0, 1], [1, 0]])
        sigma_y = sympy.Matrix([[0, -sympy.I], [sympy.I, 0]])
        sigma_z = sympy.Matrix([[1, 0], [0, -1]])

        # 3. Construct the Hamiltonian matrix H
        H = epsilon * (nx * sigma_x + ny * sigma_y + nz * sigma_z)

        # 4. Calculate the eigenvalues of H symbolically
        # The .eigenvals() method returns a dictionary of {eigenvalue: multiplicity}
        eigenvals_expr = H.eigenvals()
        
        # The eigenvalues will be expressed in terms of nx, ny, nz.
        # For example: {-varepsilon*sqrt(n_x**2 + n_y**2 + n_z**2): 1, 
        #               varepsilon*sqrt(n_x**2 + n_y**2 + n_z**2): 1}

        # 5. Apply the unit vector constraint: n_x^2 + n_y^2 + n_z^2 = 1
        unit_vector_mag_squared = nx**2 + ny**2 + nz**2
        
        simplified_eigenvals = set()
        for val_expr in eigenvals_expr.keys():
            # Substitute the magnitude squared of n with 1 and simplify
            simplified_val = sympy.simplify(val_expr.subs(unit_vector_mag_squared, 1))
            simplified_eigenvals.add(simplified_val)

        # 6. Compare with the expected answer from option A
        expected_eigenvals = {epsilon, -epsilon}

        if simplified_eigenvals == expected_eigenvals:
            return "Correct"
        else:
            return (f"Incorrect. The calculated eigenvalues after applying the unit vector "
                    f"constraint are {simplified_eigenvals}, but the expected eigenvalues "
                    f"from answer A are {expected_eigenvals}.")

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_hamiltonian_eigenvalues()
print(result)