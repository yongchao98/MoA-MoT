import sympy

def check_correctness():
    """
    Symbolically calculates the eigenvalues of the Hamiltonian H = epsilon * sigma . n
    to verify the correctness of the given answer.
    """
    # 1. Define symbols for the variables in the Hamiltonian.
    # epsilon is a real, positive energy constant.
    # nx, ny, nz are the real components of the unit vector n.
    epsilon = sympy.Symbol('varepsilon', real=True, positive=True)
    nx = sympy.Symbol('n_x', real=True)
    ny = sympy.Symbol('n_y', real=True)
    nz = sympy.Symbol('n_z', real=True)

    # 2. Define the Pauli matrices using sympy's Matrix object.
    sigma_x = sympy.Matrix([[0, 1], [1, 0]])
    sigma_y = sympy.Matrix([[0, -sympy.I], [sympy.I, 0]])
    sigma_z = sympy.Matrix([[1, 0], [0, -1]])

    # 3. Construct the Hamiltonian matrix H in symbolic form.
    # H = epsilon * (nx*sigma_x + ny*sigma_y + nz*sigma_z)
    H = epsilon * sympy.Matrix([
        [nz, nx - sympy.I * ny],
        [nx + sympy.I * ny, -nz]
    ])

    # 4. Calculate the eigenvalues of H symbolically.
    # The .eigenvals() method returns a dictionary of {eigenvalue: algebraic_multiplicity}.
    try:
        eigenvals_symbolic_map = H.eigenvals()
        eigenvals_symbolic_list = list(eigenvals_symbolic_map.keys())
    except Exception as e:
        return f"Failed to calculate eigenvalues symbolically. Error: {e}"

    # 5. Apply the unit vector constraint to simplify the result.
    # The constraint is n_x^2 + n_y^2 + n_z^2 = 1.
    # The symbolic eigenvalues will contain the term sqrt(nx**2 + ny**2 + nz**2).
    unit_vector_norm = sympy.sqrt(nx**2 + ny**2 + nz**2)
    
    # Substitute the norm of the unit vector (which is 1) into the eigenvalue expressions.
    simplified_eigenvals = {val.subs(unit_vector_norm, 1) for val in eigenvals_symbolic_list}

    # 6. Compare the derived eigenvalues with the proposed answer (B).
    # Answer B corresponds to eigenvalues {+epsilon, -epsilon}.
    expected_eigenvals = {epsilon, -epsilon}

    if simplified_eigenvals == expected_eigenvals:
        return "Correct"
    else:
        return (f"Incorrect. The symbolic calculation shows the eigenvalues are {simplified_eigenvals}, "
                f"which does not match the expected answer B: {expected_eigenvals}.")

# Run the check and print the result.
result = check_correctness()
print(result)