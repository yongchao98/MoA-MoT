import numpy as np

def check_physics_statement():
    """
    This function programmatically checks the four statements related to the
    Y-component of the intrinsic angular momentum operator (Ay) for a muon.
    """
    # The operator is Ay = c * S, where c = h / (4*pi) and S is the Pauli-Y matrix.
    # For calculation purposes, we can set the constant c = 1, as the core
    # properties (reality of eigenvalues, relationships between eigenvectors)
    # are independent of its specific value.
    c = 1.0
    
    # Define the S matrix (Pauli-Y)
    S = np.array([[0, -1j],
                  [1j,  0]], dtype=complex)

    # Define the operator Ay
    Ay = c * S

    # --- Step 1: Calculate eigenvalues and eigenvectors of Ay ---
    eigenvalues_Ay, eigenvectors_Ay = np.linalg.eig(Ay)
    
    print("--- Starting Verification ---")
    print(f"Operator Ay:\n{Ay}")
    print(f"Calculated eigenvalues of Ay: {np.round(eigenvalues_Ay, 3)}")
    print(f"Calculated eigenvector 1 (for eigenvalue {np.round(eigenvalues_Ay[0], 2)}): {np.round(eigenvectors_Ay[:, 0], 3)}")
    print(f"Calculated eigenvector 2 (for eigenvalue {np.round(eigenvalues_Ay[1], 2)}): {np.round(eigenvectors_Ay[:, 1], 3)}")
    print("-" * 40)

    # --- Step 2: Evaluate each statement ---
    
    # Statement C & D: Check properties of eigenvalues.
    # Ay is a Hermitian matrix (Ay.conj().T == Ay), so its eigenvalues must be real.
    are_eigenvalues_real = np.all(np.isclose(np.imag(eigenvalues_Ay), 0))
    
    if not are_eigenvalues_real:
        reason_C = "Incorrect. Calculation shows eigenvalues have non-zero imaginary parts, but they must be real for a Hermitian operator."
        reason_D = "Incorrect. Calculation shows eigenvalues have non-zero imaginary parts, but they must be real for a Hermitian operator."
    else:
        reason_C = "Incorrect. The eigenvalues of Ay are purely real. Statement C claims they have non-zero imaginary parts (+/- 1/2)."
        reason_D = "Incorrect. The eigenvalues of Ay are purely real. Statement D claims they have non-zero imaginary parts (+/- 2πh)."
    
    print(f"Analysis of C: {reason_C}")
    print(f"Analysis of D: {reason_D}")
    print("-" * 40)

    # Statement B: The eigenfunctions φ of Ay are the basis functions of the matrix operator Ay.
    # The matrix is given in the standard (Z-component) basis: {[1, 0], [0, 1]}.
    # We check if the calculated eigenvectors are the standard basis vectors.
    v1 = eigenvectors_Ay[:, 0]
    # An eigenvector is a standard basis vector if one component is 0 and the other is non-zero.
    is_standard_basis = (np.isclose(v1[0], 0) and not np.isclose(v1[1], 0)) or \
                        (not np.isclose(v1[0], 0) and np.isclose(v1[1], 0))
    
    if is_standard_basis:
        reason_B = "Incorrect. This would only be true if Ay was diagonal, which it is not."
    else:
        reason_B = "Incorrect. The operator Ay is written in the standard basis (eigenvectors of Az), but its own eigenvectors are different. The calculated eigenvectors are linear combinations of the standard basis vectors, not the basis vectors themselves."

    print(f"Analysis of B: {reason_B}")
    print("-" * 40)

    # Statement A: The eigenfunction of Ay can also be an eigenfunction of A^2, but not of Az.
    
    # Part 1: Is an eigenfunction of Ay also an eigenfunction of Ay^2?
    # If Ay|v> = a|v>, then Ay^2|v> = Ay(a|v>) = a(Ay|v>) = a^2|v>. So, yes.
    A_squared = Ay @ Ay
    v1 = eigenvectors_Ay[:, 0]
    a1 = eigenvalues_Ay[0]
    is_eigenvector_of_A_squared = np.allclose(A_squared @ v1, (a1**2) * v1)

    # Part 2: Is an eigenfunction of Ay also an eigenfunction of Az?
    # This is true only if [Ay, Az] = 0. For spin operators, [Sy, Sz] = i*hbar*Sx != 0.
    Sz = np.array([[1, 0], [0, -1]], dtype=complex)
    Az = c * Sz
    
    # Check if v1 is an eigenvector of Az by testing for linear dependence.
    # For 2D vectors u and w, they are linearly dependent if u[0]*w[1] - u[1]*w[0] = 0.
    u = Az @ v1
    w = v1
    is_eigenvector_of_Az = np.isclose(u[0]*w[1] - u[1]*w[0], 0)

    if is_eigenvector_of_A_squared and not is_eigenvector_of_Az:
        reason_A = "Correct. An eigenvector of Ay is also an eigenvector of Ay^2 (since [Ay, Ay^2]=0). However, Ay and Az do not commute ([Ay, Az] != 0), so they do not share eigenvectors. Numerical verification confirms this."
    else:
        reason_A = "Incorrect. The theoretical principles or the numerical calculation failed."

    print(f"Analysis of A: {reason_A}")
    print("-" * 40)

    # Final conclusion
    if reason_A.startswith("Correct"):
        print("\nConclusion: Statement A is the only correct option.")
    else:
        print("\nConclusion: The analysis did not find a correct option.")

# Execute the verification function
check_physics_statement()