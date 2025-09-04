import numpy as np
from scipy.constants import h, pi

def check_answer():
    """
    Checks the correctness of the provided answer for the quantum mechanics problem.
    """
    # Use complex numbers for matrices
    i = 1j

    # Define constants
    # h is Planck's constant from scipy.constants
    hbar = h / (2 * pi)
    c = h / (4 * pi) # as given in the problem, which is hbar/2

    # Define the Pauli matrices and the Identity matrix
    S_y = np.array([[0, -i],
                    [i,  0]], dtype=complex)
    S_x = np.array([[0, 1],
                    [1, 0]], dtype=complex)
    S_z = np.array([[1, 0],
                    [0,-1]], dtype=complex)
    I = np.identity(2, dtype=complex)

    # Define the angular momentum operators
    Ay = c * S_y
    # We need Az and A^2 to check statement B
    Az = (hbar / 2) * S_z
    # For a spin-1/2 particle, A^2 = s(s+1)ħ^2 * I, where s=1/2
    # A^2 = (1/2)*(3/2)*ħ^2 * I = (3/4)*ħ^2 * I
    A_squared = (3/4) * (hbar**2) * I

    # --- Step 1: Calculate eigenvalues of Ay to check statements A and C ---
    eigenvalues, eigenvectors = np.linalg.eig(Ay)
    
    # Expected eigenvalues are +/- hbar/2 = +/- h/(4*pi)
    expected_eig_val_1 = h / (4 * pi)
    expected_eig_val_2 = -h / (4 * pi)

    # Check Statement C: "The imaginary part of the eigenvalue of Ay are +2πh or –2πh, and the real part of that are +h/4π or –h/4π."
    # Check Statement A: "The imaginary part of the eigenvalue of Ay are +1/2 or –1/2, and the real part of that are +1 or –1."
    for val in eigenvalues:
        if not np.isclose(np.imag(val), 0):
            return f"Incorrect. The final answer states that statements A and C are incorrect because the eigenvalues are purely real. However, the calculated eigenvalue {val} has a non-zero imaginary part."
    
    # Check if the real parts match the expected values
    if not ( (np.isclose(eigenvalues[0], expected_eig_val_1) and np.isclose(eigenvalues[1], expected_eig_val_2)) or \
             (np.isclose(eigenvalues[0], expected_eig_val_2) and np.isclose(eigenvalues[1], expected_eig_val_1)) ):
        return f"Incorrect. The calculated eigenvalues {eigenvalues} do not match the expected values +/- h/(4*pi) ({expected_eig_val_1})."

    # --- Step 2: Check statement D ---
    # "The eigenfunctions φ of the operator Ay are the basis functions of the matrix operator Ay given above."
    # The basis functions are the standard basis vectors e1=[1,0] and e2=[0,1].
    # An eigenvector is parallel to a basis vector if their dot product's magnitude is 1 (for normalized vectors).
    # A simpler check: an eigenvector is a basis vector if one of its components is zero and the other is non-zero.
    e1 = np.array([1, 0])
    e2 = np.array([0, 1])
    v1 = eigenvectors[:, 0]
    v2 = eigenvectors[:, 1]

    is_v1_basis = np.allclose(np.abs(v1), e1) or np.allclose(np.abs(v1), e2)
    is_v2_basis = np.allclose(np.abs(v2), e1) or np.allclose(np.abs(v2), e2)

    if is_v1_basis or is_v2_basis:
        return f"Incorrect. The final answer states that statement D is incorrect. However, a calculated eigenvector of Ay ({v1} or {v2}) is a standard basis vector, which would make statement D correct."

    # --- Step 3: Check statement B ---
    # "The eigenfunction of the operator Ay can also be an eigenfunction of A^2, but not of the Z-component, Az."
    
    # Part 1: Check if Ay and A^2 share eigenfunctions. This is true if they commute.
    commutator_Ay_A2 = Ay @ A_squared - A_squared @ Ay
    if not np.allclose(commutator_Ay_A2, np.zeros((2,2))):
        return f"Incorrect. The final answer relies on [Ay, A^2] = 0. However, the calculated commutator is non-zero:\n{commutator_Ay_A2}"

    # Part 2: Check if Ay and Az share eigenfunctions. This is false if they do not commute.
    commutator_Ay_Az = Ay @ Az - Az @ Ay
    if np.allclose(commutator_Ay_Az, np.zeros((2,2))):
        return f"Incorrect. The final answer relies on [Ay, Az] != 0. However, the calculated commutator is zero, meaning they share eigenfunctions."

    # Direct check: Test if an eigenvector of Ay is also an eigenvector of Az.
    # (A @ v) should be parallel to v. We check if the cross product of the complex vectors is zero.
    # A simpler way is to check if (A @ v) / v gives a constant value for all non-zero elements of v.
    v_ay = eigenvectors[:, 0] # Take one eigenvector of Ay
    
    # Test with Az
    v_after_az = Az @ v_ay
    # Check for collinearity. If v_after_az = k * v_ay, then they are eigenvectors.
    # This is true if v_after_az[0]*v_ay[1] - v_after_az[1]*v_ay[0] == 0
    if np.isclose(v_after_az[0]*v_ay[1] - v_after_az[1]*v_ay[0], 0):
        return f"Incorrect. The final answer states an eigenfunction of Ay is NOT an eigenfunction of Az. However, the eigenvector {v_ay} of Ay was found to also be an eigenvector of Az."

    # Test with A_squared
    v_after_a2 = A_squared @ v_ay
    if not np.isclose(v_after_a2[0]*v_ay[1] - v_after_a2[1]*v_ay[0], 0):
        return f"Incorrect. The final answer states an eigenfunction of Ay IS an eigenfunction of A^2. However, the eigenvector {v_ay} of Ay was found to NOT be an eigenvector of A^2."

    # If all checks pass, the reasoning for selecting B is sound.
    return "Correct"

# Run the check
result = check_answer()
print(result)