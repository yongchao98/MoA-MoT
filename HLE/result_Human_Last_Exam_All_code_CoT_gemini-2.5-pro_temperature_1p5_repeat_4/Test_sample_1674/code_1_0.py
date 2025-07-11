import numpy as np

def run_simulation():
    """
    Simulates light passing through a depolarizing random medium and a
    birefringent plate to test if the system is invertible.
    """
    # The theory being tested is whether the effect of an optical system can be perfectly inverted.
    # We will show that when the system causes depolarization, this theory fails.

    # We represent the polarization state of light using a 4-element Stokes vector [I, Q, U, V],
    # where I is total intensity, and Q, U, V describe the polarization.
    # We start with a fully polarized horizontal input beam.
    S_in = np.array([1.0, 1.0, 0.0, 0.0])
    print("THEORY: An optical system's effects can be inverted.")
    print("TEST: Does this hold if we add a birefringent plate to a random medium?")
    print("-" * 50)
    print("1. INITIAL STATE")
    print("We start with horizontally polarized light. Its Stokes vector is:")
    print(S_in)
    print("Initial Degree of Polarization: 1.0 (Fully Polarized)")
    print("-" * 50)

    # We model the system's components using 4x4 Mueller matrices.
    # A random medium often acts as a depolarizer, scrambling polarization.
    # We model it as a perfect depolarizer, which transforms any input
    # into unpolarized light. Its Mueller matrix is singular (not invertible).
    # The equation for this matrix is:
    # M_T = [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    M_T = np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0]
    ])
    
    # The birefringent medium is a quarter-wave plate (QWP) at 45 degrees.
    # Its Mueller matrix is well-defined and invertible.
    # The equation for this matrix is:
    # M_B = [[1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0]]
    M_B = np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, -1.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 0.0]
    ])
    
    print("2. THE OPTICAL SYSTEM")
    print("The system consists of a random medium (depolarizer) followed by a birefringent plate (QWP).")
    print("\nMueller Matrix for the Random Medium (Depolarizer):")
    print(M_T)
    print("\nMueller Matrix for the Birefringent Plate (QWP):")
    print(M_B)

    # The total system is described by the product of the matrices.
    # The light first passes through T, then B.
    M_total = M_B @ M_T
    print("\nThe total Mueller matrix for the combined system (M_total = M_B * M_T) is:")
    print(M_total)
    print("-" * 50)
    
    # Calculate the output state of the light.
    # The final equation is: S_out = M_total * S_in
    S_out = M_total @ S_in
    # Calculate the Degree of Polarization (DoP) of the output light.
    # DoP = sqrt(Q^2 + U^2 + V^2) / I
    dop_out = np.sqrt(S_out[1]**2 + S_out[2]**2 + S_out[3]**2) / S_out[0]

    print("3. FINAL STATE")
    print(f"The final output Stokes vector is calculated as: S_out = M_total * S_in")
    print(f"{S_out} = \n{M_total}\n * \n{S_in}")
    print("\nFinal Degree of Polarization:", round(dop_out, 4))
    print("The fully polarized input has become completely unpolarized. Information is lost.")
    print("-" * 50)

    # Test the theory: Can we invert M_total?
    print("4. TESTING THE INVERSION")
    print("For the theory to hold, we must be able to find an inverse for M_total.")
    print("This would allow us to recover the input from the output: S_in = M_inverse * S_out.")
    print("We will now attempt to calculate the inverse of M_total...")

    try:
        M_inv = np.linalg.inv(M_total)
        print("\nSUCCESS: An inverse matrix was found:")
        print(M_inv)
    except np.linalg.LinAlgError:
        print("\nFAILURE: The inversion failed.")
        determinant = np.linalg.det(M_total)
        print(f"The total Mueller matrix is singular (its determinant is {determinant}).")
        print("This means no unique inverse exists, and the system's effect is irreversible.")
        print("\nCONCLUSION: The theory does not hold.")

if __name__ == '__main__':
    run_simulation()