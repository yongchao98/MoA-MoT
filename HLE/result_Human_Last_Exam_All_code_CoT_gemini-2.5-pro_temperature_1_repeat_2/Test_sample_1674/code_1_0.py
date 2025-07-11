import numpy as np

def solve():
    """
    Simulates an optical system to demonstrate how adding a new element
    breaks a previously working reversal (inversion) process.
    """
    # Set a seed for reproducibility of the random matrix
    np.random.seed(0)

    # --- 1. Define the initial state and operators ---

    # We represent the beam's polarization as a 2D vector (Jones vector).
    # Let's define an initial beam with horizontal polarization.
    initial_beam = np.array([1, 0], dtype=complex)

    # Create a model for the "random medium" T.
    # This is a random 2x2 complex matrix. We ensure it's invertible.
    T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)
    while np.linalg.det(T) == 0:
        T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)

    # Create a model for the birefringent medium B.
    # We use the Jones matrix for a Quarter-Wave Plate at a 45-degree angle.
    B = 0.5 * np.array([[1 - 1j, 1 + 1j],
                        [1 + 1j, 1 - 1j]], dtype=complex)

    # Calculate the inverse of T, which is used in the reversal process.
    T_inverse = np.linalg.inv(T)

    # --- 2. Simulate the process ---

    # In the modified system, the beam passes through T, then B.
    # The combined operator is B @ T.
    # (@ is the symbol for matrix multiplication in numpy)
    output_from_new_system = B @ T @ initial_beam

    # Now, we apply the original, but now incomplete, reversal process.
    # This process only knows about T, not B.
    final_recovered_beam = T_inverse @ output_from_new_system

    # --- 3. Print the results ---

    print("This simulation demonstrates that the theory fails after adding a birefringent medium.\n")
    print("--- The Final Equation ---")
    print("The goal is to see if the 'Final Recovered Beam' equals the 'Initial Beam'.")
    print("The calculation performed is:")
    print("Final Recovered Beam = (T_inverse) * B * T * (Initial Beam)\n")

    print(f"Initial Beam vector:")
    print(f"[{initial_beam[0]:.1f}, {initial_beam[1]:.1f}]\n")

    print("Final Recovered Beam vector (after incorrect reversal):")
    # We use np.round to make the output cleaner
    final_vector_rounded = np.round(final_recovered_beam, 4)
    print(f"[{final_vector_rounded[0]}, {final_vector_rounded[1]}]\n")

    print("--- Conclusion ---")
    print("The Initial Beam and the Final Recovered Beam do not match.")
    print("This is because the reversal operator (T_inverse) was not updated to account")
    print("for the newly added birefringent medium (B). Therefore, the original theory no longer holds.")

solve()