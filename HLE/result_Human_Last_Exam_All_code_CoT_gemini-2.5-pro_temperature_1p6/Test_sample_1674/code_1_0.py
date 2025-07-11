import numpy as np

def simulate_system(B, B_name):
    """
    Simulates the optical system and tests the inversion theory.
    For simplicity, we ignore the free-space propagation operators as they
    are invertible and do not affect the core logic.

    Args:
        B (np.ndarray): The 2x2 Jones matrix for the added medium.
        B_name (str): The name of the medium for printing.
    """
    print(f"--- Testing with an added '{B_name}' ---")

    # 1. Define the system components
    # Input beam with diagonal polarization (an equal mix of horizontal and vertical components)
    input_jones_vector = np.array([[1],
                                   [1]]) / np.sqrt(2)

    # A sample invertible transmission matrix T for the random medium
    T = np.array([[0.8, 0.1j],
                  [0.2j, 0.7]])

    # 2. Forward Propagation
    # The output is calculated by applying the matrices in order: B, then T.
    # The full equation would be Output = P2 * T * B * P1 * Input. We simplify to:
    # Output = T * B * Input
    output_jones_vector = T @ B @ input_jones_vector

    print("\nInitial Values:")
    print("---------------")
    print(f"Original Input Vector:\n{input_jones_vector}")
    print(f"\nRandom Medium Matrix (T):\n{T}")
    print(f"\nAdded Medium Matrix (B) - {B_name}:\n{B}")

    print("\nForward Propagation Result:")
    print("---------------------------")
    print(f"Final Output Vector (T * B * Input):\n{output_jones_vector}")

    # 3. Test the Inversion Theory
    # The theory states: Recovered_Input = B_inverse * T_inverse * Output
    print("\nTesting Inversion Theory:")
    print("-------------------------")
    print("Attempting to recover the input using the equation: Recovered_Input = B_inv @ T_inv @ Output")

    try:
        # Calculate the inverses of the matrices
        T_inv = np.linalg.inv(T)
        print(f"\nInverse of T (T_inv) found:\n{T_inv}")
        
        B_inv = np.linalg.inv(B)
        print(f"\nInverse of B (B_inv) found:\n{B_inv}")

        # Apply the inverse operations to the output to recover the input
        recovered_input = B_inv @ T_inv @ output_jones_vector

        print(f"\nFinal Recovered Input Vector:\n{recovered_input}")

        # Check if the recovered input is numerically close to the original input
        if np.allclose(input_jones_vector, recovered_input):
            print("\nConclusion: SUCCESS! The recovered input matches the original input.")
            print("The theory holds true for this case.")
        else:
            print("\nConclusion: FAILURE! The recovered input does not match the original.")
            print("The theory failed, likely due to numerical precision issues.")

    except np.linalg.LinAlgError:
        print(f"\nInverse of T (T_inv) found:\n{np.linalg.inv(T)}")
        print("\nERROR: Could not find the inverse of B.")
        print("\nConclusion: FAILURE! The matrix for the added medium is singular (non-invertible).")
        print("This means information was permanently lost (one polarization was deleted).")
        print("The theory does NOT hold for this case because the system is not reversible.")
    print("-" * 50 + "\n")


# --- Case A: Birefringent Medium is an INVERTIBLE Wave Plate ---
# This is a quarter-wave plate. It is a unitary and thus invertible matrix.
quarter_wave_plate = np.array([[1, 0],
                               [0, 1j]])
simulate_system(quarter_wave_plate, "Invertible Quarter-Wave Plate")


# --- Case B: Birefringent Medium is a NON-INVERTIBLE Polarizer ---
# An ideal horizontal polarizer, which has extreme polarization-dependent loss (dichroism).
# This matrix is singular (its determinant is 0) and cannot be inverted.
horizontal_polarizer = np.array([[1, 0],
                                 [0, 0]])
simulate_system(horizontal_polarizer, "Non-Invertible Horizontal Polarizer")