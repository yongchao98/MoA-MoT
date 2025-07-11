import numpy as np

def demonstrate_optical_system_invertibility():
    """
    Demonstrates how adding different types of birefringent media can affect
    the invertibility of an optical system's transfer matrix.
    """
    print("### Optical System Invertibility Demonstration ###\n")

    # --- Step 1: Define a random, invertible transmission matrix T ---
    # We create a random 2x2 complex matrix for the "random medium".
    # We ensure it's invertible by checking if its determinant is non-zero.
    T = np.zeros((2, 2), dtype=complex)
    while np.linalg.det(T) == 0:
        T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)

    print("Random Invertible Transmission Matrix (T):")
    print(T)
    print(f"Determinant of T: {np.linalg.det(T):.4f}\n")

    # --- Step 2: Case 1 - Add an Invertible Birefringent Medium (Wave Plate) ---
    print("--- Case 1: Adding an Invertible Birefringent Medium (Quarter-Wave Plate) ---")
    # Jones matrix for a quarter-wave plate with fast axis at 45 degrees.
    # This matrix is unitary and therefore invertible.
    B_invertible = np.exp(-1j * np.pi/4) * np.array([[1, -1j], [-1j, 1]]) / np.sqrt(2)
    print("Jones Matrix for Quarter-Wave Plate (B_invertible):")
    print(B_invertible)

    # The new system matrix is the product of the elements.
    S_new_invertible = B_invertible @ T
    print("\nNew System Matrix (S_new = B_invertible * T):")
    print(S_new_invertible)

    # Check if this new system is invertible by calculating its determinant.
    det_S1 = np.linalg.det(S_new_invertible)
    print(f"Determinant of S_new: {det_S1:.4f}")

    if abs(det_S1) > 1e-9:
        print("Result: The determinant is non-zero. The system IS INVERTIBLE.")
        print("The theory holds. We can find a unique input for a given output.\n")
    else:
        print("Result: The determinant is zero. The system IS NOT INVERTIBLE.\n")


    # --- Step 3: Case 2 - Add a Non-Invertible Birefringent Medium (Polarizer) ---
    print("--- Case 2: Adding a Non-Invertible Birefringent Medium (Linear Polarizer) ---")
    # Jones matrix for a horizontal linear polarizer.
    # This matrix is singular (non-invertible) because it projects the light.
    B_non_invertible = np.array([[1, 0], [0, 0]])
    print("Jones Matrix for Horizontal Polarizer (B_non_invertible):")
    print(B_non_invertible)

    # The new system matrix is the product of the elements.
    S_new_non_invertible = B_non_invertible @ T
    print("\nNew System Matrix (S_new = B_non_invertible * T):")
    print(S_new_non_invertible)

    # Check if this new system is invertible.
    det_S2 = np.linalg.det(S_new_non_invertible)
    print(f"Determinant of S_new: {det_S2:.4f}")

    if abs(det_S2) > 1e-9:
        print("Result: The determinant is non-zero. The system IS INVERTIBLE.\n")
    else:
        print("Result: The determinant is zero. The system IS NOT INVERTIBLE.")
        print("The theory fails. Information is lost, and we cannot recover the unique input.\n")

    # Demonstrate the failure by attempting to calculate the inverse.
    try:
        S_inv = np.linalg.inv(S_new_non_invertible)
        print("Inverse calculation succeeded (this should not happen).")
    except np.linalg.LinAlgError as e:
        print(f"As expected, attempting to invert the matrix failed with an error: {e}")

if __name__ == '__main__':
    demonstrate_optical_system_invertibility()