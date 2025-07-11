import numpy as np

def demonstrate_optical_theory():
    """
    Demonstrates how adding an unaccounted-for birefringent layer
    can cause a time-reversal/inversion theory to fail.
    """
    # Set a random seed for reproducibility
    np.random.seed(42)

    # 1. Define the initial optical state (e.g., a horizontally polarized beam)
    # This is a Jones vector [H, V]
    psi_input = np.array([1.0 + 0.0j, 0.0 + 0.0j])

    # 2. Define the optical elements as 2x2 Jones Matrices

    # The "random medium" T. We create a random unitary matrix to represent
    # a complex but lossless medium (like a polarization scrambler).
    Z = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)
    Q, R = np.linalg.qr(Z)
    T = Q
    # The inverse of a unitary matrix is its conjugate transpose.
    T_inv = T.conj().T

    # The birefringent medium B. It applies a different phase shift to
    # horizontal (phi_H) and vertical (phi_V) components.
    phi_H = np.pi / 4  # 45 degrees phase shift
    phi_V = -np.pi / 2 # -90 degrees phase shift
    B = np.array([
        [np.exp(1j * phi_H), 0],
        [0, np.exp(1j * phi_V)]
    ])

    print("--- Setup ---")
    print(f"Initial Input State (psi_input):\n{psi_input}\n")
    print(f"Random Medium Matrix (T):\n{np.round(T, 3)}\n")
    print(f"Birefringent Medium Matrix (B):\n{np.round(B, 3)}\n")
    print("-" * 20 + "\n")


    # --- Case 1: The Original Theory ---
    # Propagate through the random medium T and then apply its inverse T_inv.
    # We expect to recover the original input.

    print("--- Case 1: Original System (Input -> T -> T_inv -> Output) ---")
    # Forward propagation: Output1 = T * Input
    output_1 = T @ psi_input
    
    # Reverse propagation: Output2 = T_inv * Output1
    recovered_output_1 = T_inv @ output_1

    print("Equation: psi_recovered = T_inv * T * psi_input")
    print(f"T_inv * T gives:\n{np.round(T_inv @ T, 3)}")
    print("\nThis is the identity matrix, so the input should be recovered.")
    print(f"Final Recovered State:\n{np.round(recovered_output_1, 3)}")
    print(f"Is recovered state close to input? {np.allclose(psi_input, recovered_output_1)}\n")
    print("-" * 20 + "\n")


    # --- Case 2: Theory with Birefringent Layer ---
    # Add the birefringent layer B in the forward path, but do NOT
    # account for it in the reverse path.
    
    print("--- Case 2: System with Birefringent Layer (Input -> T -> B -> T_inv -> Output) ---")
    # Forward propagation: New_Output = B * T * Input
    new_output_intermediate = T @ psi_input
    new_output = B @ new_output_intermediate

    # Reverse propagation with the OLD inverse: New_Recovered = T_inv * New_Output
    recovered_output_2 = T_inv @ new_output
    
    # Full operation is: T_inv * B * T
    total_operator = T_inv @ B @ T

    print("Equation: psi_recovered_new = T_inv * B * T * psi_input")
    print(f"The total operator (T_inv * B * T) gives:\n{np.round(total_operator, 3)}")
    print("\nThis is NOT the identity matrix, so the theory fails.")
    print(f"Final Recovered State:\n{np.round(recovered_output_2, 3)}")
    print(f"Is recovered state close to input? {np.allclose(psi_input, recovered_output_2)}\n")


if __name__ == '__main__':
    demonstrate_optical_theory()
