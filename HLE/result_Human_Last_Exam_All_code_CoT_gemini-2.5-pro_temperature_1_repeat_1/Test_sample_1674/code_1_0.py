import numpy as np

def run_simulation():
    """
    Simulates light propagation through a complex optical system to test
    if a simplified "scalar" theory holds vs. a full "vectorial" theory.
    """
    # 1. Define System Parameters
    # Let's model the system with N=10 spatial modes and 2 polarizations (H and V).
    # The total number of degrees of freedom is 2 * N.
    N = 10
    total_modes = 2 * N
    np.random.seed(0) # for reproducibility

    # 2. Create the Optical Elements as Matrices
    # --- Random Medium (T) ---
    # A real random medium is depolarizing. We model this with a full 2N x 2N
    # random unitary matrix. A unitary matrix represents a lossless scatterer.
    # We create it from a random complex matrix via QR decomposition.
    random_matrix = np.random.randn(total_modes, total_modes) + 1j * np.random.randn(total_modes, total_modes)
    T, _ = np.linalg.qr(random_matrix)

    # --- Birefringent Plate (B) ---
    # This acts on polarization but not on the spatial mode.
    # Let's model a quarter-wave plate at 45 degrees.
    # Its Jones matrix J mixes H and V polarizations.
    J = 0.5 * np.array([[1, -1j], 
                        [-1j, 1]])
    # The full matrix B is the Kronecker product of J and an identity matrix for the spatial modes.
    B = np.kron(J, np.identity(N))

    # --- Total System (M_true) ---
    # The light first passes through the random medium, then the plate.
    # Note: The order doesn't change the core conclusion.
    M_true = B @ T
    
    print("--- System Setup ---")
    print(f"Simulating with {N} spatial modes and 2 polarizations.")
    print(f"The 'true' system transmission matrix M_true has shape: {M_true.shape}\n")


    # 3. Define the Goal: Focus light to a single spot
    # Our target is to have all power in the first horizontal mode, and zero everywhere else.
    # This is our desired output vector y_target.
    y_target = np.zeros(total_modes, dtype=complex)
    y_target[0] = 1.0 # All power in mode 0 (H-pol, 1st spatial mode)
    target_intensity = np.abs(y_target[0])**2

    # 4. The "Scalar Theory" Approach (Flawed)
    # This approach assumes the random medium T is not depolarizing.
    # The experimenter measures only the H-to-H block of T (T_HH) and assumes
    # the other blocks are either zero or simple. This is an incomplete model.
    T_HH_measured = T[0:N, 0:N]
    
    # Build the flawed model of the scatterer
    T_assumed = np.zeros_like(T)
    T_assumed[0:N, 0:N] = T_HH_measured # Assumes H->H part is known
    T_assumed[N:total_modes, N:total_modes] = T_HH_measured # Assumes V->V is the same
    # Crucially, assumes H->V and V->H are zero (no depolarization)

    # The flawed model of the total system is then:
    M_assumed = B @ T_assumed
    
    # Calculate the required input based on this flawed model
    x_in_scalar = np.linalg.inv(M_assumed) @ y_target
    
    # Now, pass this calculated input through the REAL system
    y_out_scalar = M_true @ x_in_scalar
    
    # Check the result
    achieved_intensity_scalar = np.abs(y_out_scalar[0])**2
    total_intensity_scalar = np.sum(np.abs(y_out_scalar)**2)

    print("--- 1. Simplified 'Scalar' Theory Result ---")
    print("This theory ignores the depolarization of the random medium.")
    print(f"Goal: Focus light to achieve an intensity of {target_intensity:.4f} in one specific mode.")
    print("Calculation: input_vector = inverse(M_assumed) * target_vector")
    print("Result: Propagate this input through the real system M_true.")
    print(f"Achieved intensity in target mode: {achieved_intensity_scalar:.4f}")
    print(f"Total output intensity: {total_intensity_scalar:.4f}")
    print("Conclusion: The scalar theory fails significantly.\n")


    # 5. The "Vectorial Theory" Approach (Correct)
    # This approach uses the full, correct transmission matrix M_true,
    # which would be acquired by a full polarization-sensitive measurement.
    x_in_vectorial = np.linalg.inv(M_true) @ y_target
    
    # Pass this correct input through the real system
    y_out_vectorial = M_true @ x_in_vectorial
    
    # Check the result
    achieved_intensity_vectorial = np.abs(y_out_vectorial[0])**2
    total_intensity_vectorial = np.sum(np.abs(y_out_vectorial)**2)

    print("--- 2. Full 'Vectorial' Theory Result ---")
    print("This theory uses the complete, correct model of the system.")
    print(f"Goal: Focus light to achieve an intensity of {target_intensity:.4f} in one specific mode.")
    print("Calculation: input_vector = inverse(M_true) * target_vector")
    print("Result: Propagate this input through the real system M_true.")
    # The output should be almost exactly the target vector due to numerical precision.
    print(f"Achieved intensity in target mode: {achieved_intensity_vectorial:.4f}")
    print(f"Total output intensity: {total_intensity_vectorial:.4f}")
    print("Conclusion: The full vectorial theory succeeds perfectly.")
    
if __name__ == '__main__':
    run_simulation()