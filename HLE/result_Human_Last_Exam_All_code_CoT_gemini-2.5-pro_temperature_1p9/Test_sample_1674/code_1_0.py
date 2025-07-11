import numpy as np
from scipy.stats import unitary_group

def run_simulation():
    """
    Simulates an optical system to test if a correction theory holds
    after adding a birefringent medium.
    """
    # --- 1. Define System Parameters ---
    
    # Let's model the spatial part of the beam with a 4-dimensional vector space.
    # This represents, for example, a basis of 4 spatial modes.
    N_spatial = 4
    
    # Polarization is a 2-dimensional space (Horizontal, Vertical).
    N_pol = 2
    
    print("--- System Setup ---")
    print(f"Modeling spatial state with {N_spatial} modes.")
    print(f"Modeling polarization state with {N_pol} components (H, V).\n")

    # --- 2. Define Initial State and Operators ---
    
    # Input Beam: A random spatial mode with pure Horizontal polarization.
    spatial_in = np.random.rand(N_spatial) + 1j * np.random.rand(N_spatial)
    spatial_in /= np.linalg.norm(spatial_in) # Normalize
    pol_in = np.array([1, 0]) # Horizontal polarization [H, V]
    
    # The full state is the tensor product of polarization and spatial parts.
    # We use np.kron to create this combined state vector.
    psi_in = np.kron(pol_in, spatial_in)

    # Random Medium (T): A random unitary matrix acting only on the spatial part.
    # Unitary means it scrambles the phase/shape but conserves energy.
    T_spatial = unitary_group.rvs(N_spatial, random_state=42)
    T_inv_spatial = T_spatial.T.conj() # Inverse of a unitary matrix is its conjugate transpose
    
    # Full operator for T: It acts on spatial part, identity on polarization.
    T = np.kron(np.identity(N_pol), T_spatial)
    T_inv = np.kron(np.identity(N_pol), T_inv_spatial)

    # Birefringent Medium (B): A quarter-wave plate at 45 degrees.
    # This operator acts only on the polarization part. It will change
    # linear horizontal polarization into right-circular polarization.
    B_pol = 0.5 * np.array([[1, -1j], [-1j, 1]])
    
    # Full operator for B: It acts on polarization part, identity on spatial.
    B = np.kron(B_pol, np.identity(N_spatial))

    # Free-Space Propagation (P): For simplicity, we model this as identity.
    # This represents a system where elements are placed directly adjacent,
    # making the core logic easier to see. The user's theory holds in this case.
    P = np.identity(N_spatial * N_pol)
    
    
    # --- 3. Scenario 1: Original System (No Birefringence) ---
    
    print("--- Scenario 1: Original System (without Birefringence) ---")
    
    # The full system operator S = P * T * P
    S_original = P @ T @ P
    
    # Propagate the input beam to get output 1
    psi_out1 = S_original @ psi_in
    
    # The theory proposes a way to find a new input (output 2) that
    # should regenerate output 1.
    # Output 2 = T_inverse * P * Output 1
    psi_out2 = T_inv @ P @ psi_out1
    
    # Test the theory: Does propagating output 2 through the system
    # give back output 1?
    # Test_Output = S_original * Output 2
    psi_test_1 = S_original @ psi_out2
    
    # Check the error between the result and the expected output
    error_1 = np.linalg.norm(psi_test_1 - psi_out1)
    
    print(f"Propagating the derived beam 'Output 2' back through the original system.")
    print("Equation to check: S_original * Output_2 == Output_1")
    print(f"Error (norm of difference): {error_1:.2e}")
    if error_1 < 1e-9:
        print("Result: The theory holds for the original system (as expected in this simplified model).\n")
    else:
        print("Result: The theory fails for the original system.\n")


    # --- 4. Scenario 2: New System (With Birefringence) ---

    print("--- Scenario 2: System with Birefringent Plate Added ---")
    
    # The new system operator includes B: S_new = P * B * T * P
    S_new = P @ B @ T @ P
    
    # Propagate the same input beam to get the new output
    psi_out1_new = S_new @ psi_in
    
    # IMPORTANT: The procedure to generate the correction beam (output 2)
    # remains the same. It only knows about T_inverse, not about B.
    # Output 2_new = T_inverse * P * Output 1_new
    psi_out2_new = T_inv @ P @ psi_out1_new
    
    # Test the theory again with the new system:
    # Does propagating output 2_new through the new system give back output 1_new?
    # Test_Output_new = S_new * Output 2_new
    psi_test_2 = S_new @ psi_out2_new
    
    # Check the error
    error_2 = np.linalg.norm(psi_test_2 - psi_out1_new)

    print(f"Propagating the derived beam 'Output 2_new' back through the new system.")
    print("Equation to check: S_new * Output_2_new == Output_1_new")
    print(f"Error (norm of difference): {error_2:.2f}")
    if error_2 > 1e-9:
        print("Result: The theory FAILS when the birefringent plate is added.\n")
    else:
        print("Result: The theory still holds.\n")
        
    print("Conclusion: The birefringent plate altered the polarization in a way that the T_inverse operator alone could not correct for, causing the theory to fail.")

run_simulation()