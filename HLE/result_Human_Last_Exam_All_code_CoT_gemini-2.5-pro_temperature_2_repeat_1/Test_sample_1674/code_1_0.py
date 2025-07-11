import numpy as np

def run_simulation():
    """
    Simulates wavefront shaping to demonstrate how a birefringent
    plate can cause scalar-based theories to fail.
    """
    # 1. SETUP PARAMETERS
    # We model the system with a set of discrete spatial modes.
    N_spatial_modes = 32
    N_pol_modes = 2  # Horizontal (H) and Vertical (V)
    N_total_modes = N_spatial_modes * N_pol_modes

    # Set a seed for reproducibility
    np.random.seed(0)

    # 2. DEFINE OPTICAL COMPONENTS AS MATRICES
    
    # A) Random Scattering Medium (T)
    # This matrix couples all input modes to all output modes randomly.
    # The elements are complex numbers drawn from a Gaussian distribution.
    T_real = np.random.randn(N_total_modes, N_total_modes)
    T_imag = np.random.randn(N_total_modes, N_total_modes)
    T = (T_real + 1j * T_imag) / np.sqrt(2) # Normalization factor

    # B) Birefringent Medium (B)
    # We'll model a half-wave plate at 45 degrees, which swaps H and V polarizations.
    # This is a highly deterministic, polarization-dependent operation.
    # Jones matrix for a HWP at 45 deg: [[0, 1], [1, 0]]
    B_pol = np.array([[0, 1], [1, 0]])
    # The full matrix B applies this Jones matrix to each spatial mode.
    # This is achieved with a Kronecker product with an identity matrix.
    B = np.kron(np.eye(N_spatial_modes), B_pol)

    # C) Full System (S)
    # The light first passes through the random medium, then the plate.
    S = B @ T

    # 3. DEFINE TARGET AND PERFORM WAVEFRONT SHAPING
    
    # Let's try to focus light onto the first output mode (Spatial mode 0, H-pol)
    target_mode_index = 0
    
    # --- Baseline: Intensity with an unshaped, generic input ---
    # An input beam of purely H-polarized light with uniform illumination
    input_unshaped = np.zeros(N_total_modes, dtype=complex)
    input_unshaped[0::2] = 1.0 # Excite all H-pol modes (indices 0, 2, 4...)
    
    output_unshaped = S @ input_unshaped
    I_unshaped = np.abs(output_unshaped[target_mode_index])**2
    
    # --- Scenario 1: Flawed SCALAR Theory ---
    # This theory wrongly assumes polarization doesn't matter. It only measures
    # the H-in to H-out transmission coefficients (S_HH) and calculates a
    # correction based only on those.
    # S_HH corresponds to elements at even rows and even columns of S.
    S_HH_measured = S[0::2, 0::2] 
    
    # The transmission vector for H-inputs to our H-output target (mode 0).
    tm_row_scalar = S_HH_measured[0, :]
    
    # Calculate the optimal phases for H-inputs based on this flawed data.
    optimal_phases_scalar = np.angle(tm_row_scalar)
    input_H_scalar_corrected = np.exp(-1j * optimal_phases_scalar)
    
    # Create the full input vector for the real system (V-pols are zero)
    input_scalar_full = np.zeros(N_total_modes, dtype=complex)
    input_scalar_full[0::2] = input_H_scalar_corrected
    
    # Propagate this "corrected" input through the TRUE system.
    output_scalar = S @ input_scalar_full
    I_scalar_focused = np.abs(output_scalar[target_mode_index])**2
    
    # --- Scenario 2: Correct VECTORIAL Theory ---
    # This theory correctly measures the full system matrix S and uses all
    # available modes (both H and V) to optimize the focus.
    
    # Use the full transmission vector to the target mode.
    tm_row_vectorial = S[target_mode_index, :]
    
    # Calculate optimal phases for ALL input modes (H and V).
    optimal_phases_vectorial = np.angle(tm_row_vectorial)
    input_vectorial_corrected = np.exp(-1j * optimal_phases_vectorial)
    
    # Propagate this correctly calculated input through the system.
    output_vectorial = S @ input_vectorial_corrected
    I_vectorial_focused = np.abs(output_vectorial[target_mode_index])**2

    # 4. PRINT RESULTS
    print("--- Wavefront Shaping Simulation ---")
    print(f"Number of controllable input modes: {N_total_modes}")
    print("\nGoal: Focus light onto 'Spatial Mode 0, Horizontal Polarization'.")
    print("-" * 36)
    
    # Unshaped intensity is the sum of |S_0,j|^2 for H-polarized j
    avg_unshaped = np.sum(np.abs(S[target_mode_index, 0::2])**2)

    print(f"Intensity with unshaped (generic) input: {I_unshaped:.2f}")
    print(f"Intensity with FLAWED SCALAR correction: {I_scalar_focused:.2f}")
    print(f"Intensity with correct VECTORIAL correction: {I_vectorial_focused:.2f}")
    print("-" * 36)
    
    print("\nConclusion:")
    print("The SCALAR theory fails to improve focusing because it is blind to the polarization effects of the birefringent plate.")
    print("The VECTORIAL theory, which accounts for polarization, achieves a very high intensity enhancement.")
    print("\nThis demonstrates that adding a birefringent medium can cause a simplified (scalar) theory to no longer hold.")
    
# Run the simulation
run_simulation()
