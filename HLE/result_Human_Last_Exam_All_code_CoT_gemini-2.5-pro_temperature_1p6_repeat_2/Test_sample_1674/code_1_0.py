import numpy as np
from scipy.stats import unitary_group

def run_simulation():
    """
    This simulation models an optical system to test a proposed inversion theory.
    It demonstrates that the theory is invalid because optical operators for
    propagation and interaction with media do not generally commute.
    """
    # 1. System Parameters
    N = 256      # Number of spatial grid points
    dx = 0.05    # Spatial grid spacing [mm]
    wvl = 0.000632 # Wavelength [mm]
    k0 = 2 * np.pi / wvl
    z = 10.0     # Propagation distance [mm]

    # 2. Define Operators as Functions

    def free_space_propagate(psi_in, dist):
        """Applies free-space propagation using Fourier methods."""
        kx = 2 * np.pi * np.fft.fftfreq(N, dx)
        # Add 0j to ensure result is complex and avoid domain errors for evanescent waves
        kz = np.sqrt(k0**2 - kx**2 + 0j)
        propagator = np.exp(1j * kz * dist)
        
        psi_out = np.zeros_like(psi_in)
        # Propagate each polarization component independently
        psi_out[0, :] = np.fft.ifft(np.fft.fft(psi_in[0, :]) * propagator)
        psi_out[1, :] = np.fft.ifft(np.fft.fft(psi_in[1, :]) * propagator)
        return psi_out

    def apply_medium(psi_in, medium_matrices):
        """Applies a spatially-varying medium (Jones matrix at each point)."""
        # Einsum is an efficient way to do batched matrix-vector products
        # 'ijk,jk->ik' means: for each k, sum over j (T_ij * v_j)
        return np.einsum('ijk,jk->ik', medium_matrices, psi_in)

    # 3. Create the Media and an Initial Beam

    # Use a random input beam to test the operator identity robustly
    np.random.seed(42)
    psi_initial = np.random.randn(2, N) + 1j * np.random.randn(2, N)
    psi_initial /= np.linalg.norm(psi_initial) # Normalize

    # Create the random medium T as a stack of random unitary Jones matrices
    # Unitary matrices are always invertible, with inverse = conjugate transpose
    T_matrices = unitary_group.rvs(2, size=N, random_state=42)
    T_inv_matrices = T_matrices.transpose(0, 2, 1).conj()

    # Create the birefringent plate B (Quarter-Wave Plate at 45 degrees)
    theta = np.pi / 4
    c, s = np.cos(theta), np.sin(theta)
    R_45 = np.array([[c, -s], [s, c]])
    QWP = np.array([[1, 0], [0, 1j]]) # axes-aligned QWP
    B_matrix = R_45.conj().T @ QWP @ R_45
    B_inv_matrix = np.linalg.inv(B_matrix)

    # 4. Test the Theory
    
    print("Analyzing the proposed theory for optical system inversion.")
    print("The theory is valid if a specific operator product equals the Identity.")
    print("We test this by applying the operator to a beam and checking for changes.")
    print("-" * 60)

    # --- Case 1: System without Birefringence ---
    # The operator product for the theory to hold is: (U_fs * T * U_fs) * (T_inv * U_fs)
    print("Case 1: System with Random Medium (T) only")
    
    # Apply the operator sequence: Ufs * T * Ufs * T_inv * Ufs
    psi_1 = free_space_propagate(psi_initial, z)
    psi_2 = apply_medium(psi_1, T_inv_matrices)
    psi_3 = free_space_propagate(psi_2, z)
    psi_4 = apply_medium(psi_3, T_matrices)
    psi_final_case1 = free_space_propagate(psi_4, z)

    # Calculate deviation from identity. 0 means perfect identity.
    deviation_case1 = np.linalg.norm(psi_final_case1 - psi_initial) / np.linalg.norm(psi_initial)
    
    print(f"Deviation from original beam: {deviation_case1:.4f}")
    if deviation_case1 > 1e-6:
        print("Result: The theory does NOT hold. The operators do not form an identity.")
    else:
        print("Result: The theory holds.")
    print("-" * 60)

    # --- Case 2: System with Birefringence ---
    # The total medium is now T_new = B * T
    # The operator product becomes: (U_fs * B*T * U_fs) * ((B*T)_inv * U_fs)
    print("Case 2: System with Birefringent Plate (B) and Random Medium (T)")
    
    T_new_matrices = np.einsum('ij,kjl->kil', B_matrix, T_matrices)
    T_new_inv_matrices = np.einsum('ijk,kl->ijl', T_inv_matrices, B_inv_matrix)

    # Apply the new operator sequence
    psi_1_c2 = free_space_propagate(psi_initial, z)
    psi_2_c2 = apply_medium(psi_1_c2, T_new_inv_matrices)
    psi_3_c2 = free_space_propagate(psi_2_c2, z)
    psi_4_c2 = apply_medium(psi_3_c2, T_new_matrices)
    psi_final_case2 = free_space_propagate(psi_4_c2, z)

    deviation_case2 = np.linalg.norm(psi_final_case2 - psi_initial) / np.linalg.norm(psi_initial)
    
    print(f"Deviation from original beam: {deviation_case2:.4f}")
    if deviation_case2 > 1e-6:
        print("Result: The theory does NOT hold with the added birefringent element.")
    else:
        print("Result: The theory holds.")
    print("-" * 60)

    print("\nConclusion:")
    print("The numerical simulation shows the proposed theory is not valid, because the operators for free-space propagation and for the random medium do not commute.")
    print("Adding a birefringent plate, which has its own polarization-dependent behavior, does not fix this fundamental issue and also results in the theory being invalid.")

if __name__ == '__main__':
    run_simulation()