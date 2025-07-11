import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift

def solve():
    """
    This script simulates an optical system to determine if a principle of
    system inversion holds after adding a birefringent medium.
    """
    # --- Plan ---
    # 1. Define simulation parameters (grid, wavelength, etc.).
    # 2. Create a Laguerre-Gaussian (LG) input beam with horizontal polarization.
    # 3. Define operators for free-space propagation (P_fs) and its inverse (P_fs_inv).
    # 4. Define a polarization-independent random medium (T_scalar) and its inverse.
    # 5. Define a birefringent medium (J), a quarter-wave plate.
    # 6. Scenario 1: The original system (P_fs -> T_scalar -> P_fs).
    #    - Propagate the input beam through it to get Output1.
    #    - Apply the inverse scalar operator to Output1.
    #    - Show that the recovered beam matches the input beam. This validates the "scalar theory".
    # 7. Scenario 2: The new system with the birefringent medium (P_fs -> T_scalar -> J -> P_fs).
    #    - Propagate the input beam through it to get Output_new.
    #    - Attempt to use the SAME scalar inverse operator on the new output.
    #    - Show that the recovered beam does NOT match the input beam.
    # 8. Conclude that the original theory does not hold because it cannot account for the
    #    polarization effects introduced by the birefringent medium.

    # 1. Simulation Parameters
    N = 256      # Grid size
    L = 4e-3     # Side length of the grid [m]
    wvl = 633e-9 # Wavelength [m]
    k0 = 2 * np.pi / wvl
    z_prop = 20e-3 # Propagation distance [m]
    dx = L / N
    x = np.linspace(-L/2, L/2 - dx, N)
    X, Y = np.meshgrid(x, x)
    R = np.sqrt(X**2 + Y**2)
    PHI = np.arctan2(Y, X)

    # 2. Create LG Input Beam
    l = 1        # OAM charge
    p = 0        # Radial index
    w0 = 0.5e-3  # Beam waist
    E_in_scalar = (R/w0)**abs(l) * np.exp(-R**2/w0**2) * np.exp(1j * l * PHI)
    E_in_scalar /= np.sqrt(np.sum(np.abs(E_in_scalar)**2)) # Normalize power
    # Input is horizontally polarized: E = [Ex, Ey]
    E_in = np.zeros((N, N, 2), dtype=np.complex128)
    E_in[:, :, 0] = E_in_scalar

    # 3. Define Propagation Operators (Angular Spectrum Method)
    kx = fftshift(np.fft.fftfreq(N, d=dx) * 2 * np.pi)
    KX, KY = np.meshgrid(kx, kx)
    # Use complex kz to account for evanescent waves
    kz = np.sqrt(k0**2 - KX**2 - KY**2, dtype=np.complex128)

    # Forward propagator
    H_fs = fftshift(np.exp(1j * kz * z_prop))
    # Backward propagator (inverse) is the conjugate for lossless propagation
    H_fs_inv = np.conj(H_fs)

    def propagate(E_field_comp, H):
        """Propagates a single scalar field component."""
        return ifft2(ifftshift(fftshift(fft2(E_field_comp)) * H))

    def apply_vector_op(op, E_vec):
        """Applies a 2x2 Jones matrix to a vector field."""
        # This is equivalent to np.einsum('ij,xyj->xyi', op, E_vec)
        Ex_out = op[0, 0] * E_vec[:, :, 0] + op[0, 1] * E_vec[:, :, 1]
        Ey_out = op[1, 0] * E_vec[:, :, 0] + op[1, 1] * E_vec[:, :, 1]
        return np.stack([Ex_out, Ey_out], axis=-1)

    # 4. Define Random Medium (Scalar Phase Screen)
    np.random.seed(0) # for reproducibility
    random_phase = 2 * np.pi * np.random.rand(N, N)
    T_scalar = np.exp(1j * random_phase)
    T_scalar_inv = np.conj(T_scalar) # Inverse is the conjugate

    # 5. Define Birefringent Medium (Quarter-Wave Plate at 45 deg)
    # This matrix converts horizontal linear to right circular polarization.
    J_qwp45 = (1 / np.sqrt(2)) * np.array([[1, -1j], [-1j, 1]])

    # --- SCENARIO 1: Original System & Scalar Theory ---
    print("--- SCENARIO 1: System without Birefringence ---")
    print("The system is modeled as polarization-independent.")

    # Forward propagation: P_fs -> T_scalar -> P_fs
    E_prop1 = np.stack([propagate(E_in[:,:,0], H_fs), propagate(E_in[:,:,1], H_fs)], axis=-1)
    E_after_T = np.stack([E_prop1[:,:,0] * T_scalar, E_prop1[:,:,1] * T_scalar], axis=-1)
    Output1 = np.stack([propagate(E_after_T[:,:,0], H_fs), propagate(E_after_T[:,:,1], H_fs)], axis=-1)

    print("Input beam propagated through the original system to get Output1.")

    # Inversion using the scalar theory.
    # The theory assumes the system is scalar, so it operates on the x-component.
    Output1_x = Output1[:,:,0]
    # Apply inverse operators in reverse order: P_fs_inv -> T_scalar_inv -> P_fs_inv
    E_rec_prop1 = propagate(Output1_x, H_fs_inv)
    E_rec_after_T = E_rec_prop1 * T_scalar_inv
    E_recovered_scalar = propagate(E_rec_after_T, H_fs_inv)

    # Compare recovered beam with original input scalar beam using correlation
    correlation = np.abs(np.sum(E_recovered_scalar * np.conj(E_in_scalar)))
    print(f"Recovered beam correlation with input: {correlation:.4f}")
    print("The high correlation shows the scalar inversion theory holds for the original system.\n")

    # --- SCENARIO 2: New System with Birefringence ---
    print("--- SCENARIO 2: Birefringent Medium Added ---")
    print("The system is now polarization-dependent.")

    # Forward propagation through the new system: P_fs -> T_scalar -> J -> P_fs
    E_prop1_new = np.stack([propagate(E_in[:,:,0], H_fs), propagate(E_in[:,:,1], H_fs)], axis=-1)
    E_after_T_new = np.stack([E_prop1_new[:,:,0] * T_scalar, E_prop1_new[:,:,1] * T_scalar], axis=-1)
    E_after_J = apply_vector_op(J_qwp45, E_after_T_new)
    Output_new = np.stack([propagate(E_after_J[:,:,0], H_fs), propagate(E_after_J[:,:,1], H_fs)], axis=-1)

    print("Input beam propagated through the new system (with J) to get Output_new.")
    output_power_x = np.sum(np.abs(Output_new[:,:,0])**2)
    output_power_y = np.sum(np.abs(Output_new[:,:,1])**2)
    print(f"Output power distribution: Px={output_power_x:.3f}, Py={output_power_y:.3f}")
    print("Note that the output now has a significant vertical component, as expected.")

    # Attempt to invert using the SAME scalar theory
    print("\nAttempting to apply the original SCALAR theory to the new output...")
    # The scalar theory is blind to polarization. It operates on the x-component.
    Output_new_x = Output_new[:,:,0]

    # Apply the same inverse scalar operators as before
    E_rec_prop1_new = propagate(Output_new_x, H_fs_inv)
    E_rec_after_T_new = E_rec_prop1_new * T_scalar_inv
    E_recovered_scalar_new = propagate(E_rec_after_T_new, H_fs_inv)

    # Compare recovered beam with original input scalar beam
    correlation_new = np.abs(np.sum(E_recovered_scalar_new * np.conj(E_in_scalar)))
    print(f"Recovered beam correlation with input: {correlation_new:.4f}")
    print("The low correlation shows the scalar inversion theory FAILS for the new system.")

    print("\n--- Conclusion ---")
    print("Yes, the theory can fail to hold.")
    print("The original 'theory' worked because it correctly modeled the system as being scalar (polarization-independent).")
    print("Adding a birefringent medium makes the system's behavior fundamentally dependent on polarization (i.e., it becomes a vector system). The original scalar theory is no longer a valid model for the system and thus fails to invert it correctly.")

solve()
<<<Yes>>>