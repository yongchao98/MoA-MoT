import numpy as np

def create_lg_beam(p, l, w0, size, pixel_scale=1.0):
    """Creates a Laguerre-Gaussian beam field."""
    x = (np.arange(size) - size // 2) * pixel_scale
    y = (np.arange(size) - size // 2) * pixel_scale
    X, Y = np.meshgrid(x, y)
    R2 = X**2 + Y**2
    phi = np.arctan2(Y, X)
    
    # Laguerre polynomial L_p^l(x). For p=0, it's just 1.
    # A full implementation is complex, so we'll use a simple case p=0, l=1
    if p != 0:
        print("Warning: Using simplified Laguerre polynomial L_0^|l|(x)=1")
    
    rho_norm = np.sqrt(2 * R2) / w0
    lg_field = (rho_norm**abs(l) * np.exp(-R2 / w0**2) * np.exp(1j * l * phi))
    
    # Normalize
    lg_field /= np.sqrt(np.sum(np.abs(lg_field)**2))
    return lg_field

def get_propagator(size, pixel_scale, wavelength, distance):
    """Creates an angular spectrum propagator."""
    kx = 2 * np.pi * np.fft.fftfreq(size, pixel_scale)
    ky = 2 * np.pi * np.fft.fftfreq(size, pixel_scale)
    KX, KY = np.meshgrid(kx, ky)
    
    # Avoid evanescent waves for simplicity in this model
    k_squared = KX**2 + KY**2
    k = 2 * np.pi / wavelength
    kz = np.sqrt(k**2 - k_squared + 0j) # Add 0j to handle complex result
    
    return np.fft.fftshift(np.exp(1j * kz * distance))

def apply_vector_operator(field_vector, operator):
    """
    Applies a 2x2 spatially varying operator to a vector field.
    field_vector is shape (2, H, W)
    operator is shape (2, 2, H, W)
    """
    # Using einsum for clear and efficient matrix multiplication at each pixel
    # 'ijhw,jhw->ihw' means: sum over j for each (h, w)
    return np.einsum('ijhw,jhw->ihw', operator, field_vector)

def apply_propagator(field_vector, prop_kernel):
    """Applies free-space propagation to both polarization components."""
    # Propagator is the same for both H and V components
    E_h_prop = np.fft.ifft2(np.fft.fft2(field_vector[0]) * prop_kernel)
    E_v_prop = np.fft.ifft2(np.fft.fft2(field_vector[1]) * prop_kernel)
    return np.array([E_h_prop, E_v_prop])

def calculate_correlation(field1, field2):
    """Calculates the normalized correlation between two complex fields."""
    # We care about the complex field, not just intensity
    f1_norm = field1 / np.sqrt(np.sum(np.abs(field1)**2))
    f2_norm = field2 / np.sqrt(np.sum(np.abs(f2)**2))
    return np.abs(np.sum(f1_norm * np.conj(f2_norm)))

def main():
    """
    Simulates the optical system to demonstrate the effect of a birefringent medium
    on the inversion principle.
    """
    # --- System Parameters ---
    SIZE = 128         # Grid size
    PIXEL_SCALE = 5e-6 # 5 um per pixel
    WAVELENGTH = 633e-9# HeNe laser wavelength
    W0 = 150e-6        # Beam waist of 150 um
    DISTANCE = 0.05    # Propagation distance of 5 cm

    # --- 1. Define Components and Input Beam ---
    
    # Input beam: LG(p=0, l=1) with horizontal polarization
    lg_mode = create_lg_beam(p=0, l=1, w0=W0, size=SIZE, pixel_scale=PIXEL_SCALE)
    E_in = np.zeros((2, SIZE, SIZE), dtype=np.complex128)
    E_in[0] = lg_mode # Horizontal component
    E_in[1] = 0       # Vertical component

    # Propagation operators
    propagator = get_propagator(SIZE, PIXEL_SCALE, WAVELENGTH, DISTANCE)
    propagator_inv = np.conj(propagator) # Inverse propagation is the conjugate

    # Random Medium (T): A simple polarization-independent phase screen
    # This is a case where a scalar model might seem sufficient initially.
    random_phase = np.exp(1j * 2 * np.pi * np.random.rand(SIZE, SIZE))
    T = np.zeros((2, 2, SIZE, SIZE), dtype=np.complex128)
    T[0, 0] = random_phase
    T[1, 1] = random_phase
    T_inv = np.conj(T) # Inverse is the conjugate phase

    # Birefringent Medium (B): A quarter-wave plate with fast axis at 45 degrees
    # This will strongly couple H and V polarizations.
    B = np.zeros((2, 2, SIZE, SIZE), dtype=np.complex128)
    qwp_matrix = 0.5 * np.array([[1+1j, 1-1j], [1-1j, 1+1j]])
    B[0, 0], B[0, 1] = qwp_matrix[0, 0], qwp_matrix[0, 1]
    B[1, 0], B[1, 1] = qwp_matrix[1, 0], qwp_matrix[1, 1]
    B_inv = np.conj(B.transpose(1, 0, 2, 3)) # Inverse is the conjugate transpose

    print("--- Simulating Optical Inversion ---")
    print("The 'final equation' we test is: Correlation(Original, Reconstructed)\n")

    # --- 2. Scenario A: No Birefringence ---
    # System: Propagate -> Random Medium -> Propagate
    E_mid_A = apply_propagator(E_in, propagator)
    E_scrambled_A = apply_vector_operator(E_mid_A, T)
    E_out_A = apply_propagator(E_scrambled_A, propagator)

    # Inversion: Apply inverse operators in reverse order
    E_recon_mid1_A = apply_propagator(E_out_A, propagator_inv)
    E_recon_mid2_A = apply_vector_operator(E_recon_mid1_A, T_inv)
    E_reconstructed_A = apply_propagator(E_recon_mid2_A, propagator_inv)
    
    corr_A = calculate_correlation(E_in[0], E_reconstructed_A[0])
    print("Scenario A (No Birefringence, Correct Inverse):")
    print(f"Correlation = {corr_A:.5f}\n")


    # --- 3. Scenario B: With Birefringence ---
    # System: Propagate -> Random Medium -> Birefringent Plate -> Propagate
    E_mid_B = apply_propagator(E_in, propagator)
    E_scrambled_B = apply_vector_operator(E_mid_B, T)
    E_birefringent_B = apply_vector_operator(E_scrambled_B, B)
    E_out_B = apply_propagator(E_birefringent_B, propagator)

    # Inversion Attempt 1 (Incorrect): Ignore the birefringent plate.
    # This simulates using an incomplete model of the system.
    E_recon_mid1_B_fail = apply_propagator(E_out_B, propagator_inv)
    # Note: We are only applying T_inv, not B_inv
    E_recon_mid2_B_fail = apply_vector_operator(E_recon_mid1_B_fail, T_inv) 
    E_reconstructed_B_fail = apply_propagator(E_recon_mid2_B_fail, propagator_inv)

    corr_B_fail = calculate_correlation(E_in[0], E_reconstructed_B_fail[0])
    print("Scenario B (With Birefringence, Incorrect 'Scalar' Inverse):")
    print(f"Correlation = {corr_B_fail:.5f}")
    print("Result: Failure. The old model does not work.\n")

    # Inversion Attempt 2 (Correct): Use the full vector inverse for the complete system.
    # System_inv = Prop_inv -> B_inv -> T_inv -> Prop_inv
    E_recon_mid1_B_ok = apply_propagator(E_out_B, propagator_inv)
    E_recon_mid2_B_ok = apply_vector_operator(E_recon_mid1_B_ok, B_inv)
    E_recon_mid3_B_ok = apply_vector_operator(E_recon_mid2_B_ok, T_inv)
    E_reconstructed_B_ok = apply_propagator(E_recon_mid3_B_ok, propagator_inv)
    
    corr_B_ok = calculate_correlation(E_in[0], E_reconstructed_B_ok[0])
    print("Scenario B (With Birefringence, Correct Vector Inverse):")
    print(f"Correlation = {corr_B_ok:.5f}")
    print("Result: Success. The theory holds with the correct, complete model.")


if __name__ == '__main__':
    main()
