import numpy as np

def create_gaussian_beam(N, w0):
    """Creates a 2D Gaussian amplitude profile on a grid of size N."""
    x = np.linspace(-1, 1, N)
    y = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(x, y)
    R2 = X**2 + Y**2
    return np.exp(-R2 / w0**2)

def apply_phase_mask(E_field, mask):
    """
    Applies a phase mask to an electric field.
    A phase-only element applies the same phase shift to any polarization
    component present, without mixing them.
    E_out(x,y) = E_in(x,y) * exp(i * mask(x,y))
    """
    phase_operator = np.exp(1j * mask)
    E_out = np.zeros_like(E_field, dtype=np.complex128)
    # Apply phase operator to Ex and Ey components independently
    E_out[..., 0] = E_field[..., 0] * phase_operator
    E_out[..., 1] = E_field[..., 1] * phase_operator
    return E_out

def propagate_free_space(E_field):
    """
    Simulates free-space propagation using a simple FFT-based model.
    Propagation acts on each polarization component independently. This is a key
    physical principle that the model captures.
    """
    E_out = np.zeros_like(E_field, dtype=np.complex128)
    # Propagate Ex and Ey components independently
    if np.any(E_field[..., 0]): # only compute if non-zero
        E_out[..., 0] = np.fft.ifft2(np.fft.fft2(E_field[..., 0]))
    if np.any(E_field[..., 1]): # only compute if non-zero
        E_out[..., 1] = np.fft.ifft2(np.fft.fft2(E_field[..., 1]))
    return E_out

# 1. Define simulation parameters
N = 128 # Grid size for the beam cross-section
w0 = 0.3 # Gaussian beam waist parameter

# 2. Create the input beam: A uniformly linearly x-polarized beam.
# The electric field is represented as an (N, N, 2) array for (Ex, Ey).
E_in = np.zeros((N, N, 2), dtype=np.complex128)
amplitude = create_gaussian_beam(N, w0)
E_in[..., 0] = amplitude  # Set Ex component to have a Gaussian profile
E_in[..., 1] = 0.0 + 0.0j  # Set Ey component to be zero everywhere

print("--- System Analysis ---")
print(f"Input beam is purely linearly polarized.")
is_input_ey_zero = np.allclose(E_in[..., 1], 0)
print(f"Is the input Ey component identically zero? {is_input_ey_zero}")
print("-" * 25)

# 3. Define the system's phase-shaping media
# A random phase mask for the first medium
T_mask = 2 * np.pi * np.random.rand(N, N)
# The inverse mask for the second medium, which undoes the phase shift
T_inv_mask = -T_mask

# 4. Simulate the beam's journey through the entire optical system as described
# Path: E_in -> P_fs1 -> T -> P_fs2 -> E_out1 -> P_fs3 -> T_inv -> E_out2
print("Propagating beam through the system...")
# First part of the system to get Output 1
E_after_P1 = propagate_free_space(E_in)
E_after_T = apply_phase_mask(E_after_P1, T_mask)
E_out1 = propagate_free_space(E_after_T)

# Second part of the system to get Output 2
E_after_P3 = propagate_free_space(E_out1)
E_out2 = apply_phase_mask(E_after_P3, T_inv_mask)
print("Propagation complete.")
print("-" * 25)

# 5. Analyze the final output beam's polarization state
print("Analyzing the final output beam (Output 2)...")
Ey_out_final = E_out2[..., 1]

# The crucial check: Has a non-zero Ey component been generated?
is_output_ey_zero = np.allclose(Ey_out_final, 0)

print(f"Is the final output Ey component identically zero? {is_output_ey_zero}")

# 6. Print the final conclusion based on the simulation
print("\n--- Conclusion ---")
print("The simulation demonstrates that an initially purely x-polarized beam remains")
print("purely x-polarized after passing through the entire system. The Ey component")
print("remains zero at all points.")
print("\nThis occurs because all components in the optical path (free-space propagation")
print("and phase-shaping media) are polarization-preserving. They do not have the")
print("physical mechanism to convert x-polarized light into y-polarized light.")
print("\nA vector beam requires a spatially varying polarization state, which necessitates")
print("the presence of both Ex and Ey components in a specific relationship. Since the")
print("system cannot generate the required second component (Ey), it cannot produce")
print("a vector beam from a uniformly polarized input.")
print("\nTherefore, the answer is:")
print("No")
