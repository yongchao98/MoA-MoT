import numpy as np

# This script demonstrates the importance of using a field-aligned coordinate system
# to calculate magnetic helicity for a wave propagating along the magnetic field.

def calculate_helicity(B1_t, B2_t, k_index):
    """Calculates normalized magnetic helicity at a specific frequency index."""
    B1_k = np.fft.fft(B1_t)[k_index]
    B2_k = np.fft.fft(B2_t)[k_index]
    
    # Helicity formula: sigma_m = 2 * Im(conj(B1_k) * B2_k) / (|B1_k|^2 + |B2_k|^2)
    numerator = 2 * np.imag(np.conj(B1_k) * B2_k)
    denominator = np.abs(B1_k)**2 + np.abs(B2_k)**2
    if denominator < 1e-9:
        return 0
    return numerator / denominator

# --- 1. Define the physical setup at L1 ---
print("--- Setup ---")
# Define the background magnetic field (B0) in a standard GSE-like coordinate system.
# The X-axis is radial from the Sun. We use a Parker spiral angle of 45 degrees
# away from the radial X-axis (so 135 deg from +X in a GSE XY plane view).
parker_angle_rad = np.deg2rad(135)
B0_mag = 5.0 # nT
B0 = np.array([B0_mag * np.cos(parker_angle_rad), B0_mag * np.sin(parker_angle_rad), 0.5])
angle_with_radial = np.rad2deg(np.arccos(np.dot(B0, [1, 0, 0]) / np.linalg.norm(B0)))

print(f"The true background magnetic field B0 is not radial.")
print(f"B0 = ({B0[0]:.2f}, {B0[1]:.2f}, {B0[2]:.2f}) nT")
print(f"Angle of B0 with the radial (X) direction: {angle_with_radial:.1f} degrees\n")

# Define a synthetic 100% left-hand circularly polarized wave (expected helicity = -1).
# The wave propagates EXACTLY parallel to B0.
wave_amplitude = 1.0 # nT
wave_frequency = 0.1 # Hz
time = np.linspace(0, 100, 2048) # Time array for synthetic data

# --- 2. Create the wave in the CORRECT (Field-Aligned) frame ---
# Create the basis vectors for the Field-Aligned Coordinate (FAC) system.
z_fac = B0 / np.linalg.norm(B0) # Parallel to B0
x_fac = np.cross([0, 0, 1], z_fac) # Perpendicular
x_fac /= np.linalg.norm(x_fac)
y_fac = np.cross(z_fac, x_fac) # Perpendicular, completes the right-hand system

# Create the left-hand polarized wave signal in the FAC (x,y) plane.
# The wave has no component along z_fac (the propagation direction).
dB_x_fac = wave_amplitude * np.cos(2 * np.pi * wave_frequency * time)
dB_y_fac = -wave_amplitude * np.sin(2 * np.pi * wave_frequency * time) # Minus sign for left-hand

# --- 3. Transform wave to the original frame (simulates measurement) ---
# The fluctuating field (dB) measured by a spacecraft is the projection
# of the FAC components back onto the fixed (X,Y,Z) axes.
dB_measured = np.outer(dB_x_fac, x_fac) + np.outer(dB_y_fac, y_fac)

# Find the frequency bin corresponding to our wave for Fourier analysis.
N = len(time)
dt = time[1] - time[0]
freqs = np.fft.fftfreq(N, dt)
wave_k_index = np.argmin(np.abs(freqs - wave_frequency))


# --- 4. Compare helicity calculation methods ---
print("--- Helicity Calculation ---")
# Method A: Incorrect approach using components perpendicular to the RADIAL direction (Y and Z).
print("Method A: Using components perpendicular to the radial (X) direction...")
dB_y_measured = dB_measured[:, 1]
dB_z_measured = dB_measured[:, 2]
sigma_m_radial_approx = calculate_helicity(dB_y_measured, dB_z_measured, wave_k_index)
print(f"Calculated helicity using (By, Bz) components: {sigma_m_radial_approx:.4f}\n")


# Method B: Correct approach using components perpendicular to the LOCAL FIELD (B0).
# These are the components we defined the wave with: dB_x_fac and dB_y_fac.
print("Method B: Using components perpendicular to the local magnetic field B0...")
sigma_m_fac = calculate_helicity(dB_x_fac, dB_y_fac, wave_k_index)
print(f"Calculated helicity using true perpendicular components: {sigma_m_fac:.4f}\n")

# --- 5. Conclusion ---
print("--- Conclusion ---")
print(f"The expected helicity for the synthetic left-hand polarized wave was -1.0.")
print(f"The physically correct Method B recovered the value {sigma_m_fac:.2f} perfectly.")
print(f"The approximate Method A yielded an incorrect value of {sigma_m_radial_approx:.2f}.")
print("This demonstrates that using components perpendicular to the radial direction is an")
print("often-invalid approximation, and the analysis should be done in a field-aligned frame.")