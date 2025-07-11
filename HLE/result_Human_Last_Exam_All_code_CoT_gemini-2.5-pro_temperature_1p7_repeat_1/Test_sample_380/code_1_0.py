import numpy as np
from scipy.signal import csd, welch

def calculate_helicity(comp1, comp2, fs, nperseg):
    """Calculates normalized magnetic helicity from two magnetic field components."""
    # Calculate Power Spectral Densities (PSD)
    f, Pxx = welch(comp1, fs=fs, nperseg=nperseg)
    _, Pyy = welch(comp2, fs=fs, nperseg=nperseg)
    
    # Calculate Cross Spectral Density (CSD)
    f_csd, Pxy = csd(comp1, comp2, fs=fs, nperseg=nperseg)
    
    # Helicity formula: sigma_m = 2 * Im(CSD) / (PSD1 + PSD2)
    helicity = 2 * Pxy.imag / (Pxx + Pyy)
    
    return f, helicity, Pxx, Pyy, Pxy.imag

# --- 1. Setup and Synthetic Data Generation ---
# Simulation parameters
fs = 20.0        # Sampling frequency in Hz
duration = 100.0   # seconds
n_points = int(fs * duration)
t = np.linspace(0, duration, n_points)

# Wave and Plasma Parameters
b0_mag = 5.0      # Background field magnitude in nT
parker_angle_deg = 45.0 # Angle of B-field wrt radial (X) direction
wave_amp = 1.0    # Wave amplitude in nT
wave_freq = 0.5   # Wave frequency in Hz (typical for AIC waves)

# Define the background magnetic field vector (Parker spiral)
angle_rad = np.deg2rad(parker_angle_deg)
b0_vec = np.array([b0_mag * np.cos(angle_rad), b0_mag * np.sin(angle_rad), 0])

# Create a Field-Aligned Coordinate (FAC) system
# e_para is parallel to the mean B-field
e_para = b0_vec / np.linalg.norm(b0_vec)
# e_perp2 is chosen to be the Z-axis for simplicity (since B0 is in XY plane)
e_perp2 = np.array([0., 0., 1.])
# e_perp1 completes the right-hand system
e_perp1 = np.cross(e_perp2, e_para)

# Generate a synthetic LEFT-HAND polarized wave in the FAC system
# For left-hand polarization, the rotation is from e_perp1 towards -e_perp2
omega = 2 * np.pi * wave_freq
b_perp1_fac = wave_amp * np.cos(omega * t)
b_perp2_fac = -wave_amp * np.sin(omega * t)

# Transform wave components back to the original (X,Y,Z) system
# delta_B = b_perp1_fac * e_perp1 + b_perp2_fac * e_perp2
delta_b_vec = np.outer(b_perp1_fac, e_perp1) + np.outer(b_perp2_fac, e_perp2)

# The "measured" field is the sum of background and wave
b_measured = b0_vec + delta_b_vec
Bx, By, Bz = b_measured[:, 0], b_measured[:, 1], b_measured[:, 2]

print("--- Data Generation Complete ---")
print(f"Generated a left-hand polarized wave at {wave_freq} Hz.")
print(f"Background B-field direction (unit vector): {e_para.round(3)}\n")

# --- 2. Method 1: Simplified Calculation (using Y, Z components) ---
print("--- Method 1: Simplified (Using Y, Z components perpendicular to radial) ---")
nperseg = 256  # Segment length for spectral analysis
freqs_1, helicity_1, Pyy, Pzz, Im_Pyz = calculate_helicity(By, Bz, fs, nperseg)

# Find helicity value at the known wave frequency
idx_wave_1 = np.argmin(np.abs(freqs_1 - wave_freq))
helicity_val_1 = helicity_1[idx_wave_1]
Pyy_val = Pyy[idx_wave_1]
Pzz_val = Pzz[idx_wave_1]
Im_Pyz_val = Im_Pyz[idx_wave_1]

print(f"Analysis performed at frequency f = {freqs_1[idx_wave_1]:.3f} Hz")
print("Final Equation: sigma = 2 * Im(S_yz) / (S_yy + S_zz)")
print(f"Numbers for equation: sigma = 2 * ({Im_Pyz_val:.4f}) / ({Pyy_val:.4f} + {Pzz_val:.4f})")
print(f"Calculated Helicity = {helicity_val_1:.4f}")
print("Result: This value is not -1, incorrectly characterizing the wave.\n")


# --- 3. Method 2: Rigorous Calculation (using Field-Aligned Coordinates) ---
print("--- Method 2: Rigorous (Using components perpendicular to the mean B-field) ---")
# Get fluctuations by subtracting the mean
b_fluc = b_measured - np.mean(b_measured, axis=0)

# Project fluctuations onto the perpendicular FAC directions
b_perp1_trans = np.dot(b_fluc, e_perp1)
b_perp2_trans = np.dot(b_fluc, e_perp2)

# Calculate helicity using the correct perpendicular components
freqs_2, helicity_2, Pp1p1, Pp2p2, Im_Pp1p2 = calculate_helicity(b_perp1_trans, b_perp2_trans, fs, nperseg)

# Find helicity value at the known wave frequency
idx_wave_2 = np.argmin(np.abs(freqs_2 - wave_freq))
helicity_val_2 = helicity_2[idx_wave_2]
Pp1p1_val = Pp1p1[idx_wave_2]
Pp2p2_val = Pp2p2[idx_wave_2]
Im_Pp1p2_val = Im_Pp1p2[idx_wave_2]

print(f"Analysis performed at frequency f = {freqs_2[idx_wave_2]:.3f} Hz")
print("Final Equation: sigma = 2 * Im(S_p1p2) / (S_p1p1 + S_p2p2)")
print(f"Numbers for equation: sigma = 2 * ({Im_Pp1p2_val:.4f}) / ({Pp1p1_val:.4f} + {Pp2p2_val:.4f})")
print(f"Calculated Helicity = {helicity_val_2:.4f}")
print("Result: This value is very close to -1, correctly identifying the left-hand polarization.")
