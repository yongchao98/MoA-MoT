import math

# --- System Parameters ---
# These parameters are based on the problem description and typical values.

# The active length of the Spin Valve sensor
sensor_length_um = 5.0
sensor_length_m = sensor_length_um * 1e-6

# Electronics bandpass filter settings
# A high-pass filter removes low-frequency drift.
high_pass_cutoff_hz = 500
# A low-pass filter removes high-frequency noise.
low_pass_cutoff_hz = 10000

# --- Scenario 1: Correct Magnet Positioning ---
print("--- Scenario 1: Correct Magnet Positioning ---")
print("In this case, particles flow near the center of the channel at optimal velocity.")

# Particle velocity in the center of the channel (e.g., 10 mm/s)
velocity_correct_m_s = 10 * 1e-3

# The signal frequency is roughly velocity / sensor length
signal_frequency_correct_hz = velocity_correct_m_s / sensor_length_m

print(f"Optimal particle velocity: {velocity_correct_m_s * 1000:.1f} mm/s")
print(f"Sensor length: {sensor_length_um} um")
print(f"Resulting Signal Frequency = ({velocity_correct_m_s} m/s) / ({sensor_length_m} m) = {signal_frequency_correct_hz:.0f} Hz")
print(f"System Bandwidth: {high_pass_cutoff_hz} Hz - {low_pass_cutoff_hz} Hz")

if high_pass_cutoff_hz < signal_frequency_correct_hz < low_pass_cutoff_hz:
    print("Result: Signal is INSIDE the system bandwidth and will be detected.\n")
else:
    print("Result: Signal is OUTSIDE the system bandwidth and may be rejected.\n")


# --- Scenario 2: Improper Magnet Positioning ---
print("--- Scenario 2: Improper Magnet Positioning ---")
print("A lateral magnetic force pushes particles to the channel wall, reducing their velocity.")

# Particle velocity near the channel wall is significantly lower (e.g., 1 mm/s)
velocity_improper_m_s = 1 * 1e-3

# Calculate the new, lower signal frequency
signal_frequency_improper_hz = velocity_improper_m_s / sensor_length_m

print(f"Reduced particle velocity: {velocity_improper_m_s * 1000:.1f} mm/s")
print(f"Sensor length: {sensor_length_um} um")
print(f"Resulting Signal Frequency = ({velocity_improper_m_s} m/s) / ({sensor_length_m} m) = {signal_frequency_improper_hz:.0f} Hz")
print(f"System Bandwidth: {high_pass_cutoff_hz} Hz - {low_pass_cutoff_hz} Hz")

if high_pass_cutoff_hz < signal_frequency_improper_hz < low_pass_cutoff_hz:
    print("Result: Signal is INSIDE the system bandwidth and will be detected.")
else:
    print("Result: Signal is OUTSIDE the system bandwidth and will be rejected by the high-pass filter.")
