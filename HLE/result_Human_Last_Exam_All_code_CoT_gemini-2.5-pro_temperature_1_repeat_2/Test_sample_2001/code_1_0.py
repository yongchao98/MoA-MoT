import math

# --- System Parameters ---
# Let's define the parameters of our simulated flow cytometry system.

# The effective length of the Spin Valve sensor in micrometers (µm).
sensor_length_um = 15.0

# The nominal fluid velocity from the syringe pump in meters per second (m/s).
fluid_velocity_ms = 0.2

# The maximum frequency the electronics can reliably process, in Hertz (Hz).
electronics_bandwidth_hz = 50000.0

# --- The Effect of the Improperly Positioned Magnet ---
# An improperly positioned magnet creates a field gradient, which exerts a force
# on the particle, adding to its velocity.
magnetic_velocity_increase_ms = 0.6 # m/s

# --- Calculations ---

# 1. Convert fluid velocity to µm/s for consistent units.
fluid_velocity_ums = fluid_velocity_ms * 1e6

print("--- Nominal Conditions (Correctly Positioned Magnet) ---")
# Calculate the time the particle takes to pass over the sensor.
# Transit Time = Sensor Length / Particle Velocity
transit_time_s = sensor_length_um / fluid_velocity_ums
# The characteristic frequency of the signal is related to the inverse of the transit time.
# A common approximation for the required bandwidth is 1 / (2 * pulse_width).
nominal_signal_freq_hz = 1.0 / (2 * transit_time_s)

print(f"Nominal particle velocity: {fluid_velocity_ms} m/s")
print(f"Equation: Transit Time = {sensor_length_um:.1f} µm / {fluid_velocity_ums:.0f} µm/s = {transit_time_s * 1e6:.2f} µs")
print(f"Equation: Estimated Signal Frequency = 1 / (2 * {transit_time_s:.8f} s) = {nominal_signal_freq_hz / 1e3:.2f} kHz")
print(f"This frequency ({nominal_signal_freq_hz / 1e3:.2f} kHz) is within the system bandwidth of {electronics_bandwidth_hz / 1e3:.0f} kHz.\n")


print("--- Aberrant Conditions (Improperly Positioned Magnet) ---")
# Calculate the new, higher velocity of the particle.
total_velocity_ms = fluid_velocity_ms + magnetic_velocity_increase_ms
total_velocity_ums = total_velocity_ms * 1e6

# Calculate the new, shorter transit time and higher frequency.
fast_transit_time_s = sensor_length_um / total_velocity_ums
fast_signal_freq_hz = 1.0 / (2 * fast_transit_time_s)

print(f"Magnetic force adds {magnetic_velocity_increase_ms} m/s to the particle's velocity.")
print(f"Equation: New Velocity = {fluid_velocity_ms:.1f} m/s + {magnetic_velocity_increase_ms:.1f} m/s = {total_velocity_ms:.1f} m/s")
print(f"Equation: New Transit Time = {sensor_length_um:.1f} µm / {total_velocity_ums:.0f} µm/s = {fast_transit_time_s * 1e6:.2f} µs")
print(f"Equation: New Signal Frequency = 1 / (2 * {fast_transit_time_s:.8f} s) = {fast_signal_freq_hz / 1e3:.2f} kHz")


# --- Conclusion ---
# Compare the new frequency to the system's bandwidth.
if fast_signal_freq_hz > electronics_bandwidth_hz:
    print(f"\nConclusion: The new signal frequency ({fast_signal_freq_hz / 1e3:.2f} kHz) EXCEEDS the system's bandwidth of {electronics_bandwidth_hz / 1e3:.0f} kHz.")
    print("The signal will be distorted or missed, leading to failed detection.")
else:
    print(f"\nConclusion: The new signal frequency ({fast_signal_freq_hz / 1e3:.2f} kHz) is still within the system's bandwidth of {electronics_bandwidth_hz / 1e3:.0f} kHz.")
    print("In this specific case, the system can still detect the particle.")
