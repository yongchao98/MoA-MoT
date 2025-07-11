import math

# Define the parameters of the magnetic flow cytometry setup
channel_height_um = 10.0  # Channel height in micrometers
sensor_length_um = 5.0    # Assumed effective length of the sensor in micrometers
avg_flow_velocity_um_s = 1000.0  # Average velocity of the fluid in micrometers per second

# --- Step 1: Model the fluid velocity profile (Poiseuille flow) ---
# In a channel, fluid flows fastest at the center and is slowest near the walls.
# For a simplified parabolic profile, the maximum velocity (at the center)
# is 1.5 times the average velocity.
v_max_um_s = 1.5 * avg_flow_velocity_um_s

def get_velocity_at_position(z_position_um, h_channel_um):
    """Calculates fluid velocity based on vertical position 'z' from the center."""
    # z_position_um is the distance from the center (-h/2 to +h/2)
    # The flow is parabolic: v(z) = v_max * (1 - (z / (h/2))^2)
    half_height = h_channel_um / 2.0
    velocity = v_max_um_s * (1 - (z_position_um / half_height)**2)
    return velocity

# --- Step 2: Calculate the "expected" signal frequency ---
# Assume with a correctly positioned magnet, gravity and magnetic force push
# particles to the bottom. A 1 um particle sits with its center 0.5 um from the wall.
# Position from center z = - (channel_height/2 - particle_radius) = -(5 - 0.5) = -4.5 um
expected_particle_pos_z = - (channel_height_um / 2.0 - 0.5)
expected_velocity = get_velocity_at_position(expected_particle_pos_z, channel_height_um)

# The signal frequency is proportional to velocity / sensor_length
expected_frequency_hz = expected_velocity / sensor_length_um

# --- Step 3: Define the electronics bandwidth based on the expected frequency ---
# Assume the electronics are designed for the expected frequency +/- 20%
bandwidth_margin = 0.20
min_freq_hz = expected_frequency_hz * (1 - bandwidth_margin)
max_freq_hz = expected_frequency_hz * (1 + bandwidth_margin)

# --- Step 4: Simulate the effect of an improperly positioned magnet ---
# A misaligned magnet can create gradients that lift the particle or alter its
# equilibrium position. Let's simulate a case where the force pushes the particle
# to flow at the channel's center (z = 0), where velocity is maximum.
improper_pos_particle_z = 0.0
new_velocity = get_velocity_at_position(improper_pos_particle_z, channel_height_um)
new_frequency_hz = new_velocity / sensor_length_um

# --- Step 5: Output the results and conclusion ---
print("--- System Analysis ---")
print(f"Channel Height: {channel_height_um} µm")
print(f"Average Flow Velocity: {avg_flow_velocity_um_s} µm/s")
print(f"Maximum (Center) Velocity: {v_max_um_s:.2f} µm/s\n")

print("--- Case 1: Correctly Positioned Magnet ---")
print(f"Expected particle vertical position (from center): {expected_particle_pos_z} µm")
print(f"Resulting particle velocity: {expected_velocity:.2f} µm/s")
print(f"Expected signal frequency: {expected_frequency_hz:.2f} Hz")
print(f"System electronics bandwidth: {min_freq_hz:.2f} Hz to {max_freq_hz:.2f} Hz\n")

print("--- Case 2: Improperly Positioned Magnet ---")
print("An improper magnet position can create new forces, altering the particle's path.")
print(f"Simulated new vertical position (from center): {improper_pos_particle_z} µm")
print(f"Resulting new particle velocity: {new_velocity:.2f} µm/s")
print(f"Resulting new signal frequency: {new_frequency_hz:.2f} Hz\n")

print("--- Conclusion ---")
if new_frequency_hz < min_freq_hz or new_frequency_hz > max_freq_hz:
    print(f"The new frequency ({new_frequency_hz:.2f} Hz) is OUTSIDE the system's operating bandwidth ({min_freq_hz:.2f}-{max_freq_hz:.2f} Hz).")
    print("This can lead to signal loss or distortion.")
else:
    print(f"The new frequency ({new_frequency_hz:.2f} Hz) is INSIDE the system's operating bandwidth.")
