def solve():
    """
    This script simulates the effect of a permanent magnet's position
    on a Spin Valve (SV) sensor in a magnetic flow cytometry setup.
    """

    # --- Define a model for the Spin Valve (SV) sensor ---
    # The sensor has a linear response to a magnetic field up to a certain point,
    # after which it saturates and the output no longer increases.
    def get_sensor_output(field, saturation_field, sensitivity):
        """
        Calculates the sensor's output voltage based on the applied field.
        """
        if field >= saturation_field:
            # If the field is at or above saturation, the output is maxed out.
            return saturation_field * sensitivity
        else:
            # Otherwise, the output is proportional to the field.
            return field * sensitivity

    # --- System Parameters (in arbitrary units) ---
    # The magnetic field produced by a single paramagnetic particle as it passes the sensor.
    particle_signal_field = 2.0

    # The maximum field the Spin Valve sensor can measure before it saturates.
    sensor_saturation_field = 100.0

    # The sensitivity of the sensor (e.g., Volts per unit of magnetic field).
    sensor_sensitivity = 0.05

    print("--- System Simulation: Detecting a Paramagnetic Particle ---")
    print(f"Spin Valve sensor saturation field: {sensor_saturation_field}")
    print(f"Magnetic field from a single particle: {particle_signal_field}\n")

    # --- Scenario 1: Correctly Positioned Magnet ---
    # The biasing field is strong enough to magnetize particles but well within the sensor's operating range.
    correct_bias_field = 50.0

    print("--- Scenario 1: Correct Magnet Position ---")
    print(f"Biasing field from magnet: {correct_bias_field}")

    # Calculate the sensor's baseline output from the bias field alone.
    output_baseline_1 = get_sensor_output(correct_bias_field, sensor_saturation_field, sensor_sensitivity)

    # Calculate the total field when a particle is present.
    total_field_with_particle_1 = correct_bias_field + particle_signal_field
    output_with_particle_1 = get_sensor_output(total_field_with_particle_1, sensor_saturation_field, sensor_sensitivity)

    # The detected signal is the *change* in output voltage.
    detected_signal_1 = output_with_particle_1 - output_baseline_1

    print(f"Total field at sensor: {correct_bias_field} (bias) + {particle_signal_field} (particle) = {total_field_with_particle_1}")
    print(f"Sensor output change (detected signal): {detected_signal_1:.4f} V")
    print("Result: A clear signal is detected.\n")


    # --- Scenario 2: Improperly Positioned Magnet (Too Close) ---
    # The magnet is too close, creating a bias field that exceeds the sensor's saturation point.
    incorrect_bias_field = 110.0

    print("--- Scenario 2: Improper Magnet Position (Too Close) ---")
    print(f"Biasing field from magnet: {incorrect_bias_field}")

    # Calculate the sensor's baseline output from the strong bias field.
    # The sensor is already saturated by the bias field alone.
    output_baseline_2 = get_sensor_output(incorrect_bias_field, sensor_saturation_field, sensor_sensitivity)

    # Calculate the total field when a particle is present.
    total_field_with_particle_2 = incorrect_bias_field + particle_signal_field
    output_with_particle_2 = get_sensor_output(total_field_with_particle_2, sensor_saturation_field, sensor_sensitivity)

    # The detected signal is the *change* in output voltage.
    detected_signal_2 = output_with_particle_2 - output_baseline_2

    print(f"Total field at sensor: {incorrect_bias_field} (bias) + {particle_signal_field} (particle) = {total_field_with_particle_2}")
    print(f"Sensor output change (detected signal): {detected_signal_2:.4f} V")
    print("Result: No signal is detected because the sensor was already saturated by the magnet's bias field.\n")

    print("Conclusion: An improperly positioned magnet creating too strong a field leads to Spin Valve Saturation, preventing particle detection.")

solve()