import numpy as np

def simulate_sv_sensor_response():
    """
    This simulation models the effect of the permanent magnet's position on the
    Spin Valve (SV) sensor's ability to detect a paramagnetic particle.
    """

    # --- Define System Parameters (in arbitrary units) ---

    # SV Sensor Characteristics
    # The sensor has a linear response to a magnetic field within a certain range.
    # Outside this range, it saturates.
    sv_linear_range_max = 5.0  # Max field for linear response (e.g., in mT)
    sv_sensitivity = 10.0      # Output signal per unit of field (e.g., in V/mT)
    sv_saturation_output = sv_linear_range_max * sv_sensitivity # Max output signal

    # Particle's Magnetic Field
    # The small magnetic field produced by a magnetized particle at the sensor.
    particle_field = 0.5  # (e.g., in mT)

    def get_sv_output(total_field):
        """Calculates the sensor's output signal based on the total magnetic field."""
        if abs(total_field) >= sv_linear_range_max:
            # Sensor is saturated
            return np.sign(total_field) * sv_saturation_output
        else:
            # Sensor is in the linear range
            return total_field * sv_sensitivity

    # --- Scenario 1: Properly Positioned Magnet ---
    # The bias field is set to place the sensor in its most sensitive (linear) region.
    # A bias of 0.0 is optimal for detecting positive and negative field changes.
    bias_field_proper = 0.0

    print("--- Scenario 1: Magnet is Properly Positioned ---")
    # Calculate the sensor's output signal from the bias field alone
    output_no_particle_proper = get_sv_output(bias_field_proper)
    print(f"Equation (no particle): Output = get_sv_output({bias_field_proper:.1f})")
    print(f"Sensor output (no particle): {output_no_particle_proper:.2f}")

    # A particle passes by, adding its field to the bias field
    total_field_proper = bias_field_proper + particle_field
    output_with_particle_proper = get_sv_output(total_field_proper)
    print(f"Equation (with particle): Output = get_sv_output({bias_field_proper:.1f} + {particle_field:.1f})")
    print(f"Sensor output (with particle): {output_with_particle_proper:.2f}")

    # The detected signal is the change in output
    detected_signal_proper = abs(output_with_particle_proper - output_no_particle_proper)
    print(f"Final Equation: Detected Signal = |{output_with_particle_proper:.2f} - {output_no_particle_proper:.2f}|")
    print(f"Result: Detected Signal = {detected_signal_proper:.2f}\n")


    # --- Scenario 2: Improperly Positioned Magnet (Too Close) ---
    # The bias field is too strong, pushing the sensor into saturation.
    bias_field_improper = 8.0  # This field is > sv_linear_range_max

    print("--- Scenario 2: Magnet is Improperly Positioned (Too Close) ---")
    # Calculate the sensor's output signal from the bias field alone
    output_no_particle_improper = get_sv_output(bias_field_improper)
    print(f"Equation (no particle): Output = get_sv_output({bias_field_improper:.1f})")
    print(f"Sensor output (no particle): {output_no_particle_improper:.2f} (Status: Saturated)")

    # A particle passes by, adding its field to the already strong bias field
    total_field_improper = bias_field_improper + particle_field
    output_with_particle_improper = get_sv_output(total_field_improper)
    print(f"Equation (with particle): Output = get_sv_output({bias_field_improper:.1f} + {particle_field:.1f})")
    print(f"Sensor output (with particle): {output_with_particle_improper:.2f} (Status: Saturated)")

    # The detected signal is the change in output
    detected_signal_improper = abs(output_with_particle_improper - output_no_particle_improper)
    print(f"Final Equation: Detected Signal = |{output_with_particle_improper:.2f} - {output_no_particle_improper:.2f}|")
    print(f"Result: Detected Signal = {detected_signal_improper:.2f}\n")

    print("Conclusion: When the magnet is too close, the strong bias field saturates the Spin Valve sensor.")
    print("The sensor's output is already at its maximum and does not change when a particle passes by.")
    print("This makes particle detection impossible, demonstrating the negative effect of 'Spin Valve Saturation'.")

if __name__ == '__main__':
    simulate_sv_sensor_response()
<<<D>>>