def simulate_sensor_response(biasing_field, particle_field=0.1):
    """
    Simulates the response of a Spin Valve sensor to a particle's magnetic field,
    under the influence of a larger biasing field from a permanent magnet.

    Args:
        biasing_field (float): The strength of the field from the permanent magnet.
        particle_field (float): The small magnetic field from the particle to be detected.
    """
    # These values represent the sensor's characteristics
    optimal_bias_min = 5.0
    optimal_bias_max = 10.0
    saturation_threshold = 15.0

    print(f"--- Testing with Biasing Field: {biasing_field:.1f} ---")

    # Determine the sensor's state and sensitivity based on the biasing field
    if biasing_field > saturation_threshold:
        sensitivity = 0.0  # Sensor is saturated, no sensitivity
        state = "SATURATED"
    elif optimal_bias_min <= biasing_field <= optimal_bias_max:
        sensitivity = 100.0  # High sensitivity in the optimal linear region
        state = "OPTIMAL"
    else:
        sensitivity = 10.0  # Low sensitivity outside the optimal range
        state = "NON-OPTIMAL"
    
    print(f"Sensor State: {state}. Biasing field is in the {state.lower()} range.")

    # The detected signal is the particle's field multiplied by the sensor's sensitivity
    detected_signal = particle_field * sensitivity
    
    print(f"The magnetic field from the particle is {particle_field}.")
    print(f"The sensor's sensitivity in this state is {sensitivity}.")
    print(f"Resulting Detected Signal: {detected_signal:.2f}\n")
    return detected_signal

def main():
    """
    Main function to demonstrate the effect of magnet positioning on sensor performance.
    """
    print("This simulation shows how an improperly positioned magnet can cause Spin Valve saturation.\n")

    # Case 1: Magnet is correctly positioned, providing an optimal biasing field.
    print("Case 1: Correctly Positioned Magnet")
    correct_biasing_field = 8.0  # An optimal field strength
    simulate_sensor_response(correct_biasing_field)

    # Case 2: Magnet is improperly positioned (too close), providing a very strong biasing field.
    print("Case 2: Improperly Positioned Magnet (Too Close)")
    improper_biasing_field = 20.0  # A field strong enough to cause saturation
    simulate_sensor_response(improper_biasing_field)

    print("Conclusion:")
    print("As the simulation shows, when the biasing field from the magnet is too high (Case 2),")
    print("the sensor enters saturation. Its sensitivity drops to zero, and it can no longer produce")
    print("a signal from the passing particle. This means the system fails to detect particles,")
    print("which is a significant negative effect directly caused by the magnet's improper position.")

if __name__ == "__main__":
    main()
<<<D>>>