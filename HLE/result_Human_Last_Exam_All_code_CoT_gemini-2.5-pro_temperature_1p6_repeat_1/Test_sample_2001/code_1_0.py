import sys

def simulate_sensor_response():
    """
    This script simulates the effect of magnet positioning on a Spin Valve sensor.
    It calculates the total magnetic field experienced by the sensor and checks if it
    exceeds the sensor's saturation limit.
    """

    # --- System Parameters (values are illustrative, in Oersted units) ---

    # The small magnetic field from a single paramagnetic particle to be detected.
    H_particle = 0.5  # Oe

    # The Spin Valve sensor has a limited linear operating range.
    # Beyond this range, it saturates.
    H_sensor_max = 5.0  # Oe
    H_sensor_min = -5.0 # Oe

    # --- Scenario 1: Correctly Positioned Magnet ---
    # The stray magnetic field from a well-aligned magnet in the sensor's
    # sensitive direction is very small.
    H_magnet_ideal = 0.1 # Oe

    # --- Scenario 2: Improperly Positioned Magnet ---
    # A misaligned magnet creates a large stray field in the sensor's
    # sensitive direction, even if the vertical component is correct.
    H_magnet_misaligned = 10.0 # Oe

    print("--- Analyzing Magnet Positioning Effects on a Spin Valve Sensor ---")

    # --- Calculation for Ideal Case ---
    print("\n[SCENARIO 1: IDEAL POSITIONING]")
    H_total_ideal = H_magnet_ideal + H_particle
    is_saturated_ideal = not (H_sensor_min < H_total_ideal < H_sensor_max)

    print("The total field is the sum of the background field and the particle's field.")
    print("Equation: H_total = H_background + H_particle")
    # Outputting each number in the final equation
    print(f"Calculation: {H_total_ideal:.2f} Oe = {H_magnet_ideal:.2f} Oe + {H_particle:.2f} Oe")
    print(f"Sensor Operating Range: [{H_sensor_min:.2f} Oe, {H_sensor_max:.2f} Oe]")
    print(f"Result: The total field is within the operating range. Sensor is NOT saturated: {is_saturated_ideal}")


    # --- Calculation for Misaligned Case ---
    print("\n[SCENARIO 2: MISALIGNED POSITIONING]")
    H_total_misaligned = H_magnet_misaligned + H_particle
    is_saturated_misaligned = not (H_sensor_min < H_total_misaligned < H_sensor_max)

    print("The large background field from the misaligned magnet now dominates.")
    print("Equation: H_total = H_background + H_particle")
    # Outputting each number in the final equation
    print(f"Calculation: {H_total_misaligned:.2f} Oe = {H_magnet_misaligned:.2f} Oe + {H_particle:.2f} Oe")
    print(f"Sensor Operating Range: [{H_sensor_min:.2f} Oe, {H_sensor_max:.2f} Oe]")
    print(f"Result: The total field exceeds the operating range. Sensor IS saturated: {is_saturated_misaligned}")

    print("\n--- CONCLUSION ---")
    print("An improperly positioned magnet can create a large background field that drives the sensor")
    print("into saturation, making it impossible to detect the smaller signals from the particles.")
    print("This corresponds to 'Spin Valve Saturation'.")


if __name__ == '__main__':
    simulate_sensor_response()
    # Adding the final answer tag to stdout for parsing
    sys.stdout.write("<<<D>>>\n")
