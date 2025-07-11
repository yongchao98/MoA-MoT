import math

def demonstrate_temperature_effect():
    """
    Demonstrates that molecular motion depends on temperature.
    """
    # Constants
    BOLTZMANN_CONSTANT = 1.380649e-23  # J/K
    # Mass of a single Nitrogen molecule (N2) in kg
    # Molar mass of N2 is approx 28 g/mol or 0.028 kg/mol
    # Avogadro's number is approx 6.022e23 /mol
    MASS_OF_N2_MOLECULE = 0.028 / 6.022e23 # kg

    def calculate_rms_velocity(temperature_kelvin, mass_kg):
        """Calculates the root-mean-square velocity of a gas particle."""
        if temperature_kelvin < 0:
            print("Temperature cannot be negative in Kelvin.")
            return 0
        # The formula for RMS velocity is v = sqrt(3 * k * T / m)
        # where k is Boltzmann constant, T is temperature, m is mass
        v_rms_squared = (3 * BOLTZMANN_CONSTANT * temperature_kelvin) / mass_kg
        return math.sqrt(v_rms_squared)

    print("Analyzing the effect of temperature on molecular motion for the Maxwell's demon apparatus.")
    print("-" * 80)

    # Case 1: Absolute Zero Temperature
    temp_zero = 0  # Kelvin
    velocity_zero = calculate_rms_velocity(temp_zero, MASS_OF_N2_MOLECULE)
    print(f"At a temperature of T = {temp_zero} K (absolute zero):")
    print(f"The root-mean-square velocity of gas molecules is {velocity_zero:.2f} m/s.")
    print("Conclusion: Without temperature, there is no molecular motion. The gas cannot move between chambers.\n")

    # Case 2: Room Temperature
    temp_room = 298.15  # Kelvin (approx 25°C)
    velocity_room = calculate_rms_velocity(temp_room, MASS_OF_N2_MOLECULE)
    print(f"At a room temperature of T = {temp_room:.2f} K:")
    print(f"The root-mean-square velocity of gas molecules is {velocity_room:.2f} m/s.")
    print("Conclusion: At non-zero temperatures, molecules are in constant motion, allowing them to pass through the one-way door.")
    print("-" * 80)
    print("Therefore, Temperature is the essential experimental parameter required for the process to occur.")

demonstrate_temperature_effect()