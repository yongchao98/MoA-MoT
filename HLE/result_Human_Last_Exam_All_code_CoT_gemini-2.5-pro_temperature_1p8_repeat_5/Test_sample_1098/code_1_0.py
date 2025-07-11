import math

def calculate_rms_velocity(temperature, gas_name, molar_mass):
    """
    Calculates the root-mean-square velocity of a gas at a given temperature.
    Prints the step-by-step equation and the result.
    """
    # Constants
    R = 8.314  # Ideal gas constant in J/(mol·K)
    M = molar_mass  # Molar mass in kg/mol

    print(f"Calculating the velocity of {gas_name} molecules at {temperature} K...")

    # Ensure we don't divide by zero, although the formula handles T=0 gracefully.
    if M == 0:
        print("Molar mass cannot be zero.")
        return

    # Calculation
    # The term inside the square root is (3 * R * T) / M
    term = (3 * R * temperature) / M
    # RMS velocity is the square root of the term
    v_rms = math.sqrt(term) if term >= 0 else 0

    # Output the full equation with numbers
    print(f"v_rms = sqrt(3 * {R} * {temperature} / {M})")
    print(f"Result: The average molecular velocity is {v_rms:.2f} m/s.\n")
    return v_rms

def main():
    """
    Main function to demonstrate the role of temperature.
    """
    print(
        "For the 'Maxwell's demon' apparatus to work, gas molecules must be in motion.\n"
        "This motion is a direct result of the gas's temperature. Let's see how\n"
        "molecular velocity changes with temperature for Nitrogen gas (N2).\n"
    )

    gas_name_N2 = "Nitrogen (N2)"
    molar_mass_N2 = 0.028014  # Molar mass of N2 in kg/mol

    # Case 1: Room Temperature (T > 0 K)
    temp_room = 298.15  # approx. 25°C or 77°F
    calculate_rms_velocity(temp_room, gas_name_N2, molar_mass_N2)

    # Case 2: Absolute Zero (T = 0 K)
    temp_zero = 0.0
    velocity_at_zero = calculate_rms_velocity(temp_zero, gas_name_N2, molar_mass_N2)

    if velocity_at_zero == 0:
        print("At absolute zero (0 K), molecules stop moving. Therefore, the gas cannot\n"
              "pass through the door. This shows that a non-zero TEMPERATURE is the\n"
              "required parameter for the experiment to function.")

if __name__ == "__main__":
    main()