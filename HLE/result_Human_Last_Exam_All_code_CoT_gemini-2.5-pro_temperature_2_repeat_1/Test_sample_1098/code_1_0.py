import scipy.constants as const

def demonstrate_thermal_motion(temperature_kelvin):
    """
    Explains why temperature is the key parameter by calculating
    the average kinetic energy (KE) of a gas particle.
    The formula is: KE_avg = (3/2) * k * T
    """
    # The Boltzmann constant, k, relates temperature to energy.
    k = const.k

    # Calculate the average kinetic energy for a single particle.
    avg_kinetic_energy = 1.5 * k * temperature_kelvin

    print(f"--- Analyzing for Temperature = {temperature_kelvin} K ---")

    if temperature_kelvin <= 0:
        print("At absolute zero, particles have no kinetic energy and are motionless.")
        print("Without motion, no particle can cross the door.")
    else:
        print("Above absolute zero, particles have kinetic energy and are in constant motion.")
        print("This motion is required for particles to eventually pass through the one-way door.")

    # Per the instructions, outputting each number in the final equation.
    print("The final equation is: KE_avg = 1.5 * k * T")
    print("Calculated values:")
    print(f"{avg_kinetic_energy:.4g} (J) = 1.5 * {k:.4g} (J/K) * {temperature_kelvin} (K)")
    print("-" * 55)

print("The apparatus requires the gas molecules to be in motion to work.")
print("The motion of gas molecules is a direct result of their temperature.")
print("\nLet's calculate the kinetic energy at different temperatures:")

# Case 1: At absolute zero, no motion occurs.
demonstrate_thermal_motion(0)

# Case 2: At room temperature, motion occurs, allowing the process to proceed.
demonstrate_thermal_motion(298.15)
