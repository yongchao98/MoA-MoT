def check_condition_for_gas_accumulation(temperature_kelvin):
    """
    This function demonstrates the relationship between temperature and the
    molecular motion required for the one-way door mechanism to work.
    """
    # The Boltzmann constant (Joules/Kelvin)
    k = 1.380649e-23

    # The formula for the average kinetic energy (KE) of a molecule in an ideal gas is:
    # KE = (3/2) * k * T
    # where k is the Boltzmann constant and T is the temperature in Kelvin.

    # Calculate the average kinetic energy
    average_kinetic_energy = (3/2) * k * temperature_kelvin

    print(f"Testing for temperature: {temperature_kelvin} K")
    print("Equation: Average KE = (3/2) * k * T")
    print(f"Calculation: Average KE = (3/2) * ({k:.6e}) * ({temperature_kelvin}) = {average_kinetic_energy:.6e} Joules")

    # The one-way door requires molecules to be in motion to pass through it.
    # Molecular motion only exists if there is kinetic energy.
    if average_kinetic_energy > 0:
        print("Result: Molecules have kinetic energy and are in motion. The process will occur, and gas will accumulate on one side.")
    else:
        print("Result: Molecules have zero kinetic energy and are stationary. The process cannot occur.")
    print("-" * 50)


# Scenario 1: A typical non-zero temperature (e.g., Room Temperature)
room_temperature = 293.15
check_condition_for_gas_accumulation(room_temperature)

# Scenario 2: The absolute minimum temperature
absolute_zero = 0
check_condition_for_gas_accumulation(absolute_zero)
