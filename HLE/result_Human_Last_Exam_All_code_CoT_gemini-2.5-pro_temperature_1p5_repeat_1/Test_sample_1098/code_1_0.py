import math

def analyze_molecular_motion():
    """
    Analyzes the effect of temperature on the kinetic energy of gas molecules.
    """
    # The Boltzmann constant (Joules per Kelvin)
    BOLTZMANN_CONSTANT = 1.380649e-23

    # The factor from the equipartition theorem for a monoatomic ideal gas
    FACTOR = 3/2

    # --- Case 1: Absolute Zero ---
    temp_absolute_zero = 0  # in Kelvin

    # Calculate kinetic energy at absolute zero
    ke_zero = FACTOR * BOLTZMANN_CONSTANT * temp_absolute_zero

    print("This experiment relies on the random motion of gas molecules to pass through a one-way door.")
    print("The motion of molecules is determined by their kinetic energy (KE).")
    print("Let's analyze the effect of Temperature (T) on KE using the formula: KE = (3/2) * k * T\n")

    print(f"--- Analyzing at Absolute Zero Temperature ---")
    print(f"The equation is: KE = {FACTOR} * {BOLTZMANN_CONSTANT:.4e} J/K * T")
    print(f"Substituting T = {temp_absolute_zero} K:")
    # Showing each number in the final equation for calculation
    print(f"KE = {FACTOR} * {BOLTZMANN_CONSTANT:.4e} * {temp_absolute_zero}")
    print(f"Resulting Kinetic Energy = {ke_zero:.4e} Joules")
    print("At absolute zero (0 K), molecules have no kinetic energy. Without energy, they do not move, so they cannot pass through the door.\n")

    # --- Case 2: Non-Zero Temperature ---
    temp_non_zero = 298  # Room temperature in Kelvin

    # Calculate kinetic energy at a non-zero temperature
    ke_non_zero = FACTOR * BOLTZMANN_CONSTANT * temp_non_zero

    print(f"--- Analyzing at a Non-Zero Temperature (Room Temp) ---")
    print(f"The equation is: KE = {FACTOR} * {BOLTZMANN_CONSTANT:.4e} J/K * T")
    print(f"Substituting T = {temp_non_zero} K:")
    # Showing each number in the final equation for calculation
    print(f"KE = {FACTOR} * {BOLTZMANN_CONSTANT:.4e} * {temp_non_zero}")
    print(f"Resulting Kinetic Energy = {ke_non_zero:.4e} Joules")
    print("At a non-zero temperature, molecules have kinetic energy and are in constant random motion.")
    print("This motion will eventually cause them to pass through the door, accumulating on one side.\n")

    print("Conclusion: A non-zero Temperature is the essential experimental parameter required for the gas molecules to have motion and for the experiment to proceed.")

analyze_molecular_motion()