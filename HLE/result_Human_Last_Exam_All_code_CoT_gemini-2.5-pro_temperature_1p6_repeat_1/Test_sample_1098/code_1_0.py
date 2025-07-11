import math

def demonstrate_ratchet_temperature_dependence():
    """
    This function demonstrates the principle behind the Maxwell's demon thought experiment
    using a one-way door (Brownian ratchet) analogy.

    A physical one-way door can only successfully trap gas on one side if there is
    a temperature difference between the gas and the door itself. This is because the
    door's own thermal energy causes it to randomly fluctuate, sometimes allowing
    molecules to "leak" through in the wrong direction. This leakage probability
    is governed by the Boltzmann factor.

    This code calculates a "leakage factor" for different door temperatures to show that
    only as its temperature approaches absolute zero does the door become a perfect
    one-way valve.
    """

    # --- Constants and Parameters ---
    # Boltzmann constant in Joules/Kelvin
    k = 1.380649e-23
    # A hypothetical energy barrier the door's latch must overcome to fail (in Joules)
    E_barrier = 1.0e-20
    # A list of door temperatures to test (in Kelvin)
    door_temperatures = [300.0, 77.0, 4.0, 1.0, 0.1]

    print("This demonstration calculates the 'leakage factor' of a one-way door.")
    print("This factor represents the probability of failure due to the door's own thermal energy.")
    print("Equation: Leakage Factor = exp(-E_barrier / (k * T_door))")
    print("-" * 60)
    print(f"Assumed Energy Barrier (E_barrier): {E_barrier:.2e} J")
    print(f"Boltzmann Constant (k): {k:.2e} J/K\n")

    # --- Calculation for each temperature ---
    for T_door in door_temperatures:
        # To trap all the gas, this factor must be effectively zero.
        
        # Calculate the value of the exponent in the Boltzmann factor.
        # This is the final equation whose numbers will be printed.
        exponent_val = -E_barrier / (k * T_door)
        
        # Calculate the resulting leakage factor.
        leakage_factor = math.exp(exponent_val)

        print(f"--- Case: Door Temperature (T_door) = {T_door} K ---")
        # Print each number in the final equation as requested.
        print(f"Calculation: exp( -({E_barrier}) / (({k}) * ({T_door})) )")
        print(f"Result: exp({exponent_val:.4f})")
        print(f"Leakage Factor: {leakage_factor:.4e}\n")

    print("-" * 60)
    print("Conclusion:")
    print("If the door has the same temperature as the gas (e.g., ~300 K), the leakage is significant.")
    print("This allows the system to remain in equilibrium, with no net flow of gas to one side.")
    print("Only when the door's temperature approaches absolute zero (0 K) does the leakage factor")
    print("approach zero, allowing the door to perfectly rectify motion and trap the gas.")
    print("\nTherefore, Temperature is the key experimental parameter.")

# Execute the demonstration
demonstrate_ratchet_temperature_dependence()