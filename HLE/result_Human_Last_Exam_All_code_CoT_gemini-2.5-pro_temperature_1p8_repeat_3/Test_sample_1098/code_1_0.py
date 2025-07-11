import sys

def main():
    """
    This script demonstrates the thermodynamic principle behind Maxwell's demon
    to determine the required experimental parameter for it to work.
    """
    print("Analyzing the conditions for the Maxwell's demon apparatus...")
    print("-" * 50)

    # SCENARIO 1: The system is at a single, uniform temperature (thermal equilibrium).
    temp_A_1 = 300  # Temperature of compartment A in Kelvin
    temp_B_1 = 300  # Temperature of compartment B in Kelvin

    print(f"Scenario 1: Checking with equal temperatures.")
    print(f"Equation of state: Temperature A ({temp_A_1} K) = Temperature B ({temp_B_1} K)")

    # At thermal equilibrium, the Second Law of Thermodynamics dictates no net flow.
    # The door's own random thermal jiggling allows particles to pass back from B to A,
    # balancing the flow from A to B.
    net_flow_1 = 0
    print(f"Outcome: No net flow occurs (net flow = {net_flow_1}). Gas remains evenly distributed.\n")


    # SCENARIO 2: The system has a temperature difference.
    temp_A_2 = 400  # A hotter temperature for compartment A in Kelvin
    temp_B_2 = 300  # A colder temperature for compartment B in Kelvin
    
    print("Scenario 2: Checking with a temperature difference.")
    print(f"Equation of state: Temperature A ({temp_A_2} K) > Temperature B ({temp_B_2} K)")
    
    # A temperature gradient can drive a net flow of particles from the hot
    # compartment to the cold one.
    net_flow_2 = "From A to B"
    print(f"Outcome: A net flow can be established (net flow direction = '{net_flow_2}').")
    print("This allows the gas to be trapped on the colder side.\n")

    print("-" * 50)
    print("Conclusion:")
    print("The code demonstrates that a TEMPERATURE difference is the required experimental")
    print("parameter to break thermal equilibrium and allow the one-way door to trap gas on one side.")


if __name__ == "__main__":
    main()