def solve_maxwells_demon():
    """
    Analyzes the Maxwell's demon thought experiment with a Brownian ratchet
    to determine the key experimental parameter.
    """
    print("Analyzing the thought experiment:")
    print("1. The apparatus is designed with a one-way door (a Brownian ratchet) separating two chambers of gas.")
    print("2. Naively, random motion would suggest all gas molecules eventually pass through the one-way door, getting trapped on one side.")
    print("3. However, this would decrease the entropy of the gas, seemingly violating the Second Law of Thermodynamics.")
    print("\nThe resolution, as explained by Richard Feynman, lies in the thermal energy of the system components themselves.")
    print("\nCondition for Failure (Equilibrium):")
    print("If the entire system (gas and the one-way door) is at a single, uniform TEMPERATURE, the door itself will experience random thermal motion (Brownian motion). This thermal jiggling will cause the door to open randomly, allowing gas molecules to pass back. At thermal equilibrium, the flow of gas in both directions will be equal, and no net transfer will occur.")
    
    print("\nCondition for Success (Non-Equilibrium):")
    print("For the door to successfully act as a one-way gate and trap all the gas, it must be at a different TEMPERATURE than the gas (specifically, it must be colder).")
    print("This TEMPERATURE difference means the system is not in equilibrium. The flow of heat from the hotter gas to the colder door provides the energy to operate the ratchet, allowing it to rectify the random motion of the gas molecules.")

    print("\nConclusion:")
    print("The essential experimental parameter required to create the necessary non-equilibrium condition is TEMPERATURE.")

solve_maxwells_demon()