import scipy.constants

def demonstrate_gas_kinetics():
    """
    This script explains why temperature is the crucial parameter for the
    Maxwell's demon thought experiment apparatus. It uses the formula for
    the average kinetic energy of gas particles.
    """
    # Introduction to the principle
    print("Principle: The movement of ideal gas particles is driven by their thermal energy.")
    print("This relationship is described by the Kinetic Theory of Gases.")
    print("The formula for average kinetic energy (KE_avg) is: KE_avg = (3/2) * k_B * T\n")
    
    # Define constants for the equation
    coefficient = 1.5  # This is the 3/2 part
    boltzmann_constant = scipy.constants.k
    
    print("Let's analyze two cases:\n")
    
    # --- Case 1: Absolute Zero Temperature ---
    temp_case_1 = 0.0  # in Kelvin
    ke_case_1 = coefficient * boltzmann_constant * temp_case_1
    
    print(f"Case 1: Temperature is absolute zero (T = {temp_case_1} K).")
    # Print the equation with its numbers
    print(f"Equation: KE_avg = {coefficient} * k_B * {temp_case_1}")
    print(f"Result: The average kinetic energy is {ke_case_1} Joules.")
    print("Interpretation: At T=0 K, particles have no motion. They cannot pass through the door.\n")
    
    # --- Case 2: Non-Zero Temperature ---
    temp_case_2 = 293.15  # Room temperature in Kelvin
    ke_case_2 = coefficient * boltzmann_constant * temp_case_2

    print(f"Case 2: Temperature is above absolute zero (e.g., T = {temp_case_2} K).")
    # Print the equation with its numbers
    print(f"Equation: KE_avg = {coefficient} * k_B * {temp_case_2}")
    print(f"Result: The average kinetic energy is approximately {ke_case_2:.2e} Joules.")
    print("Interpretation: At T>0 K, particles are in constant, random motion. This motion allows them to pass through the one-way door until all are trapped on one side.\n")

    # --- Final Conclusion ---
    print("Conclusion: The experimental parameter required for the gas to become trapped is a non-zero Temperature.")

# Execute the demonstration
demonstrate_gas_kinetics()