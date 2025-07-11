import math

def design_liquid_crystal():
    """
    Designs a liquid crystal by finding the optimal alkyl chain length 'n'
    to achieve a target transition temperature.
    """
    # A. Key requirements are defined.
    target_temperature = 25  # degrees Celsius

    # E. The general structure is CnH2n+1-Ph-CN.
    # F. We start with n=5 and tune from there.
    # To make this problem solvable, we will create a predictive model based on the rules.
    # Let's assume for n=5, the transition temperature is 40°C.
    # Let's assume each additional CH2 group decreases the temperature by 5°C.
    # This creates the model: T(n) = 40 - 5 * (n - 5)
    
    n_initial = 5
    temp_initial = 40
    temp_change_per_n = 5

    print("Liquid Crystal Design Simulation")
    print("------------------------------")
    print(f"Goal: Achieve a transition temperature of {target_temperature}°C.")
    print(f"Molecular Template: C(n)H(2n+1)-Ph-CN")
    print(f"Model: T(n) = {temp_initial} - {temp_change_per_n} * (n - {n_initial})")
    print("------------------------------\n")

    # Start the search for the optimal 'n'
    current_n = n_initial
    
    # Using a simple algebraic solution for efficiency
    # target = initial_temp - change * (n - initial_n)
    # (target - initial_temp) / -change = n - initial_n
    # n = (target - initial_temp) / -change + initial_n
    
    optimal_n_float = (target_temperature - temp_initial) / -temp_change_per_n + n_initial
    # Since n must be an integer, we round it.
    optimal_n = round(optimal_n_float)

    print(f"Solving for n where T(n) = {target_temperature}:")
    print(f"{target_temperature} = {temp_initial} - {temp_change_per_n} * (n - {n_initial})")
    print(f"{target_temperature - temp_initial} = -{temp_change_per_n} * (n - {n_initial})")
    print(f"({target_temperature - temp_initial}) / (-{temp_change_per_n}) = n - {n_initial}")
    print(f"{(target_temperature - temp_initial) / -temp_change_per_n} = n - {n_initial}")
    print(f"{(target_temperature - temp_initial) / -temp_change_per_n} + {n_initial} = n")
    print(f"n = {optimal_n_float}\n")


    print("--- Final Design ---")
    print(f"The optimal integer chain length is n = {optimal_n}.")

    # D. Final specific structure
    final_molecule_formula = f"C{optimal_n}H{2*optimal_n+1}-Ph-CN"
    print(f"Proposed Molecular Structure: {final_molecule_formula}")

    # Verify the result with our model equation
    final_temperature = temp_initial - temp_change_per_n * (optimal_n - n_initial)
    print("\nVerification using the final equation:")
    # Print each number in the final equation as requested
    print(f"T({optimal_n}) = {temp_initial} - {temp_change_per_n} * ({optimal_n} - {n_initial})")
    print(f"Calculated Temperature = {final_temperature}°C")

design_liquid_crystal()