import math

def simulate_liquid_crystal_design():
    """
    This script simulates the design of a single-ring liquid crystal molecule
    based on the rules provided. It aims to find the optimal alkyl chain length 'n'
    for a transition temperature near room temperature.
    """
    
    # --- Part 1: Define Key Requirements and Model (based on prompt) ---
    # General Structure: C(n)H(2n+1)-Ph-CN
    TARGET_TEMP = 22.0  # Target room temperature in °C
    TOLERANCE = 2.0     # How close we need to get to the target temperature

    # Start with the suggested molecule (n=5), inspired by 5CB.
    # We will assume a hypothetical transition temperature for this starting molecule.
    # Let's assume T_NI for n=5 is 35.0°C.
    initial_n = 5
    base_temp_for_n5 = 35.0
    
    # Per rule F, increasing chain length lowers the temperature. We'll model this
    # as a linear change. Let's assume a 5.0°C change per carbon atom in the chain.
    temp_change_per_carbon = 5.0

    def calculate_transition_temp(n):
        """
        Calculates the estimated transition temperature based on alkyl chain length 'n'
        using our defined model.
        Equation: T(n) = T_at_n5 - change_per_carbon * (n - n_initial)
        """
        return base_temp_for_n5 - temp_change_per_carbon * (n - initial_n)

    # --- Part 2: Iteratively Find the Best Design ---
    current_n = initial_n
    max_iterations = 10

    print(f"Goal: Design a C(n)H(2n+1)-Ph-CN molecule with a transition temp near {TARGET_TEMP}°C.")
    print(f"Starting with n={initial_n}, which has a model transition temperature of {base_temp_for_n5}°C.\n")

    for i in range(max_iterations):
        current_temp = calculate_transition_temp(current_n)
        error = current_temp - TARGET_TEMP
        
        print(f"--- Iteration {i + 1} ---")
        print(f"  Trying molecule: C{current_n}H{2*current_n+1}-Ph-CN")
        
        # Display the calculation for the current step
        print(f"  Calculation: {base_temp_for_n5} - {temp_change_per_carbon} * ({current_n} - {initial_n}) = {current_temp:.1f}°C")

        if abs(error) <= TOLERANCE:
            print(f"\nSUCCESS: The estimated temperature {current_temp:.1f}°C is within {TOLERANCE}°C of the target {TARGET_TEMP}°C.")
            break

        if error > 0:
            print(f"  Result: Temperature {current_temp:.1f}°C is too high. Increasing chain length 'n' to lower it.")
            current_n += 1
        else:
            print(f"  Result: Temperature {current_temp:.1f}°C is too low. Decreasing chain length 'n' to raise it.")
            current_n -= 1
        print("") # Newline for readability
    else:
         print("\nWarning: Maximum iterations reached without converging on the target temperature.")

    # --- Part 3: Final Answer ---
    final_temp = calculate_transition_temp(current_n)
    
    print("\n------------------------------")
    print("--- Final Proposed Design ---")
    print("------------------------------")
    print(f"Optimal chain length 'n': {current_n}")
    print(f"Final Proposed Molecule: C{current_n}H{2*current_n+1}-Ph-CN")
    print(f"Final Estimated Temperature: {final_temp:.1f}°C")
    print("\nFinal Equation:")
    # Here we output the numbers in the final equation as requested
    print(f"{final_temp:.1f} = {base_temp_for_n5} - {temp_change_per_carbon} * ({current_n} - {initial_n})")

simulate_liquid_crystal_design()
<<<8>>>