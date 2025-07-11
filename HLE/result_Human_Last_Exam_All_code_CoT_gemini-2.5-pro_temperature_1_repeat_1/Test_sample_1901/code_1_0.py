import sys

def design_liquid_crystal():
    """
    Designs a liquid crystal molecule based on a set of rules.

    The script aims to find the optimal alkyl chain length 'n' for a molecule
    of the family CnH2n+1-Ph-CN to have a transition temperature
    near room temperature (20-25°C).
    """

    # --- 1. Define Design Parameters ---
    TARGET_TEMP_MIN = 20.0
    TARGET_TEMP_MAX = 25.0
    
    # Starting parameters based on the problem description
    # Start with a pentyl chain (n=5)
    current_n = 5
    
    # Model parameters for T(n) = base_temp - slope * (n - base_n)
    # Based on the idea that 5CB's single-ring analogue is too high
    # and that increasing n decreases the temperature.
    BASE_TEMP = 35.0
    SLOPE = 5.0
    BASE_N = 5

    # Helper for naming the final molecule
    alkyl_names = {
        1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl",
        5: "pentyl", 6: "hexyl", 7: "heptyl", 8: "octyl",
        9: "nonyl", 10: "decyl"
    }
    
    print("--- Liquid Crystal Design Process ---")
    print(f"Goal: Find 'n' for CnH2n+1-Ph-CN with transition temp in [{TARGET_TEMP_MIN}, {TARGET_TEMP_MAX}]°C.\n")
    print("Using model: T(n) = 35.0 - 5.0 * (n - 5)")
    print("-------------------------------------\n")

    # --- 2. Iteratively Search for the Best Molecule ---
    # Set a limit to prevent infinite loops
    max_iterations = 10
    for i in range(max_iterations):
        # Calculate transition temperature using the model
        calculated_temp = BASE_TEMP - SLOPE * (current_n - BASE_N)

        print(f"Iteration {i+1}: Trying n = {current_n}")
        print(f"   Calculated Temperature = {calculated_temp:.1f}°C")

        # Check if the temperature is in the target range
        if TARGET_TEMP_MIN <= calculated_temp <= TARGET_TEMP_MAX:
            print("\n--- Design Found! ---")
            alkyl_name = alkyl_names.get(current_n, f'{current_n}-carbon alkyl')
            molecule_name = f"4-{alkyl_name}benzonitrile"
            molecule_formula = f"C{current_n}H{2*current_n+1}-Ph-CN"
            
            print(f"Suggested Molecule: {molecule_name}")
            print(f"Formula: {molecule_formula}")
            
            # Print the final equation with all numbers, as requested
            print("\nFinal Equation:")
            # We use repr() to ensure floats are printed with the decimal, e.g., 35.0
            print(f"Transition Temperature = {repr(BASE_TEMP)} - {repr(SLOPE)} * ({current_n} - {BASE_N}) = {calculated_temp:.1f} °C")
            return
            
        # Adjust 'n' based on the rules
        elif calculated_temp > TARGET_TEMP_MAX:
            print(f"   Result: Temperature is too high. Increasing chain length (n).\n")
            current_n += 1
        else: # calculated_temp < TARGET_TEMP_MIN
            print(f"   Result: Temperature is too low. Decreasing chain length (n).\n")
            current_n -= 1
        
        # Safety check for chain length
        if not (1 < current_n < 15):
             print("\n--- Design Failed ---")
             print("Search resulted in an unrealistic chain length. The model may need refinement.")
             break

    if i == max_iterations -1:
        print("\n--- Design Failed ---")
        print("Could not find a suitable molecule within the iteration limit.")

# Execute the design function
design_liquid_crystal()