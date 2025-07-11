def design_liquid_crystal():
    """
    Suggests a liquid crystal molecule based on a simple predictive model
    derived from the user's design rules.
    """
    # 1. Key requirements and assumptions from the prompt
    TARGET_T = 22.5  # Target room temperature in degrees Celsius
    
    # General Structure: C(n)H(2n+1)-Ph-CN (Alkyl-Cyanobenzene)
    
    # Model based on Rule F:
    # "Start with a pentyl chain (n=5)"
    BASE_N = 5
    
    # A single-ring molecule has a lower transition temp than a multi-ring one.
    # We'll assume a hypothetical base temperature for our n=5 molecule.
    BASE_T = 10.0  # Hypothetical transition temp for C5H11-Ph-CN
    
    # "If transition temperature is too high, increase chain length"
    # "If too low, decrease chain length"
    # This implies temperature increases with chain length (for this model).
    # Let's assume a 4.0°C increase per carbon atom added to the chain.
    T_CHANGE_PER_N = 4.0

    print("--- Liquid Crystal Design Assistant ---")
    print(f"Goal: Design a C(n)H(2n+1)-Ph-CN molecule with a transition temperature near {TARGET_T}°C.")
    print(f"Using a model based on your rules (A-F).\n")

    # 2. Calculation to find the optimal chain length 'n'
    # We want to solve for n in the equation: Target_T = Base_T + (n - Base_n) * T_change
    # n = (Target_T - Base_T) / T_change + Base_n
    
    # Calculate the ideal 'n' as a floating-point number first
    n_float = (TARGET_T - BASE_T) / T_CHANGE_PER_N + BASE_N
    
    # Since 'n' must be an integer, we round to the nearest whole number
    n_ideal = round(n_float)
    if n_ideal < 1: # Ensure chain length is at least 1
        n_ideal = 1

    # 3. Predict the transition temperature for the suggested molecule
    predicted_t = BASE_T + (n_ideal - BASE_N) * T_CHANGE_PER_N
    hydrogen_count = 2 * n_ideal + 1
    molecular_formula = f"C{n_ideal}H{hydrogen_count}-Ph-CN"

    print("--- Suggested Molecular Design ---")
    print(f"To achieve the target temperature, the suggested alkyl chain length is: n = {n_ideal}")
    print(f"Suggested Molecule Structure: {molecular_formula}")
    
    print("\n--- Predicted Transition Temperature Calculation ---")
    # As requested, printing the final equation with all numbers:
    print(f"Predicted T_ni = Base_T + (n - Base_n) * T_change_per_n")
    print(f"Predicted T_ni = {BASE_T}°C + ({n_ideal} - {BASE_N}) * {T_CHANGE_PER_N}°C = {predicted_t:.1f}°C")

# Execute the design function
if __name__ == '__main__':
    design_liquid_crystal()